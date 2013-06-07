#!/usr/bin/env python

''' extract a frame and run a simulation '''

import re, os, sys, glob, getopt, shutil, subprocess, math, random

verbose = 0
nsteps_warmup = 10000
nsteps_enemin = 10000



def run(prjname, id, tm):
  srcdir, gro, top, mdp = extract_frame(prjname, id, tm)
  run_simul(prjname, srcdir, gro, top, mdp)



def extract_frame(prjname, id, tm0):
  ''' extract a pdb file from id and tm0 '''
  if not id.startswith("data"): id = "data" + id
  print "extracting %s %s" % (id, tm0)

  mdp = glob.glob("md*.mdp")[0]
  trrdt = get_trrdt(mdp)

  tpr = glob.glob("*.tpr")[0]
  trr = glob.glob(id+"/*.trr")[0]
  gro = prjname + ".gro"
  tm = round(tm0 / trrdt) * trrdt
  print "time %g --> %g" % (tm0, tm)
  t0 = tm - trrdt*.5
  t1 = tm + trrdt*.5
  cmd = 'echo "0" | trjconv -s %s -f %s -o %s -b %s -e %s' % (
      tpr, trr, gro, t0, t1)
  runcmd(cmd, system = 1)

  curdir = os.path.abspath(os.curdir)
  top = "topol.top"
  print "file %s, %s, %s, %s" % (curdir, gro, top, mdp)
  return curdir, gro, "topol.top", mdp



def get_trrdt(mdp):
  ''' determine trr interval '''
  dt = 2e-3  # a guess
  xout = 20000
  for s in open(mdp).readlines():
    s = s.strip()
    if s.startswith("nstxout"):
      i = s.find("=")
      xout = int( s[i+1:].strip() )
    if s.startswith("dt"):
      i = s.find("=")
      dt = float( s[i+1:].strip() )
  print "%s: dt = %g, nstxout = %d, trrdt = %g" % (mdp, dt, xout, dt*xout)
  return dt*xout



def run_simul(prjname, srcdir, gro, top, mdp):
  ''' given a .gro file, prepare system for room temperature '''
  dir, prep, prjroot = mkprjdir(prjname, srcdir) # make a new directory

  # A. prepare necessary files
  os.rename(srcdir+"/"+gro, prep+"/input.gro") # move .gro
  shutil.copy2(srcdir+"/"+top, prep+"/"+top) # copy top
  # prepare mdrun.mdp
  src = mkmdp(srcdir+"/"+mdp, nsteps_warmup)
  open("mdrun.mdp", "w").write(src)
  open("em.mdp", "w").write(file_em_mdp % nsteps_enemin)

  capture = (verbose == 0)

  # B. make a naked structure a.gro
  tpr = glob.glob(srcdir+"/md*.tpr")[0]
  cmd = 'echo "1" | trjconv -s %s -f input.gro -o a.gro' % (tpr)
  runcmd(cmd, system = 1)

  # C. energy minimize the naked
  emgro = "em0.gro"
  # make a naked.top
  arr = open(top, "r").readlines()
  src = ""
  for i in range(len(arr)):
    s = arr[i]
    if i > len(arr)-4 and (s.startswith("SOL")
        or s.startswith("Na") or s.startswith("Cl")):
      continue
    src += s
  open("naked.top", "w").write(src)
  ret = runcmd("grompp -v -f em.mdp -o em0.tpr -c a.gro -p naked.top", capture)
  die_if(ret[0] != 0, "failed to make em topology")
  ret = runcmd("mdrun -v -s em0.tpr -o em0 -c " + emgro + " -x em0 -g em0 -e em0", capture)
  die_if(ret[0] != 0, "failed energy minimization")
  # create a naked tpr
  ret = runcmd("grompp -v -f em.mdp -o naked.tpr -c em0.gro -p naked.top", capture)
  die_if(ret[0] != 0, "failed to create naked.tpr")


  # D. warmup
  ret = runcmd("grompp -v -f mdrun.mdp -o 300.tpr -c input.gro -p topol.top", capture)
  die_if(ret[0] != 0, "failed to make warmup topology")
  ret = runcmd("mdrun -v -s 300.tpr -o 300 -c 300.gro -x 300 -g 300 -e 300", capture)
  die_if(ret[0] != 0, "failed to warmup")

  # E. cleanup
  os.system("rm -f \#* mdout.mdp")

  # present
  os.chdir(dir)
  shutil.copy2("prep/300.gro", "init.gro")
  shutil.copy2("prep/topol.top", "topol.top")
  shutil.copy2("prep/em0.gro", "%s.gro" % prjname)
  shutil.copy2("prep/naked.tpr", "%s.tpr" % prjname)
  shutil.copy2("prep/naked.top", "naked.top")
  src = mkmdp(srcdir+"/"+mdp, 2000000000)
  open("mdrun.mdp", "w").write(src)
  src = mkpbs(prjroot, prjname)
  open("foo.pbs", "w").write(src)



def mkprjdir(prjname, srcdir):
  ''' build a new directory for prjname '''
  s = srcdir
  i = s.find("trj")
  if i < 0:
    print "cannot determine where to put a directory"
    raise Exception
  dir = s[:i] + prjname
  if os.path.exists(dir):
    ret = raw_input("[%s] exists, about to erase it, okay?" % dir)
    ret = ret.strip()
    if len(ret) > 0 and ret in ("n", "N"): exit(1)
    shutil.rmtree(dir)
  root = s[:i-1]
  j = root.rfind("/")
  if j >= 0:
    root = root[j+1:]
  else:
    print "no project name root = %s" % root
    raise Exception

  os.mkdir(dir)
  dirprep = dir + "/prep"
  os.mkdir(dirprep)
  os.chdir(dirprep)
  return dir, dirprep, root



def mkmdp(model, nsteps, temp = "300"):
  src = ""
  temp = str(temp)
  for s in open(model).readlines():
    s = s.strip()
    if s.startswith("ref_t"):
      i = s.find("=")
      tarr = s[i+1:].strip().split()
      s = s[:i+1] + (" "+temp)*len(tarr)
    elif s.startswith("gen_temp"):
      s = "gen_temp = "+temp
    elif s.startswith("nsteps"):
      s = "nsteps = %s" % nsteps
    elif s.startswith("gen_seed"):
      s = "gen_seed = %d" % random.randint(0, 100000000)
    src += s + '\n'
  return src



def runcmd(input, capture = 0, system = 0, verbose = 1):
  ''' run a system command and optionally capture standard output/error
  return a tuple (code, stdout, stderr)
  1. the latter two are None unless `capture' is set, and `system' is 0
  2. to use os.system() instead of subprocess, set `system' to 1
  3. if `verbose' is set, the command is echoed before executed '''
  if capture:
    pipe = subprocess.PIPE
  else:
    pipe = None

  # detect if the input is a string or not
  if type(input) == str:
    cmdstr = input
    cmd = cmdstr.split()
  else:
    cmd = input
    cmdstr = ''.join([s+' '  for s in cmd])

  if verbose >= 1:
    print "CMD:", cmdstr
    if verbose >= 2: raw_input("proceed?")

  if system:
    retcode = os.system(cmdstr)
    oe = ["", ""]  # no stdout or stderr
  else:
    p = subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe = p.communicate()
    retcode=p.returncode
  return (retcode, oe[0], oe[1])



def die_if(cond, message = ""):
  ''' die if `cond' is true '''
  if cond:
    print "fatal error:", message
    raise Exception



def usage():
  ''' print usage and die '''
  print "%s prjname dataid time" % sys.argv[0]
  print '''Options
  -W: time to warmup
  -E: time to energy minimization
  -v: verbose
  '''
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvdW:E:",
        ["help", "verbose=", "warmup=", "enemin="])
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global verbose, nsteps_warmup, nsteps_enemin

  for o, a in opts:
    if o in ("-v", "--verbose",):
      verbose = 1
    elif o in ("-W", "--warmup"):
      nsteps_warmup = int(a)
    elif o in ("-E", "--enemin"):
      nsteps_enemin = int(a)
    elif o in ("-h", "--help"):
      usage()

  if len(args) < 3:
    usage()
  return args[0], args[1], args[2]

file_em_mdp = '''
define              =  -DFLEX_SPC
constraints         =  none
integrator          =  steep
nsteps              =  %s
emtol               =  100
emstep              =  0.01
nstcomm             =  1
ns_type             =  grid
rlist               =  1.2
rcoulomb            =  1.2
rvdw                =  1.0
rvdw_switch         =  0.8
vdwtype             =  shift
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
coulombtype = PME
fourierspacing = 0.12
pme_order = 4
ewald_rtol = 1e-5
'''

def mkpbs(prjroot, prjname):
  return r'''#!/bin/bash
#PBS -N %s
#PBS -l nodes=1:ppn=128,walltime=8:00:00
#PBS -W x=NACCESSPOLICY:SINGLEJOB
#PBS -V
#PBS -j oe

export GMX_MAXCONSTRWARN=-1

username=`whoami`
prjroot=%s
# which program to use
prog=mdrun
# name of the project
prj=mdrun
jid=%s
# the parent directory of gmx and gmxmpi
workhome=$HOME
prjhome=$PROJECTS/jpma/`whoami`
homedir=$prjhome/$prjroot/$jid
# setup the corresponding running directory
#scratch1=$TMPDIR/$username
scratch1=$SHARED_SCRATCH/$username
scratch2=$scratch1/$prjroot
scratch3=$scratch2/$jid
if ! [ -d $scratch1 ] ; then mkdir $scratch1 ; fi
if ! [ -d $scratch2 ] ; then mkdir $scratch2 ; fi
if ! [ -d $scratch3 ] ; then mkdir $scratch3 ; fi
rundir=$scratch3

# clean up the rundir directory
if [ -d $rundir ] ; then rm -f $rundir/* ; fi

# now run simulation
echo "My job ran on: "
cat $PBS_NODEFILE
cd $homedir
$homedir/genrank.sh $PBS_NODEFILE > rankfile
pwd

# construct the backup dir
dataid=1
echo "I now create directory $backdir"
while [ -d data$dataid ] ; do dataid=$(($dataid+1))  ; done
backdir=$homedir/data$dataid
# $dataid is the backup directory
mkdir $backdir

if [ $dataid -gt 1 ] ; then
  echo "a continued run."
  OPTCPT="-cpi state.cpt"
  isnew=0
else
  echo "fresh new run."
  OPTCPT=" "
  # create topology file
  cd $homedir
  $workhome/gmx/src/kernel/grompp -f $prj.mdp -o $prj.tpr -c init.gro -p topol.top
  isnew=1
fi

# back up the initial conditions
if [ $isnew -eq 1 ] ; then
  cp $homedir/$prj.tpr  $backdir/
  mv $homedir/mdout.mdp $backdir/
else
  if [ -f $homedir/state.cpt ] ; then cp $homedir/state.cpt $backdir/state0.cpt ; fi
fi
cp $homedir/$prj.mdp $backdir/

# copy necessary files to the rundir directory
if [ -f $homedir/$prj.tpr ] ; then
  cp $homedir/$prj.tpr  $rundir/
else
  exit 1
fi

if [ $isnew -eq 0 ] ; then # we use checkpoint only if start from an existing trajectory
  if [ -f $homedir/state.cpt ] ; then cp $homedir/state.cpt $rundir/  ; fi
fi

# run the program in the rundir directory
cd $rundir
echo "running $prog ..."
mpirun -np 60 -rf $homedir/rankfile \
  $workhome/gmxmpi/src/kernel/$prog -npme 10 \
  -maxh 7.9 -s $prj $OPTCPT -e $prj -o $prj -x $prj -g $prj -c final.gro &> $prj.out
echo "done"

cd $rundir
cp $rundir/${prj}*.trr   $backdir/
cp $rundir/${prj}*.edr   $backdir/
cp $rundir/${prj}*.log   $backdir/
cp $rundir/${prj}*.xtc   $backdir/
if [ -f final.gro  ] ; then cp $rundir/final.gro  $backdir/  ; fi
if [ -f grompp.out ] ; then cp $rundir/grompp.out $backdir/  ; fi
if [ -f $prj.out   ] ; then cp $rundir/$prj.out   $backdir/  ; fi
if [ -f state.cpt  ] ; then cp $rundir/state.cpt  $backdir/  ; fi

# for rerun we copy the average file and final structure to homedir as well
if [ -f state.cpt  ] ; then cp $rundir/state.cpt  $homedir/  ; fi
''' % (prjname, prjroot, prjname)

if __name__ == "__main__":
  prjname, id, tm = doargs()
  run(prjname, id, float(tm))


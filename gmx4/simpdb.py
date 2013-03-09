#!/usr/bin/env python

''' run simulation from a pdb (amberized) '''
from math import *
from copy import copy
import sys, os, subprocess, getopt, shutil, re, random, glob

d2g = pi/180   # degree to radian

pdbname = None
initgro = "init.gro"
boxsize = 6.5
nsteps_warmup = 1000 # number of steps during 300K warm up simulation
nsteps_runtime = 2000000000 # number of steps in the actual simulation
donaked = 0
verbose = 0

# modify for GROMACS path here
gmxbuild = gmxsrc = kernel = tools = top = mdrun = pdb2gmx = grompp = editconf = genbox = genion = ""

def main():
  doargs()
  globpath() # determine path

  # write a few files
  open("em.mdp", "w").write(file_em_mdp)

  boxgro, charge = mkgro(pdbname, boxsize)  # 1. make a PDB file, call pdb2gmx
  open("300.mdp", "w").write(mkmdp(nsteps_warmup, charge))
  if donaked:
    enemin0(boxgro)
  solgro = addwater(boxgro, charge) # 2. add water
  emgro = enemin(solgro)  # 3. energy minimization
  roomgro = warmup(emgro) # 4. warm up to room temperature

  # copy files from prep. to project dir.
  shutil.copy2("300.gro",  initgro)
  open("mdrun.mdp", "w").write(mkmdp(nsteps_runtime, charge))
  #open("foo.pbs", "w").write(mkpbs(molname))

def mkgro(pdbname, boxsize):
  ''' make an initial pdb structure of a single chain
  then call pdb2gmx so a topology of a single chain is made '''

  groname = "a.gro"
  ret, out, err = runcmd(pdb2gmx + " -f " + pdbname + " -o " + groname + " -ff amber03star -ignh",
      capture = 1)

  # search for charge
  chargekey = "Total charge "
  for s in err.splitlines():
    if s.startswith(chargekey):
      s = s.strip()[len(chargekey):-1]
      q = int(round(float(s)))
      break
  else:
    print "cannot determine charge"
    raise Exception
  print "converted to GROMACS format %s, charge is %d" % (groname, q)

  # change spc to ffamber_tip3p
  top = "topol.top"
  s = open(top).read()
  s = re.sub('"spc.itp"', '"ffamber_tip3p.itp"', s)
  open(top, "w").write(s)

  # fit to a box
  boxgro = "box.gro"
  ret = runcmd('echo 0 | ' +
      editconf + " -princ -f " + groname + " -o " + boxgro + " -box %s " % boxsize + " -bt dodec",
      system = 1)
  return boxgro, q

def addwater(boxgro, charge, top = "topol.top"):
  ''' add water and ion'''
  solgro = "solve.gro"
  capture = (verbose == 0)
  ret = shrun(genbox + " -cp " + boxgro + " -cs ffamber_tip3p.gro -o " + solgro + " -p " + top, capture)

  if charge:
    # handle ions
    iongro = "ion.gro"
    if charge > 0:
      opt = "-nname Cl -nn %d" % charge
      iontop = ["...",]
      iontop = ["[ moleculetype ]",
                "Cl %d" % charge,
                "",
                "[ atoms ]",
                "1       amber99_30      1      Cl-     Cl        1     -1    35.45",
                ""]
    else:
      opt = "-pname Na -np %d" % (-charge)
      iontop = ["[ moleculetype ]",
                "Na %d" % (-charge),
                "",
                "[ atoms ]",
                "1       amber99_31      1      Na+      Na        1     1    22.99",
                ""]
    ret = shrun(grompp + " -f em.mdp -o ion.tpr -c " + solgro, capture)
    ret = runcmd("echo 12 | " +
        genion + " -s ion.tpr -o " + iongro + " -g ion -p topol.top " + opt,
        system = 1)

    # modify the topology file
    s = open(top).read()
    s = re.sub("Cl\s*0\n", "", s)
    s = re.sub("Na\s*0\n", "", s)
    arr = s.splitlines()
    for i in range(len(arr)):
      if arr[i].startswith("; Include generic topology"):
        arr = arr[:i] + iontop + arr[i:]
        break
    else:
      print "cannot find where to add ion topology"
      raise Exception
    open(top, "w").write('\n'.join(arr) + '\n')
  else:
    iongro = solgro

  return iongro

def enemin(solgro):
  ''' energy minimization '''
  emgro = "em.gro"
  capture = (verbose == 0)
  ret = shrun(grompp + " -v -f em.mdp -o em.tpr -c " + solgro, capture)
  ret = runcmd(mdrun + " -v -s em.tpr -o em -c " + emgro + " -x em -g em -e em", capture)
  die_if(ret[0] != 0, "failed energy minimization")
  return emgro

def enemin0(boxgro):
  ''' energy minimization '''
  em0gro = "em0.gro"
  capture = (verbose == 0)
  ret = shrun(grompp + " -v -f em.mdp -o em0.tpr -c " + boxgro, capture)
  ret = runcmd(mdrun + " -v -s em0.tpr -o em0 -c " + em0gro + " -x em0 -g em0 -e em0", capture)
  die_if(ret[0] != 0, "failed energy minimization")
  return em0gro

def warmup(emgro):
  ''' warm-up simulation '''
  roomgro = "300.gro"
  capture = (verbose == 0)
  ret = shrun(grompp + " -v -f 300.mdp -o 300.tpr -c " + emgro, capture)
  ret = runcmd(mdrun + " -v -s 300.tpr -o 300 -c " + roomgro + " -x 300 -g 300 -e 300", capture)
  die_if(ret[0] != 0, "failed in room-temperature warmup")
  return roomgro

def usage():
  ''' print usage and die '''
  print "Usage:\n  %s [OPTIONS] your.pdb\n" % (sys.argv[0])
  print '\n'.join(["OPTIONS:",
        " -b: --box=:     followed by the box size",
        " -o: output file",
        " -W: --warmup=:  followed by the number of steps in the warm-up simulation",
        " -R: --runtime=: followed by the actual simulation runtime",
        " --donaked: energy minimize the bare structure",
        " -h: help"])
  exit(1)

def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvW:d:b:o:",
        ["help", "verbose=", "box=", "output=", "warmup=", "donaked"])
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global verbose, donaked
  global pdbname, initgro
  global boxsize, nsteps_warmup, nsteps_runtime

  pdbname = "out.pdb"
  for o, a in opts:
    if o in ("-v", "--verbose",):
      verbose = 1
    elif o in ("-b", "--box",):
      boxsize = float(a)
    elif o in ("-o", "--output="):
      initgro = a
    elif o in ("-W", "--warmup",):
      nsteps_warmup = int(a)
    elif o in ("-R", "--runtime",):
      nsteps_runtiem = int(a)
    elif o in ("--donaked",):
      donaked = 1
    elif o in ("-h", "--help"):
      usage()

  if len(args) > 0:
    pdbname = args[0]
  else:
    arr = glob.glob('*.pdb')
    if len(arr) > 0:
      pdbname = arr[0]
      raw_input("I assume we use %s?" % pdbname)
    else:
      print "cannot find any pdb under the current directory"
      usage()


def globpath():
  ''' set up global paths '''
  build = "$HOME/work/gromacs/build"
  src = "$HOME/work/gromacs"

  global gmxbuild, gmxsrc, kernel, tools, top, mdrun, pdb2gmx, grompp, editconf, genbox, genion
  gmxbuild  = mkp(build)
  gmxsrc    = mkp(src)
  kernel    = mkp(gmxbuild, "src/kernel")
  tools     = mkp(gmxbuild, "src/tools")
  top       = mkp(gmxsrc, "share/top")
  mdrun     = mkp(kernel, "mdrun")
  pdb2gmx   = mkp(kernel, "pdb2gmx")
  grompp    = mkp(kernel, "grompp")
  editconf  = mkp(tools, "editconf")
  genbox    = mkp(tools, "genbox")
  genion    = mkp(tools, "genion")

def mkp(prefix, path = None, make = 0):
  ''' return a new path as prefix + path, expand $HOME
  create the path as a directory, if `make' is true '''
  if path:
    newpath = prefix
  else:  # prefix is treated as path if path is missing
    if prefix[0] == '/':
      newpath = '/'
      path = prefix[1:]
    else:
      newpath = ""
      path = prefix
  for s in path.split('/'):
    if s and s[0] == "$": # handle $HOME
      s = os.getenv(s[1:])
    newpath = os.path.join(newpath, s)
  if make == 1: # create the path
    if os.path.exists(newpath): # remove existing
      shutil.rmtree(newpath)
    os.mkdir(newpath)
    print "made directory: ", newpath
  return newpath

def die_if(cond, message = ""):
  ''' die if `cond' is true '''
  if cond:
    print "fatal error:", message
    raise Exception

def shrun(input, capture = 0, verbose = 1):
  ''' write a shell script to run the program,
  so environment variables can be defined '''
  script = "tmp.sh"
  s = "#!/bin/bash\nexport GMXLIB=%s\n%s\n" % (top, input)
  open(script, "w").write(s)
  os.system("chmod 755 %s" % script) # make the script runnable
  if verbose:
    print "CMD:", input
  ret = runcmd("./" + script, capture, verbose = 0)
  die_if(ret[0] != 0, "error occurred when running %s" % input)
  return ret

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

def getcmdout(cmd):
  ''' capture output of a command, like whoami '''
  tmpnm = "TMP"+str(random.randint(1, 99999))
  ret = runcmd(cmd+" > "+tmpnm, system = 1, verbose = 0)[0]
  die_if(ret != 0, "command [%s] failed, ret = %d" % (cmd, ret))
  out = open(tmpnm).read().strip()
  os.remove(tmpnm)
  if out: return out.split('\n')[0].strip()
  else: return ""

def mkmdp(nsteps, charge):
  ''' make a .mdp file '''
  if charge == 0:
    iongrp = taut = reft = ""
  else:
    if charge > 0:
      iongrp = "Cl-"
    else:
      iongrp = "Na+"
    taut = "0.1"
    reft = "300"
  return '''constraints         =  hbonds
integrator          =  md
dt                  =  1e-3
nsteps              =  %d
nstcomm             =  1
nstxtcout           =  1000
xtc_grps            =  System
nstxout             =  0
nstvout             =  0
nstfout             =  0
nstlog              =  1000
nstenergy           =  1000
nstlist             =  10
ns_type             =  grid
rlist               =  1.2
rcoulomb            =  1.2
vdwtype             =  shift
rvdw                =  1.0
rvdw_switch         =  0.8
Tcoupl              =  v-rescale
tc-grps	            =  Protein SOL %s
tau_t               =  0.1   0.1 %s
ref_t               =  300   300 %s
energygrps          =  Protein SOL %s
Pcoupl              =  no
gen_vel             =  yes
gen_temp            =  300.0
gen_seed            =  %d
coulombtype = PME
fourierspacing = 0.144
pme_order = 4
ewald_rtol = 1e-5
''' % (nsteps,
    iongrp, taut, reft, iongrp,
    random.randint(0, 1000000000))

file_em_mdp = '''
define              =  -DFLEX_SPC
constraints         =  none
integrator          =  steep
nsteps              =  100
emtol               =  2000
emstep              =  0.01
nstcomm             =  1
ns_type             =  grid
rlist               =  1
rcoulomb            =  1.0
rvdw                =  1.0
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
'''

def mkpbs(molname):
  ''' make foo.pbs for room temperature simulation '''
  advisor = "jpma"
  exeker = "$HOME/gmxmpi/src/kernel"
  gmxtop = "$HOME/gromacs/share/top"
  #exeker = "$HOME/gmx45/build/src/kernel"
  #gmxtop = "$HOME/gmx45/gromacs/share/top"
  return '''#!/bin/bash
#PBS -N %s
#PBS -l nodes=1:ppn=8,walltime=24:00:00
#PBS -V
#PBS -j oe

nthreads=8
# additional PBS   -W x=NACCESSPOLICY:SINGLEJOB

prjroot=pegda
# which program to use
prog=mdrun
# name of the project
prj=mdrun
jid=%s
# the parent directory of gmx and gmxmpi
workhome=$HOME
# where to find executables
exeker=%s
# make sure to find the correct gromacs
export GMXLIB=%s

prjhome=$PROJECTS/%s/`whoami`
homedir=$prjhome/$prjroot/$jid
username=`whoami`
# setup the corresponding running directory
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
pwd

# construct the backup dir
cd $homedir
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
  $exeker/grompp -f $prj.mdp -o $prj.tpr -c init.gro -p topol.top
  isnew=1
fi

# back up the initial conditions
if [ $isnew -eq 1 ] ; then
  cp $homedir/$prj.tpr  $backdir/
  mv $homedir/mdout.mdp $backdir/
else
  if [ -f $homedir/state.cpt ] ; then cp $homedir/state.cpt $backdir/state0.cpt ; fi
fi
if [ -f $homedir/MTSEED    ] ; then cp $homedir/MTSEED    $backdir/MTSEED0    ; fi
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
mpiexec -np $nthreads $exeker/$prog -maxh 24.0 -s $prj $OPTCPT -e $prj -o $prj -x $prj -g $prj -c final.gro &> $prj.out
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
''' % (molname, molname, exeker, gmxtop, advisor)

if __name__ == "__main__":
  main()


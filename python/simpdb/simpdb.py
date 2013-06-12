#!/usr/bin/env python

''' prepare a system for MD simulation from a pdb for GROMACS 4.0/4.5
    Copyright (c) 2010-2013, Cheng Zhang '''

import sys, os, subprocess, getopt, shutil, re, random, glob, math
import zcom
import gmxcom


d2g = math.pi/180   # degree to radian
fntop = "topol.top" # default topology file

boxsize = 6.5
nsteps_warmup = 100 # number of steps during 300K warm up simulation
nsteps_runtime = 2000000000 # number of steps in the actual simulation
emnaked = False
dopbs = False
verbose = 0

# GROMACS paths
gmxexe = gmxsrc = gmxtopdir = gmxver = None
mdrun = "mdrun"
pdb2gmx = "pdb2gmx"
grompp = "grompp"
editconf = "editconf"
genbox = "genbox"
genion = "genion"



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input.pdb""" % sys.argv[0]

  print """
  prepare a system for GROMACS MD simulation from a .pdb file
  Copyright (c) 2010-2013, Cheng Zhang

  OPTIONS:
   -o               followed by output .gro file
   -b, --box=       followed by the box size
   -W, --warmup=    followed by the number of steps in the 300K warm-up simulation
   -R, --runtime=   followed by the actual simulation runtime
   -B, --gmxexe=    root directory for GROMACS building tree
   -S, --gmxsrc=    root directory for GROMACS source code tree
   --ver=           GROMACS version, if --gmxsrc is not given
   --ff=            force field
   --water=         water model
   -i, --implicit   implicit solvent
   --mdp=           string of .mdp options separated by commas, like "a=1, b=2"
   --emnaked        energy minimize the naked structure
   --dopbs          produces a .pbs script
   -v               be verbose
   --verbose=       set verbocity
   -h, --help       help

  Note, for GROMACS 4.0.x, the input .pdb must go through `amberize' first
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvW:d:b:o:B:S:i",
        [ "help", "verbose=",
          "gmxexe=", "gmxsrc=", "ver=", "version=",
          "ff=", "water=", "mdp=", "implicit=",
          "box=", "output=", "warmup=", "emnaked"])
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global gmxexe, gmxsrc, gmxver, verbose, emnaked, dopbs
  global boxsize, nsteps_warmup, nsteps_runtime

  ff = "amber99sb-ildn"
  water = "tip3p"
  fngro = "init.gro"
  sver = None
  solvent = "explicit"
  mdparams = {}

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-b", "--box",):
      boxsize = float(a)
    elif o in ("-o", "--output",):
      fngro = a
    elif o in ("-W", "--warmup",):
      nsteps_warmup = int(a)
    elif o in ("-R", "--runtime",):
      nsteps_runtime = int(a)
    elif o in ("--ff",):
      ff = a
    elif o in ("--water",):
      water = a
    elif o in ("-i", "--implicit",):
      solvent = "implicit"
      emnaked = True
    elif o in ("--mdparams",):
      # parse the .mdp string
      for s in a.strip(", ").split():
        k, v = s.split("=")
        mdparams[k.strip()] = v.strip()
    elif o in ("--emnaked",):
      emnaked = True
    elif o in ("--dopbs",):
      dopbs = True
    elif o in ("-B", "--gmxexe",):
      gmxexe = a
    elif o in ("-S", "--gmxsrc",):
      gmxsrc = a
    elif o in ("--ver", "--version",):
      sver = a
    elif o in ("-h", "--help"):
      usage()

  if sver:  # GROMACS version
    if "." in sver:
      sver = sver.replace(".", "0")
    if len(sver) < 5: sver += "0" * (5 - len(sver))
    gmxver = int(sver)

  if len(args) > 0:
    fnpdb = args[0]
  else:
    print "need a pdb file"
    usage()

  return fnpdb, fngro, ff, water, solvent, mdparams



def main():
  global gmxexe, gmxsrc, dopbs

  fnpdb, fngro, ff, water, solvent, mdp = doargs()
  # use particle decomposition for implicit solvent simulation
  pd = (solvent == "implicit")
  setuppaths(gmxexe, gmxsrc) # set up paths
  if gmxver < 40500 and solvent == "implicit":
    print "implicit solvent are not supported for GROMACS 4.0 or lower (%s)" % gmxver
    raise Exception
  capture = (verbose == 0) # if not verbose, we suppress the output

  # convert `fnpdb' to the absolute path, if it is a path
  if os.path.exists(fnpdb):
    fnpdb = os.path.realpath(fnpdb)

  # create a working directory and switch into it
  prepdir = "prep"
  zcom.mkpath(prepdir, force = True)
  os.chdir(prepdir)

  # 1. convert .pdb to .gro file
  wfnpdb = mkpdb(fnpdb)
  boxgro, charge = mkgro(wfnpdb, boxsize, ff, water, fntop)

  # write a few .mdp files
  open("em.mdp", "w").write(mkemmdp(pd))
  open("300.mdp", "w").write(
      mkmdmdp(nsteps_warmup, charge, solvent, params = mdp))
  if emnaked:
    boxgro = runmd(boxgro, "em0", "em.mdp", fntop, pd, capture)

  if solvent == "explicit":
    # 2. add water
    solgro = addwater(boxgro, charge, fntop, capture)

    # 3. energy minimization to remove the friction between protein and water
    emgro = runmd(solgro, "em", "em.mdp", fntop, pd, capture)
  else:
    emgro = boxgro

  # 4. warm up to room temperature
  roomgro = runmd(emgro, "300", "300.mdp", fntop, pd, capture)

  # 5. switch out of the working directory
  os.chdir(os.pardir)
  # copy the initial coordinates and topology files
  shutil.copy2(os.path.join(prepdir, roomgro), fngro)
  shutil.copy2(os.path.join(prepdir, fntop), fntop)
  open("mdrun.mdp", "w").write(
      mkmdmdp(nsteps_runtime, charge, solvent, params = mdp))
  if dopbs:
    open("foo.pbs", "w").write(mkpbs("foo", molname, fngro, fntop))



def mkpdb(fnpdb):
  ''' copy PDB to the running directory '''

  if os.path.exists(fnpdb):
    # copy `fnpdb' to the current directory
    fnpdb0 = os.path.basename(fnpdb)
    shutil.copy2(fnpdb, fnpdb0)
    print "copying %s to %s" % (fnpdb, fnpdb0)
    return fnpdb0
  else:
    # `fnpdb' is the standard 4-letter PDB code
    if not fnpdb.endswith(".pdb") and len(fnpdb) == 4:
      fnpdb = fnpdb.upper() + ".pdb"

    # try to download from RCSB
    if zcom.runcmd("wget http://www.rcsb.org/pdb/files/" + fnpdb,
        system = True)[0] != 0:
      print "cannot download", fnpdb
      raise Exception
    return fnpdb



def mkgro(fnpdb, boxsize, forcefield, sol, fntop):
  ''' convert the initial .pdb to a .gro file, with an initial topology file '''

  global gmxver

  fngro = "a.gro"
  cmd = pdb2gmx + " -f " + fnpdb + " -o " + fngro + " -ff " + forcefield + " -ignh -p " + fntop
  if gmxver >= 40500:
    cmd +=  " -water " + sol
  ret, out, err = zcom.runcmd(cmd, capture = 1)

  # search for charge
  chargekey = "Total charge "
  for s in err.splitlines():
    if s.startswith(chargekey):
      s = s.strip()[len(chargekey):-1]
      q = int(round(float(s)))
      break
  else:
    print "cannot determine the charge, error report:\n%s" % err
    raise Exception
  print "converted to GROMACS format %s, charge is %d" % (fngro, q)

  if gmxver < 40500:
    # for GROMACS 4.0, write .itp file for tip3p
    s = open(fntop).read()
    s = re.sub('"spc.itp"', '"ffamber_tip3p.itp"', s)
    open(fntop, "w").write(s)

  # fit to a box
  boxgro = "box.gro"
  ret = zcom.runcmd('echo 0 | ' +
      editconf + " -princ -f " + fngro + " -o " + boxgro
               + " -box " + str(boxsize) + " -bt dodec",
      system = 1)
  return boxgro, q



def addwater(boxgro, charge, fntop, capture):
  ''' add water molecules and ions'''

  global gmxtopdir, gmxtop

  solgro = "solve.gro"
  ret = gmxshrun(genbox + " -cp " + boxgro +
      " -cs spc216.gro -o " + solgro + " -p " + fntop,
      capture)

  if charge:
    # handle ions
    iongro = "ion.gro"
    if charge > 0:
      opt = " -nn %d " % charge
      if gmxver < 40500: opt += " -nname Cl "
    else:
      opt = " -np %d " % (-charge)
      if gmxver < 40500: opt += " -pname Na "

    ret = gmxshrun(grompp + " -f em.mdp -o ion.tpr -c " + solgro,
        capture)

    SOLgroup = 13
    if gmxver < 40500: SOLgroup = 12

    ret = zcom.runcmd("echo %s | " % SOLgroup +
        genion + " -s ion.tpr -o " + iongro + " -p " + fntop + opt,
        system = 1)

    if gmxver < 40500:
      # for GROMACS 4.0
      iontop = '''
[ moleculetype ]
Cl 1

[ atoms ]
1       amber99_30      1      Cl-     Cl        1     -1    35.45

[ moleculetype ]
Na 1

[ atoms ]
1       amber99_31      1      Na+     Na        1      1    22.99
'''
      s = open(fntop).read()
      s = re.sub("Cl\s*0\n", "", s)
      s = re.sub("Na\s*0\n", "", s)
      arr = s.splitlines(True)
      for i in range(len(arr)):
        if arr[i].startswith("; Include generic topology"):
          arr = arr[:i] + iontop.splitlines(True) + arr[i:]
          break
      else:
        print "cannot find where to add ion topology %s" % fntop
        raise Exception
      open(fntop, "w").write(''.join(arr))
  else:
    iongro = solgro
  return iongro



def runmd(fningro, fnout, fnmdp, fntop, pd, capture):
  ''' run mdrun with the structure in `fningro' '''

  if fnout.endswith(".gro"):
    fnout = os.path.splitext(fnout)[0]

  ret = gmxshrun(grompp + " -v -f " + fnmdp + " -o " + fnout + ".tpr "
                 + "-c " + fningro + " -p " + fntop,
                 capture)

  print "\n\nRunning mdrun with %s, may take some time..." % fnout
  cmd = mdrun + " -v -deffnm " + fnout
  # if `ns_type' is `simple' instead of `grid'
  # we have to use particle decomposition
  if pd: cmd += " -pd"
  ret = zcom.runcmd(cmd, capture)
  zcom.die_if(ret[0] != 0, "mdrun failed with %s" % fnout)
  return fnout + ".gro"



def myfind(fn, root):
  ''' find `fn' under directory `root' '''
  ls = zcom.pathglob([fn], root, True, links = True)
  if len(ls): return ls[0]
  else: return fn



def setuppaths(buildroot, srcroot):
  ''' set up global paths
      if `buildroot' is empty, use programs under the default path
      if `srcroot' is empty, deduce from the current path '''

  global gmxexe, gmxsrc, gmxtopdir, gmxver
  global mdrun, pdb2gmx, grompp, editconf, genbox, genion

  if buildroot: buildroot = os.path.expanduser(buildroot)

  # we prefer to use tools under a development tree
  # for they are more up-to-date, but if it is impossible
  # use the default programs
  if buildroot and os.path.exists(buildroot):
    gmxexe  = os.path.realpath(buildroot)
    print "setting the executable directory as", gmxexe
    mdrun     = myfind("mdrun",     gmxexe)
    pdb2gmx   = myfind("pdb2gmx",   gmxexe)
    grompp    = myfind("grompp",    gmxexe)
    editconf  = myfind("editconf",  gmxexe)
    genbox    = myfind("genbox",    gmxexe)
    genion    = myfind("genion",    gmxexe)

  # gmxsrc and gmxtopdir are not necessary
  gmxsrc = gmxcom.setsrcroot(srcroot)
  if gmxsrc:
    # use the topology directory under the source tree for it is more up-to-date
    gmxtopdir = os.path.join(gmxsrc, "share", "top")
    # if gmxsrc is obtained successful, we can determine the version
    if not gmxver: gmxver = gmxcom.version()

  if not gmxver:
    gmxver = 40507
    print "don't GROMACS version, assuming %s" % gmxver



def gmxshrun(cmd, capture = 0, verbose = 1):
  ''' write a shell script to run the program, s.t. GMXLIB is defined '''
  global gmxtopdir
  envars = {}
  if gmxtopdir and os.path.exists: envars["GMXLIB"] = gmxtopdir
  return zcom.shrun(cmd, capture, verbose, envars = envars)



def mkmdmdp(nsteps, charge, solvent = None,
    Tref = 300, taut = 0.1, params = {}):
  ''' make a .mdp file suitable for molecular dynamics simulation
      if `charge' is not an integer, we use both ions
      `solvent' can be "implicit" or "explicit" '''

  d = {}
  for k in params: d[k] = params[k]

  if solvent == "implicit":
    d["implicit_solvent"] = "GBSA"
    d["gb_algorithm"] = "OBC"

    d["comm-mode"] = "Angular"
    d["ns_type"] = "Simple"
    d["nstlist"] = 0
    d["pbc"] = "no"

    # 0 cut-off mean no cut-off
    d["rlist"] = 0.0
    d["rcoulomb"] = 0.0
    d["rgbradii"] = 0.0
    d["nstgbradii"] = 1

    d["rvdw_switch"] = None
    d["vdwtype"] = "Cut-off"
    d["coulombtype"] = "Cut-off"
    d["optimize_fft"] = "yes"
    # 0.0054 kcal/A^2/mol = 2.259 kJ/nm^2/mol
    d["sa_surface_tension"] = 2.25936 # may be -1
  elif solvent == "explicit":
    # handle charges
    Na, Cl = "NA", "CL"
    global gmxver
    if gmxver < 40500:
      Na, Cl = "Na+", "Cl-"
    if charge > 0:
      iongrp = Cl
    elif charge < 0:
      iongrp = Na
    else:
      iongrp = ""

    # create seperate groups
    for grp in ("xtc_grps", "tc-grps", "energygrps"):
      d[grp] = "Protein SOL " + iongrp
      ngrp = len( d[grp].strip().split() )

    d["ref_t"] = ' '.join( [str(Tref), ] * ngrp )
    d["taut"]  = ' '.join( [str(taut), ] * ngrp )

  return gmxcom.mkmdp(d)



def mkemmdp(pd):
  d = {
    "nsteps" : 1000,
    "define" : "-DFLEX_SPC",
    "constraints" : "none",
    "integrator" : "steep",
    "emstep" : 0.01,
    "emtol" : 2000,
    "gen_vel" : "no",
    "coulombtype" : "Cut-off",
  }
  if pd: d["ns_type"] = "simple"
  return gmxcom.mkmdp(d)



def mkpbs(prjroot, molname, fngro, fntop):
  ''' make foo.pbs for room temperature simulation
      this script is provided for convenience '''

  global gmxver

  if gmxver >= 40600:
    srcroot   = "$HOME/gmx/gromacs"
    buildroot = "$HOME/gmx/gromacs/build"
    mdrundir  = buildroot + "/src/programs/mdrun"
    gromppdir = buildroot + "/src/programs/grompp"
  else: # version 4.5 or eariler
    srcroot   = "$HOME/gmx/gromacs45"
    buildroot = "$HOME/gmx/gromacs45/build"
    mdrundir  = buildroot + "src/kernel"
    gromppdir = mdrundir
  gmxtopdir = srcroot + "/share/top"

  d = {
      "prjroot" : prjroot,
      "molname" : molname,
      "mdrundir" : mdrundir,
      "gromppdir" : gromppdir,
      "gmxtopdir" : gmxtopdir,
      "fngro" : fngro,
      "fntop" : fntop,
      }
  return zcom.templrepl('''#!/bin/bash
#PBS -N {{molname}}
#PBS -l nodes=1:ppn=8,walltime=24:00:00
#PBS -V
#PBS -j oe

nthreads=8
# additional PBS -W x=NACCESSPOLICY:SINGLEJOB

prjroot={{prjroot}}
# which program to use
prog=mdrun
# name of the project
prj=mdrun
jid={{molname}}
# the parent directory of gmx and gmxmpi
workhome=$HOME
# where to find executables
mdrun={{mdrundir}}
gromppdir={{gromppdir}}
# make sure to find the correct gromacs
export GMXLIB={{gmxtopdir}}

prjhome=$HOME/work
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
  $gromppdir/grompp -f $prj.mdp -o $prj.tpr -c {{fngro}} -p {{fntop}}
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
''', d)

if __name__ == "__main__":
  main()


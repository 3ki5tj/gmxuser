#!/bin/usr/env python

''' common routines for GROMACS

    export function list:
    * srcroot()           get GROMACS root directory
    * version()           get GROMACS version
    * isver4()            test if it is version 4.0
    * sopenmm()           get OPENMM signature
    * getf(fn)            get path of `fn' from the source tree
    * mkmdp()             make a standard .mdp file
'''


import os, sys, re, random

# global variables
gmxsrcroot = None
gmxver = None
gmxisver4 = None
gmxsopenmm = None



def setsrcroot(root, child = None):
  ''' set up or determine the GROMACS source root directory
      if `root' is not given, it is deduced from the subdirectory `child',
      which is presumably the current directory '''
  global gmxsrcroot, gmxver, gmxisver4, gmxsopenmm

  if root:
    root = os.path.expanduser(root)

  if not root or not os.path.exists(root):
    # guess the GROMACS root from the current directory
    if child == None: child = os.getcwd()
    root = child
    # keep searching the parent directory
    while not (root == "/" or root.endswith(":\\")):
      root = os.path.abspath( os.path.join(root, os.pardir) )
      # there is a file "AUTHORS" under the root
      if os.path.exists( os.path.join( root, "AUTHORS" ) ):
        readme = os.path.join(root, "README")
        # there is a file "README" with "GROMACS" in it
        if ( os.path.exists(readme)
            and "GROMACS" in open(readme).read() ):
          break
    else:
      print "cannot determine the root of GROMACS source, current '%s'" % child
      return None
  gmxsrcroot = os.path.realpath( os.path.expanduser(root) )
  # set other variables to empty to recompute them
  gmxver = gmxisver = gmxsopenmm = None
  return gmxsrcroot



def srcroot():
  ''' get the GROMACS source root directory '''
  global gmxsrcroot
  # try to guess the root
  if not gmxsrcroot and not setsrcroot(None):
    raise Exception
  return gmxsrcroot



def version():
  ''' get GROMACS version from CMakeLists.txt or configure.ac
      should work for GROMACS v4.0+ '''
  global gmxver
  if gmxver != None: return gmxver

  root = srcroot()
  cfgac = os.path.join(root, "configure.ac")
  cmake = os.path.join(root, "CMakeLists.txt")
  vinfo = os.path.join(root, "cmake", "gmxVersionInfo.cmake")

  if os.path.exists(vinfo): # version 5.1
    # we look for GMX_VERSION_MAJOR, GMX_VERSION_MINOR, and GMX_VERSION_PATCH
    myver = 0
    for s in open(vinfo):
      m = re.search(r'\(GMX_VERSION_MAJOR\s+(.*?)\)', s.strip())
      if m:
        myver += int( m.group(1) ) * 10000
        continue

      m = re.search(r'\(GMX_VERSION_MINOR\s+(.*?)\)', s.strip())
      if m:
        myver += int( m.group(1) ) * 100
        continue

      m = re.search(r'\(GMX_VERSION_PATCH\s+(.*?)\)', s.strip())
      if m:
        myver += int( m.group(1) )
        continue
    gmxver = myver

  elif os.path.exists(cmake): # version 4.5 to 5.0
    for s in open(cmake):
      # we look for a line that looks like
      # set(PROJECT_VERSION "4.5.6-dev")
      # here, `*?' in `(.*?)' means a non-greedy search
      #       `(-dev)?' means `-dev' may or may not present
      m = re.search(r'PROJECT_VERSION\s+"(.*?)(-dev)?"', s.strip())
      if m:
        # convert 4.5.6 to an integer 40506
        sver = m.group(1).replace('.', '0')
        if len(sver) < 5: sver += '0' * (5 - len(sver))
        gmxver = int(sver)
        break
    else:
      print "cannot determine version from %s" % cmake
      raise Exception

  elif os.path.exists(cfgac): # v4.0 or earlier
    for s in open(cfgac):
      if s.startswith("AC_INIT"):
        # a typical string is like
        #   AC_INIT(gromacs, 4.5.6-dev, [gmx-users@gromacs.org])
        # so we take the second number in the parentheses ()
        sver = s.split(",")[1].strip()
        if not sver.startswith("4.0"):
          print "bad gromacs version %s" % sver
          raise Exception
        gmxver = int( sver[:5].replace(".", "0") )
        break
    else:
      print "cannot determine version number from %s" % cfgac
      raise Exception

  else:  # there is neither CMakeLists.txt nor configure.ac
    raise Exception

  print gmxver
  return gmxver



def isver4():
  global gmxisver4
  if gmxisver4 != None: return gmxisver4

  gmxisver4 = (version() < 40100)
  return gmxisver4



def sopenmm():
  global gmxsopenmm
  if gmxsopenmm != None: return gmxsopenmm

  gmxsopenmm = "GMX_OPENMM"
  if isver4(): gmxsopenmm = "USE_OPENMM"
  return gmxsopenmm



def getf(fn, paths = [], root = None, panic = True):
  ''' return the full path of file `fn'
      some suggestion are collected in `paths', in which
      the path separator is always `/' no matter the operating system '''

  # default set-up as the source
  if root == None: root = srcroot()
  if type(paths) == str: paths = [ paths ]
  # add common places to paths
  paths = paths + ["src/kernel", "src/mdlib", "src/gmxlib",
      "src/gromacs/gmxlib", "src/gromacs/mdlib",
      "programs/mdrun", "program/pdb2gmx", "program/grompp",
      "include"]
  # find 'fn' in the suggested and common places
  for p in paths:
    path = os.path.join(root, p.replace("/", os.sep), fn)
    if os.path.exists(path): return path

  # glob the GROMACS tree to find the file
  try:
    import zcom
    ls = [a for a in zcom.pathglob([fn], root, recur = True) if os.path.exists(a)]
    if len(ls) > 0: fn = ls[0]
  except ImportError: pass

  if not os.path.exists(fn) and panic:
    print "cannot find path for", fn
    raise Exception
  return fn



def mkmdp(params = {}):
  ''' return a string of a standard .mdp file
     ``params' is a dictionary that adds/overrides parameters
     if the value of a key is set to None, it is removed '''

  # default parameters
  d = {
    "nsteps" : 1000,
    "dt" : 0.002,
    "integrator" : "md",
    "constraints" : "hbonds",
    "nstxtcout" : 1000,
    "nstxout" : 0,
    "nstvout" : 0,
    "nstfout" : 0,

    "nstlist" : 10,
    "nstcalcenergy" : 10,
    "nstcomm" : 10,

    "nstlog"  : 1000,
    "nstenergy" : 1000,
    "ns_type" : "grid",

    "xtc_grps" : "System",
    "tc-grps" : "System",
    "energygrps" : "System",

    "rlist" : 1.2,
    "rcoulomb" :  1.2,
    "vdwtype" : "shift",
    "rvdw" : 1.0,
    "rvdw_switch" : 0.8,

    "taut" : "0.1",
    "ref_t" : "300",
    "gen_seed" : random.randint(0, 1000000000),

    "Tcoupl" : "v-rescale",
    "Pcoupl"  : "no",
    "gen_vel" : "yes",
    "gen_temp" : "300",

    "coulombtype" : "PME",
    "fourierspacing" : 0.144,
    "pme_order" : 4,
    "ewald_rtol" : 1e-5,
  }

  # override standard parameters,
  # add additional parameters are added here
  for k in params:
    if not params[k]: # remove the key
      d.pop(k, None)
    else: # add the key
      d[k] = params[k]

  s = ["; machine-generated GROMACS parameter file\n", ]
  for k in sorted(d):
    s += [ k + " = " + str(d[k]) + "\n" ]
  return ''.join(s)





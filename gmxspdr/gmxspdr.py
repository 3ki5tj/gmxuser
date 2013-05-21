#!/usr/bin/env python

'''
build a self-contained C file as a GROMACS engine
basically it grabs several key files from the GROMACS source tree and
writes an output `mdfoo.c' based on an optional input template `foo.h'
to see more, run the script with `-h', the function usage()
'''

import re, sys, os, glob, getopt
from cc import CC
from ccgmx import CCGMX
from ccaux import CCAUX
import ccmdxv40
import ccmdxv45

# default project name
prjname = "foo"
version = 1

# default input template file
finp = '''
/* $OBJ_DECL */
/* $OBJ_FUNCS */
/* $BONDFREE_C */
/* $FORCE_C */"
/* $SIMUTIL_C */
/* $MD_C */"
/* $RUNNER_C */
/* $MDRUN_C */
'''.splitlines(True)


def usage():
  ''' print usage and die '''
  print "%s [OPTIONS] project" % sys.argv[0]

  print '''
  The program grabs several core GROMACS files, and writes
  a self-contained C file, which can be compiled to an independent engine
  Optionally, it reads an input template file `foo.h',
  and then outputs `mdfoo.c', in which the functional tags, such as
  $MDRUN_C are replaced by the relevant functions in the corresponding
  GROMACS source code (`mdrun.c').
  '''

  print '''
  OPTIONS:
  -v: version number (larger means a more complex template)
  '''
  exit(1)


def doargs():
  ''' Handle args '''
  global prjname, version, finp

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hv:f:",
        ["help", "version=", "input=="])
  except getopt.GetoptError, err:
    # print help information and exit
    print str(err) # will print something like "option -a not recognized"
    usage()

  fninp = None # unset the default file name
  for o, a in opts:
    if o in ("-v", "--version"):
      version = int(a)
    elif o in ("-f", "--input"):
      fninp = a
    elif o in ("-h", "--help"):
      usage()

  if len(args) > 0:
    prjname = args[0]

  # process the project name
  if prjname.startswith("md"):  # remove the prefix `md'
    prjname = prjname[2:]
  i = prjname.find(".") # in case the project name is given as a file name
  if i >= 0:
    prjname = prjname[:i]

  # construct the input file name
  if fninp == None:
    fninp = glob.glob(prjname + ".[hc]")
    if len(fninp):
      fninp = fninp[0]
    else:
      fninp = None

  # read the input file as template, otherwise use the default
  if fninp != None:
    finp = open(fninp).readlines()

  print "project %s, input file %s, version %s" % (
      prjname, fninp, version)


def main():
  ''' main function '''
  doargs()

  obj = prjname
  hdrs = {}

  co = CCOBJ(obj)

  # add basic code
  basic = '''
#ifdef hAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <signal.h>
#include "typedefs.h"
#include "physics.h"
#include "network.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifndef RETSIGTYPE
#define RETSIGTYPE void
#endif
'''.splitlines(True)

  ca = CCGMX(basic, obj, hdrs)
  ca.s += ["#define GMXVERSION %s\n" % CCGMX.version, "\n", ]
  ca.shdr()

  if ca.version < 40010:
    cm = ccmdxv40
  elif ca.version < 40510:
    cm = ccmdxv45
  else:
    print "this version %s of GROMACS is not supported" % co.version
    raise Exception

  # set filenames to `None' to be safe
  bondfree_c = cm.get_bondfree_c(obj, None, hdrs)
  force_c = cm.get_force_c(obj, None, hdrs)
  simutil_c = cm.get_simutil_c(obj, None, hdrs)

  md_c = cm.get_md_c(obj, None, hdrs)
  runner_c = cm.get_runner_c(obj, None, hdrs)
  mdrun_c = cm.get_mdrun_c(obj, None, hdrs)


  # template replacement
  dic = {
    "/* $OBJ_DECL    */" : co.getdecl(),
    "/* $OBJ_FUNCS   */" : co.getfuncs(),
    "/* $BONDFREE_C  */" : bondfree_c.s,
    "/* $FORCE_C     */" : force_c.s,
    "/* $SIMUTIL_C   */" : simutil_c.s,
    "/* $MD_C        */" : md_c.s,
    "/* $RUNNER_C    */" : runner_c.s,
    "/* $MDRUN_C     */" : mdrun_c.s,
    }
  code = CC.tagrepl0(finp, dic)
  # add a version line
  code = ca.s + code

  # rewrite the output
  CCGMX.save0("md%s.c" % prjname, code)


if __name__ == "__main__":
  # try to accelerate the process
  try:
    import psyco
    psyco.full()
  except ImportError: pass

  main()


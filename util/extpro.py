#!/usr/bin/env python

'''
extract the protein structure in pdb format
in a frame (for Gromacs 4.0 only)
'''

import re, os, sys, glob, getopt, shutil, subprocess, math, random

verbose = 0

def extract_frame(prjname, id, tm0):
  ''' extract a pdb file from id and tm0 '''
  if not id.startswith("data"): id = "data" + id
  print "extracting %s %s" % (id, tm0)

  mdp = glob.glob("md*.mdp")[0]
  trrdt = get_trrdt(mdp)

  tpr = glob.glob("*.tpr")[0]  # search for .tpr topology file
  trr = glob.glob(id+"/*.trr")[0]
  pdb = prjname + ".pdb"
  tm = round(tm0 / trrdt) * trrdt
  print "time %g --> %g" % (tm0, tm)
  t0 = tm - trrdt*.5
  t1 = tm + trrdt*.5
  cmd = 'echo "1" | trjconv -pbc whole -s %s -f %s -o tmp.pdb -b %s -e %s' % (
      tpr, trr, t0, t1)
  runcmd(cmd, system = 1)
  # fit to the structure
  cmd = 'echo "1 1" | trjconv -s %s -f tmp.pdb -o %s -fit rot+trans' % (tpr, pdb)
  runcmd(cmd, system = 1)

  # rm
  runcmd("rm tmp.pdb", system = 1)

  # de-amberize
  cmd = "amberize -D %s" % pdb
  runcmd(cmd, system = 1)
  # remove the resulting backup
  cmd = "rm %s.bak" % pdb
  runcmd(cmd, system = 1)

  curdir = os.path.abspath(os.curdir)
  top = "topol.top"
  print "file %s, %s, %s, %s" % (curdir, pdb, top, mdp)
  return curdir, pdb, "topol.top", mdp

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
      dt = s[i+1:].strip()
      dt = dt.split(";")[0].strip()
      dt = float(dt)
  print "%s: dt = %g, nstxout = %d, trrdt = %g" % (mdp, dt, xout, dt*xout)
  return dt*xout

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

  global verbose

  for o, a in opts:
    if o in ("-v", "--verbose",):
      verbose = 1
    elif o in ("-h", "--help"):
      usage()

  if len(args) < 3:
    usage()
  return args[0], args[1], args[2]


if __name__ == "__main__":
  prjname, id, tm = doargs()
  srcdir, gro, top, mdp = extract_frame(prjname, id, float(tm))


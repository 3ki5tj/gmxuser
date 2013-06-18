#!/usr/bin/env python

"""
create different versions md1.c, md2.c, ... from the template md_t.c
  mb_t.h        -->   mb1.h, mb2.h
  md_t.c        -->   md1.c, md2.c
  md_tcore.h    -->   md1core.h, md2core.h
  md_tutil.h    -->   md1util.h, md2util.h
  zcom.h        -->   zcom1.h, zcom2.h
"""

import sys, os, getopt
import genver as gv
import zcompick

# global variables, accessible by all functions
rmlegacy = True
verbose = 0
nver = 2

class Module:
  included = 0
  def __init__(self, fn, name, opt, pfx):
    self.fn = fn
    self.name = name
    self.opt = opt
    self.pfx = pfx

# files needs to be updated
modules = [
    Module("zcom.h",       "zcom",    "z", "ZCOM"),
    Module("mb_t.h",       "mb",      "b", "MB"),
    Module("md_t.c",       "md",      "d", "AT"),
    Module("md_tcore.h",   "mdcore",  "c", "AT"),
    Module("md_tutil.h",   "mdutil",  "u", "AT"),
    Module("md2spb.h",     "spb",     "s", "XXX"),
    Module("md2bb.h",      "bb",      "r", "XXX"),
    ]

def usage():
  """
  print usage and die
  """
  print sys.argv[0], "[Options]"
  print "Options:"
  print " -a:  --all,     touch all modules"
  for i in range(len(modules)):
    print (" -%s:  --%-8s %s" %
      (modules[i].opt, modules[i].name+",", modules[i].fn))
  print " -v:  --verbose, verbose"
  print " -h:  --help,    help"
  exit(1)

def doargs():
  ''' Handle arguments from command line options '''

  global verbose, modules

  opshort = "av:h"
  oplong = ["verbose=", "all", "help"]
  for i in range(len(modules)):
    opshort +=   modules[i].opt;
    oplong  += [ modules[i].name ];

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], opshort, oplong);
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-a", "--all"):
      for i in range(len(modules)):
        modules[i].included = 1
    elif o in ("-h", "--help"):
      usage()
    else:
      for i in range(len(modules)):
        if "-" + modules[i].opt == o or "--" + modules[i].name == o:
          modules[i].included =  1
          break

  for i in range(len(modules)):
    if modules[i].included:
      break
  else:
    usage()

  for i in range(len(modules)):
    print "%s: %d;" % (modules[i].name, modules[i].included),
  print ""

def call_og(file):
  pt = os.path.splitext(file)
  proto = pt[0] + ".0" + pt[1]
  if os.path.exists(proto):
    #print "run object generator..."
    ret = os.system("python og.py %s" % file)
    if ret != 0:
      print "error occurs in og"

def apply_genver(fninp, title):
  gv.fninp    = fninp
  gv.title    = title
  gv.verbose  = verbose
  for ver in range(1, nver + 1):
    gv.genver(ver)

def dozcom(fn):
  ''' update zcom '''
  keys = [ [] ]*3
  keys[1] = ['def', 'util', 'rng', 'opt', 'cfg', 'log', 'endn', 'ss', 'argopt']
  keys[2] = keys[1] + ['rc', 'rv3', 'mds', 'eig', 'clus', 'specfunc', 'distr', 'hist']

  for i in [1, 2]:
    pt = os.path.splitext(fn)
    out = "%s%s%s" % (pt[0], i, pt[1])
    zcompick.mksmall(fn, out, keys[i])

def main():
  doargs()

  for i in range(len(modules)):
    if not modules[i].included: continue
    if modules[i].name == "zcom": # generate zcom1.h, zcom2.h
      dozcom(modules[i].fn)
    elif modules[i].name == "spb":
      call_og("md2spb.h")
    elif modules[i].name == "bb":
      call_og("md2bb.h")
    else:
      call_og(modules[i].fn)
      apply_genver(modules[i].fn, modules[i].pfx)

if __name__ == "__main__":
  main()


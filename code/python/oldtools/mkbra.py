#!/usr/bin/env python
''' make folded brackets '''

import os, sys, glob, getopt
from math import *

fnin = "RMSBB"
cutoff = 0.6
dtmin = 0.0
fnbracket = "fold.bra"
verbose = 1

def usage():
  ''' print usage and die '''
  print "%s yourfile" % sys.argv[0]
  print '''  -o: output file
  -c: cutoff (nm)
  -t: dtmin (ps)
  -v; verbose
  -h: help'''
  exit(1)

def doargs():
  global fnin, cutoff, verbose, dtmin, fnbracket

  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hc:vt:o:",
         ["help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  for o, a in opts:
    if o == "-v":
      verbose = 1
    elif o in ("-c",):
      cutoff = float(a)
    elif o in ("-t",):
      dtmin = float(a)
    elif o in ("-o",):
      fnbracket = a
    elif o in ("-h", "--help"):
      usage()

  if len(args) > 0:
    fnin  = args[0]

def main():
  doargs()
  state = 0
  str = ""
  cnt = 0
  lines = open(fnin).readlines()
  n = len(lines)
  for i in range(n):
    s = lines[i]
    arr = s.strip().split()
    if len(arr) < 2:
      print "not enough columns, s = %s" % s
      raise Exception
    t = float(arr[0])
    rmsd = float(arr[1])
    if state == 0:
      if rmsd < cutoff: # start bracket
        bra = floor(t)
        state = 1
    elif state == 1:
      if rmsd > cutoff or i == n-1:
        ket = ceil(t)
        state = 0
        if ket - bra > dtmin:
          str += "%.0f %.0f\n" % (bra, ket)
          cnt += 1
  open(fnbracket, "w").write(str)
  print "%d brackets, written to %s" % (cnt, fnbracket)

if __name__ == "__main__":
  main()

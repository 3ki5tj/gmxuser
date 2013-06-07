#!/usr/bin/env python
''' average RMSBB and others '''

import os, sys, glob
from math import *

fnin = "RMSBB"

def main():
  global fnin

  if len(sys.argv) > 1:
    fnin = sys.argv[1]
  else:
    usage()

  sum = 0.
  xsm = 0.
  xsm2 = 0.
  xmin = 1e10
  imin = -1
  xmax = -xmin
  imax = -1
  i = 1
  for s in open(fnin).readlines():
    arr = s.strip().split()
    if len(arr) < 2:
      print "not enough columns, s = %s" % s
      raise Exception
    x = float(arr[1])
    sum += 1.
    xsm += x
    xsm2 += x*x
    if x > xmax:
      xmax = x
      imax = i
    if x < xmin:
      xmin = x
      imin = i
    i += 1
    t = arr[0]
  xav = xsm/sum;
  var = xsm2/sum - xav*xav;
  print "x = %g +/- %g, samples = %g, t = %s" % (xav, sqrt(var), sum, t)
  print "xmin = %g, %d " % (xmin, imin),
  print "xmax = %g, %d " % (xmax, imax)

def usage():
  ''' print usage and die '''
  print "%s yourfile" % sys.argv[0]
  exit(1)

if __name__ == "__main__":
  main()

#!/usr/bin/env python

import sys,os,shlex,subprocess,shutil
from pyzcom import run_cmd

target = "../../../md2.tar.gz"
filelist = ['md1_README', 'md2_README',
    'md2.c', 'md2util.h', 'md2core.h', 'md2bb.h', 'md2spb.h',
    'mb2.h', 'zcom2.h', 'gmx_unused.h',
    'md2conv.c', 'md2spbx.h',
    'clusmb.c', 'clusmc.h', 'zeutil.h', 'ztoputil.h',
    'vdist.c', 'dihmb2.c', 'mkspx.c',
    'Makefile.am', 'Makefile.in',
    'CMakeLists.txt']
# files in the ala12 test folder
alatests  = ['at.cfg', 'init.gro', 'md2.mdp', 'md2.tpr', 'README', 'topol.top',
    'md2pre.mdp', 'pre.bin', 'pre.txt', 'spbfit.log',
    'test', 'clean',]
filelist += ['md2_test/a3d/'+f for f in alatests]

print "The full file list:", filelist

# kill dead entries
def prune_filelist(list):
  newlist=[]
  for f in list:
    if os.path.exists(f):
      newlist += [f]
  print "I will include these files:", newlist
  return newlist

def create_phony(files):
  for fnm in files:
    fnmbak=fnm+'b'
    phony=os.path.join('oldver', fnm+'4')

    print "Backing up", fnm, "to",fnmbak
    os.rename(fnm, fnmbak)
    print "Using the PHONY", phony, "instead of",fnm
    shutil.copy2(phony, fnm)

def restore_orig(files):
  for fnm in files:
    fnmbak=fnm+'b'
    print "Restoring the original "+fnm
    os.rename(fnmbak, fnm)

devver=False
# this is the folder where the script resides, not where we run it
# mypath=sys.path[0]
runpath=os.getcwd()
#print mypath, mypath.find("gromacs-4.0."),os.getcwd()
if (runpath.find("gmxgit") >= 0 or
    runpath.find("gromacs-4.0.") >=0): devver=True

if not devver:
  print "This is not a testing system, do you really want to continue?"
  if raw_input().startswith("n"): exit(1);
  create_phony(['Makefile.am', 'Makefile.in'])

run_cmd(['tar', 'cvphzf', target]+prune_filelist(filelist))

if not devver:
  restore_orig(['Makefile.am', 'Makefile.in'])

print "%s size %d" %(target, len(open(target).read()))

run_cmd(['gnome-open',target])


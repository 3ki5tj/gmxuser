#!/usr/bin/env python

import sys,os,shlex,subprocess,shutil
from zcom import runcmd

target = "../../../md1.tar.gz"
filelist = ['md1_README', 'md1_CHANGES', 'md1.c', 'md1util.h', 'md1core.h',
    'mb1.h', 'zcom1.h', 'gmx_unused.h',
    'mhex.c', 'zeutil.h', 'mhex_README',
    'mhidx.c',
    'Makefile.am', 'Makefile.in',
    'CMakeLists.txt']

# files in the ala12 test folder
alatests  = ['at.cfg', 'init.gro', 'md1.mdp', 'md1.tpr', 'README', 'topol.top']
filelist += ['md1_test/ala12/'+f for f in alatests]

if os.getcwd().find("gromacs-4.5") >= 0:
  # files in the GB1 test folder
  gbtests  = ['at.cfg', '2GB1.pdb', 'init.gro', 'run.mdp', 'run.tpr', 'README', 'topol.top']
  filelist += ['md1_test/2GB1/'+f for f in gbtests]

runpath = os.getcwd()
#print mypath, mypath.find("gromacs-4.0."),os.getcwd()
if runpath.find("gmxfoo") >= 0: devver = 1

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

devver = 0
# this is the folder where the script resides, not where we run it
# mypath=sys.path[0]
runpath=os.getcwd()
#print mypath, mypath.find("gromacs-4.0."),os.getcwd()
if runpath.find("gmxfoo") >= 0: devver = 1

if not devver:
  print "This is not a testing system, do you really want to continue?"
  if raw_input().startswith("n"): exit(1);
  create_phony(['Makefile.am', 'Makefile.in'])

runcmd(['tar', 'cvphzf', target]+prune_filelist(filelist))

if not devver:
  restore_orig(['Makefile.am', 'Makefile.in'])

runcmd(['gnome-open',target])


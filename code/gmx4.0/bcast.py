#!/usr/bin/env python

'''
this is a python clone of the original shell script
can remotely compile and build
'''

import sys, os, shlex, subprocess

prog = "bcast.py"

def runcmd(input, capture = False, system = False):
  if capture: pipe = subprocess.PIPE
  else: pipe = None

  # detect if the input is a string or not
  if type(input) == str:
    cmdstr = input
    cmd = cmdstr.split()
  else:
    cmd = input
    cmdstr = ''.join([s+' '  for s in cmd])

  print "CMD:", cmdstr
  if system:
    retcode = os.system(cmdstr)
    oe = ["", ""]
  else:
    p = subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe = p.communicate()
    retcode = p.returncode
  return (retcode, oe[0], oe[1])


def buildgmx(addr):
  ''' build gromacs on remote computers through ssh command '''
  # search the mkgmx command in bash script
  cmd = ['ssh', addr, 'grep', '"alias mkgmx="', '~/.bashrc']
  # capture the output from grep
  out = runcmd(cmd, capture=1)
  if out[0] != 0: return

  # extract the part between two quotes
  i0 = out[1].find("'")
  i1 = out[1].rfind("'")
  mkgmx_cmd = out[1][i0+1:i1]

  # parse the mkgmx command
  # and attach to the ssh command
  cmd = ['ssh', addr] + mkgmx_cmd.split()
  runcmd(cmd)

def help():
  ''' print usage and die '''
  print ("Usage:\n"
      + "  " + prog + " TARGET [ACTION]\n\n"
      + "TARGET could be tigger,sugar,stic,biou,tacc,archive,...\n"
      + "ACTION could be put (default), get, all (tigger only), make (not for tigger)\n");
  exit(1)

def main():
  global prog
  prog = sys.argv[0]

  argc = len(sys.argv)
  if argc <= 1: help()

  # determine target and action
  target = sys.argv[1]
  if target.find("help") >= 0: help();
  action = sys.argv[2] if argc >= 3 else "put"

  # check the validity target
  gmxname="gromacs"
  if target == "tigger" :
    addr = "czhang@tigger.rice.edu"
  elif target == "tacc":
    addr = "tg458475@tg-login.ranger.tacc.teragrid.org"
  elif target == "cave" :
    addr = "czhang@10.12.18.201"
  elif target == "archive" :
    addr = "../archive"
  elif target == "gpu" :
    addr = "cz1@badlands.rcsg.rice.edu"
  elif target.startswith("biou") or target in ("sugar", "stic", "ada"):
    addr = "cz1@"+target+".rice.edu"
  else:
    print "unknown target %s" % target
    exit(1)

  if target == "archive":
    dest = addr + "/" + gmxname + "4"
  else:
    dest = addr + ":" + gmxname

  if action == "get":
    afrom = dest + "/"
    ato = "./"
  else:
    afrom = "./"
    ato   = dest + "/"

  cmd = ('rsync -ravz --copy-links --exclude=".*" --exclude="build*" --exclude="*~" '
      +'--exclude="*.pyc" --exclude="*.bak" --exclude="*.spr" --exclude="*.orig" '
      +'--exclude="*.def[0-9]" --exclude="*.bak*" --exclude="tmp*" --exclude="old" '
      +'--exclude="autom4te.cache" --exclude="test" --exclude="oldver" --exclude="dev" '
      + afrom + '* ' + ato)

  # this rsync command can only be used by calling system
  runcmd(cmd, system = 1)

  if target == "tigger":
    if action == "all": # bcast.py tigger all
      '''
      since the target is tigger,
      the caller of the script is still the local computer
      we use ssh to send command from local to tigger
      to ask it to
      1. switch to its ~/gromacs directory
      2. run this script like
          bcast.py sugar make
      '''
      cmd = ["ssh", addr,
             'cd', 'gromacs;',
             prog, 'sugar', 'make;',
             prog, 'stic',  'make;',
             prog, 'biou',  'make'];
      runcmd(cmd)
  else:
    '''
    since the target is no longer tigger
    tigger is the caller, now tigger is running
      bcast.py sugar make
    since the action is make, tigger tries to
    send ssh command to e.g. sugar
      and execute mkgmx there
    '''
    if action == "make": # bcast.py xxx make
      buildgmx(addr)

if __name__ == "__main__":
  main()


#!/usr/bin/env python

'''
this is a python clone of the original shell script
can remotely compile and build
'''

import sys
import os, shlex
import subprocess

def run_cmd(input, capture=False, old_fashion=False):
  if capture:
    pipe=subprocess.PIPE
  else:
    pipe=None

  # detect if the input is a string or not
  if type(input)==type(""):
    cmdstr=input
    cmd=cmdstr.split()
  else:
    cmd=input
    cmdstr=''.join([s+' '  for s in cmd])

  print "CMD:",cmdstr
  if old_fashion:
    retcode=os.system(cmdstr)
    oe=["",""]
  else:
    p=subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe=p.communicate()
    retcode=p.returncode

  return (retcode, oe[0], oe[1])


def buildgmx(addr):
  '''
  optionally build gromacs on remote computers
  through ssh command
  '''
  # search the mkgmx command in bash script
  cmd=['ssh', addr, 'grep', '"alias mkgmx45="', '~/.bashrc']
  # capture the output from grep
  out=run_cmd(cmd, capture=True)
  if out[0] != 0: return

  # extract the part between two quotes
  i0 = out[1].find("'")
  i1 = out[1].rfind("'")
  mkgmx_cmd = out[1][i0+1:i1]

  # parse the mkgmx command
  # and attach to the ssh command
  cmd = ['ssh', addr] + mkgmx_cmd.split()
  run_cmd(cmd)


def help_n_quit():
  '''
  print usage and die
  '''
  print ("Usage:\n"
      +"  "+myname+" TARGET [ACTION]\n\n"
      +"TARGET could be tigger,sugar,stic,biou,archive,...\n"
      +"ACTION could be put(default), get, all(tigger only), make (not for tigger)\n");
  exit(1)


def main():
  myname=sys.argv[0]

  argc = len(sys.argv)

  # determine target and action
  target = sys.argv[1] if argc > 1 else ""
  if argc <= 1 or target=="help":
    help_n_quit();
  action = sys.argv[2] if argc >= 3 else "put"

  # check the validity target
  gmxname="gmx45/gromacs"
  if target == "tigger" :
    addr="czhang@tigger.rice.edu"
  elif target == "cave" :
    addr="czhang@10.12.18.201"
  elif target == "archive" :
    addr="../archive"
  elif target == "gpu" :
    addr="cz1@badlands.rcsg.rice.edu"
  elif target.startswith("biou") or target in ("sugar", "stic", "ada"):
    addr="cz1@"+target+".rice.edu"
  else:
    print "unknown target $1"
    exit(1)

  if target=="archive":
    dest=addr+"/"+gmxname+"4"
  else:
    dest=addr+":"+gmxname

  if action == "get":
    afrom = dest+"/"
    ato = "./"
  else:
    afrom = "./"
    ato   = dest+"/"

  cmd=('rsync -ravzL --exclude=".*" --exclude="build*" --exclude="*~" '
      +'--exclude="*.pyc" --exclude="*.bak" --exclude="*.spr" --exclude="*.orig" '
      +'--exclude="*.1" --exclude="*.bak*" '
      +'--exclude="autom4te.cache" --exclude="test" --exclude="oldver" --exclude="dev" '
      + afrom + '* ' + ato)

  # this rsync command can only be used in old fashion
  run_cmd(cmd, old_fashion=True)

  if target == "tigger":
    if action == "all":
      '''
      since the target is tigger,
      the caller of the script is still the local computer
      we use ssh to send command from local to tigger
      to ask it to
      1. switch to its ~/gmx45/gromacs directory
      2. run this script like
          bcast.py sugar make
      '''
      cmd =["ssh", addr];
      cmd+=['cd', 'gmx45/gromacs;',
          myname, 'sugar', 'make;',
          myname, 'stic',  'make;',
          myname, 'biou',  'make'];

      run_cmd(cmd)
  else:
    '''
    since the target is no longer tigger
    tigger is the caller, now tigger is running
      bcast.py sugar make
    since the action is make, tigger tries to
    send ssh command to e.g. sugar
      and execute mkgmx there
    '''
    if action == "make":
      buildgmx(addr)


if __name__ == "__main__":
  main()


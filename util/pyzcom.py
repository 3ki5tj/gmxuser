'''
common python  modules
'''

import os,sys,shutil,getopt,re,subprocess

def run_cmd(input, capture = 0, old_fashion = 0, verbose = 1):
  '''
  run a system command
  and optionally capture standard output/error
  return tuple of return code, stdout, stderr
  the latter two are empty if capture is not set
  to use os.system() instead of subprocess,
  set old_fashion to 1
  if `verbose' is set, the command is echoed before executed
  '''
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

  if verbose >= 1:
    print "CMD:",cmdstr
    if verbose >= 2:
      raw_input("proceed?")

  if old_fashion:
    retcode=os.system(cmdstr)
    oe=["", ""]  # no stdout or stderr
  else:
    p=subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe=p.communicate()
    retcode=p.returncode

  return (retcode, oe[0], oe[1])


def backup_file(file, ext = "", verbose = 1):
  '''
  backup file to a nonexisting name,
  lead by the original file name `file'
  with optional extension `ext', e.g., ".bak"
  if `verbose' is >= 2, confirmation is needed after the backup
  '''
  fn = file + ext
  i = 1
  if ext == "": fn += str(i)
  while os.path.exists(fn):
    fn = file + ext + str(i)
    i += 1
  shutil.copy2(file, fn)
  if verbose > 0:
    print file, "is backed up to", fn
    if verbose >= 2:
      raw_input("press Enter to proceed ...")


#!/usr/bin/env python

'''
create a specific version, e.g. abc1.h, abc2.h, from the template abc_t.h
this file is used by updmd.py

To use it,
1. write a module source code, and name it as xxx_t.yyy
2. structure the code with MACROS
    #if AAA_VER == 0
      ...
    #elif AAA_VER > 2
      ...
    #else
      ...
    #endif
3. setup input as "abc_t.h", title as "AAA", call genver(ver) with
   the version ver.
'''

import os,sys,shutil,getopt,re,filecmp

# module attributes
input     = ""
title     = ""
rmlegacy  = True
verbose   = 2
nver      = 2
backup    = False


# the version pattern like
#   (strver >= 2)
version_pattern = "(\!?)\s*\(*\s*(%s)"  + r"\s*([=<>!]=|[<>])\s*" + "([0-9]+)\s*\)*"

def redver(cond, strver, ver):
  '''
  reduce a condition that contains the version_pattern
  '''
  pat = version_pattern % strver
  newpat = ".*" + pat + ".*"
  m = re.search(pat, cond)
  if not m:
    print "no pattern is found, cannot simplify"
    raise SyntaxError

  # extract the version part
  vexpr = m.group(0)
  # translate to python syntax
  vexpr1 = "(%s %s %s)" % (ver, m.group(3), m.group(4))
  if m.group(1) == "!": vexpr1 = "not "+vexpr1

  try:
    val = eval(vexpr1)
  except SyntaxError:
    print "fuck [%s] --> [%s], ver = %s" % (vexpr, vexpr1, ver)
    raise SyntaxError

  if val not in (True, False):
    print "Error: [%s] --> %s" % (vexpr1, val)
    return
  val = int(val)

  # construct a new condition
  m0 = m.start(0)
  m1 = m.end(0)
  condnew = (cond[:m0] + str(val) + cond[m1:]).strip()

  p01 = r"\(*\s*[01]\s*\)*"
  if re.match(p01+"$", condnew):
    return val  # simple expression

  #print "compound expression\n[%s]\n=>\n[%s]" % (cond, condnew)

  orpat0  = p01 + r"\s*\|\|(.*)"
  orpat1  = r"(.*)\|\|\s*" + p01 + "$"
  andpat0 = p01 + r"\s*\&\&(.*)"
  andpat1 = r"(.*)\&\&\s*" + p01 + "$"

  c0 = getcond(condnew, orpat0)
  if c0:
    if val: return 1
    else: return c0.strip()
  c0 = getcond(condnew, orpat1)
  if c0:
    if val: return 1
    else: return c0.strip()
  c0 = getcond(condnew, andpat0)
  if c0:
    if val: return c0.strip()
    else: return 0
  c0 = getcond(condnew, andpat1)
  if c0:
    if val: return c0.strip()
    else: return 0

  print "don't know how to handle condnew: [%s]" % condnew
  raw_input()


def getcond(cond, pat):
  m = re.match(pat, cond)
  if m: return m.group(1);
  else: return None

def tpl2ver(fn, ver, pattern = r"_t"):
  '''
  replace pattern in file fn by ver
  '''
  if re.search(pattern, fn):
    return re.sub(pattern, repr(ver), fn)
  else:
    print ("WARNING: no pattern %s is found in %s"
        % (pattern, fn))
    # change abc.ext to abc1.ext
    pt = os.path.splitext(fn)
    return pt[0] + repr(ver) + pt[1]

def genver(ver):
  """
  create a specific version, e.g. mb2.h, from the template mb_t.h
    * input file is given by input
    * generate a special version, such that the program thinks
      `strver' is defined as `ver'
      limited to one level nesting
    * eliminate legacy code, if rmlegacy is specified
    * change zcom.h to zcom1.h, etc
    * change xxx_txxx.h to xxx1xxx.h in #include
  """

  plevel        =  0
  lineno        =  1
  ignore_at     = -1            # the plevel of the if
  ignore_cnt    =  0
  lnumout       =  1            # number of lines written to file
  last_nonblank =  0            # last non-blank line written to file

  cond_at       = -1            #  conditional starts
  cond_cnt      =  0            #  number of conditional
  cond_on       =  True
  skip_line     =  False

  global title

  if title == "":
    print "cannot proceed, neither title nor strver is specified"
    return

  if title[len(title) - 1] != '_': # add an underscore, if it is missing
    title += '_'
  strver = title + "VER"
  strlegacy = title + "LEGACY"
  verpat = version_pattern % strver

  # check version number
  if not ver in range(10):
    print "version:", ver, "is not supported yet"
    return

  # create a name for automatic file name
  output = tpl2ver(input, ver)
  # open the output file ready to write
  output_tmp = output + ".tmp"
  fo = open(output_tmp, 'w')

  # read the input file line by line
  for line in open(input, 'r'):
    lin = line.lstrip()
    skip_line = False

    # search for #include "zcom.h"
    # replace by zcom1.h for version 1
    pattern =  r'(\s*#include\s*"zcom)(\.h"\s*)';
    if re.search(pattern, line):
      i = line.find('.')
      line = line[:i] + repr(ver) + line[i:]
      if verbose > 0:
        print "Replace zcom.h:", line

    pattern = r'\s*#include\s*"\w*_t\w*.\w*"\s*'
    if re.search(pattern, line):
      line1 = re.sub("_t", repr(ver), line)
      if verbose >= 1:
        print "\t", line, "is changed to\n\t", line1
      line = line1

    '''
    increase the level if a preprocessor conditional is met
    this includes three cases #if, #ifdef, #ifndef
    '''
    if lin.startswith("#if"):
      plevel+=1
      if verbose >= 3:
        print "LEVEL +",plevel,"lineno=",lineno, ":", lin

      # test if it belongs to an ignored group
      # that is, legacy code
      # we do not support #else _LEGACY yet
      if ignore_at < 0:
        if (rmlegacy and
            lin.startswith("#ifdef") and
            lin.find(strlegacy) >= 0):
          ignore_at = plevel
          if verbose >= 1:
            print ("[%2d: %2d] %6d: start ignoring due to #ifdef %s,"
                % (ignore_cnt+1, plevel, lineno, strlegacy))

        elif (rmlegacy and
              lin.startswith("#if ") and
              # \s means space
              re.search("[^!]defined\(\s*" + strlegacy + "\s*\)", lin)):
          ignore_at = plevel
          if verbose >= 1:
            print ("[%2d: %2d] %6d: start ignoring due to #if defined(%s),"
                % (ignore_cnt+1, plevel, lineno, strlegacy))
        # just started ignoring
        if ignore_at >= 0:
          ignore_cnt += 1
          if verbose >= 2:
            print "%12s%s" % ("", line),

      # test if a conditional lead by #if is true or not
      # we only care if it's not inside an ignored block
      if ignore_at < 0:
        cond = lin[3:].rstrip()
        if ( lin.find(strver) >= 0 and lin.startswith("#if ")
            and re.search(verpat, cond) ):
          '''
          generally, we should not assume cond_at < 0, because we can have
          something nested like
            #if strver >= 2
            #if strver == 3
            #endif
            #else
            #endif
          but for the moment, we just assume one level
          '''
          if cond_at < 0:
            cond_on = redver(cond, strver, ver)
            if cond_on not in (0, 1): # complex expression
              line1 = line[:3] + " " + cond_on + '\n'
              if verbose > 0:
                print "\t%s\n-->\n\t%s\n" % (line.strip(), line1)
                #raw_input()
              line = line1
            else:
              cond_at   = plevel
              skip_line = True
              if verbose >= 1:
                print ("[%2d: %2d] %6d: start conditional due to #if%s is %r,"
                    % (cond_cnt+1, cond_at, lineno, cond, cond_on))
              if verbose >= 3: raw_input()

        # just started conditionals
        if cond_at >= 0:
          cond_cnt += 1
          if verbose >= 2:
            print "%12s%s" % ("", line),

    elif lin.startswith("#else"):
      if cond_at >= 0 and cond_at == plevel:
        cond_on = not cond_on
        skip_line = True
        if verbose >= 1:
          print ("[%2d: %2d] %6d: switch conditional due to #else to %r,"
             % (cond_cnt, cond_at, lineno, cond_on))
        if verbose >= 2:
          print "%12s%s" % ("", line),
        if verbose >= 3: raw_input()

    # note, #elif can act like an #else or act like an #if
    elif lin.startswith("#elif"):
      has_cond = 0
      # note: we want to match the current existing level
      if lin.find(strver) != 0 and plevel == cond_at:
        cond = lin[5:].rstrip()
        if re.search(verpat, cond):
          has_cond = 1

      if has_cond:
        cond_on = redver(cond, strver, ver)
        if cond_on not in (0, 1): # complex expression
          line1 = line[:5] + " " + cond_on + '\n'
          if (verbose > 0):
            print "\t%s\n-->\n\t%s\n" % (line.strip(), line1)
          line = line1
        else:
          cond_at   = plevel
          skip_line = True
          #print "a match is found!, pattern:", pattern, "  cond:", cond
          if verbose >= 1:
            print ("[%2d: %2d] %6d: start conditional due to #if%s is %r,"
                % (cond_cnt+1, cond_at, lineno, cond, cond_on))
          if verbose >= 3: raw_input()

      elif cond_at >= 0 and cond_at == plevel:
        ''' #elif act like #else '''
        cond_on = not  cond_on
        skip_line = True
        if verbose >= 1:
          print ("[%2d: %2d] %6d: switch conditional due to #elif%s to %r,"
                % (cond_cnt, cond_at, lineno, cond, cond_on))
        if verbose >= 2:
          print "%12s%s" % ("", line),
        if verbose >= 3: raw_input()

    # note: there is a second part for handling #endif
    elif lin.startswith("#endif"):
      if cond_at >= 0 and cond_at == plevel:
        skip_line = True
        if verbose >= 3: raw_input()


    # output string
    writeit = False
    if (ignore_at < 0 and
        not skip_line and
        cond_on):
      if lin != "":
        last_nonblank = lnumout
      # we allow two successive blank lines, but no more
      if lnumout - last_nonblank <= 2:
        fo.write(line)  # write output file
        lnumout += 1
        writeit  = True

    if verbose >= 3 and (verbose >= 4 or writeit):
      print ("%5s: %6d, pL: %2d, ign: %2d, cond: %2d, line: %s"
          % (repr(writeit), lineno, plevel, ignore_at, cond_at, lin.rstrip()) )

    # we update plevel change due #endif here
    if lin.startswith("#endif"):
      if ignore_at >= 0 and ignore_at == plevel:
        # stop ignoring
        if verbose >= 1:
          print ("[%2d: %2d] %6d: stop ignoring due to #endif"
              % (ignore_cnt, plevel, lineno))
          if verbose >= 2:
            print ""
        ignore_at = -1
      elif cond_at >= 0 and cond_at == plevel:
        if verbose >= 1:
          print ("[%2d: %2d] %6d: stop conditional due to #endif"
              % (cond_cnt, plevel, lineno))
          if verbose >= 2:
            print ""
        cond_at = -1
        cond_on = True

      # decrease depth
      plevel -= 1
      if verbose >= 3:
        print "LEVEL -", plevel, "lineno=", lineno, ":", lin
      if plevel < 0:
        print "plevel is negative at line:", lineno
        exit(1)

    # increase the line number
    lineno += 1

  # sync to the true output
  fo.close()
  if (not os.path.exists(output)
      or not filecmp.cmp(output_tmp, output)):
    if verbose > 0 : print "-" * 64
    print "Generating version %2d,  Input: %-12s  Output: %-12s" % (ver, input, output)
    if verbose > 0 : print "-" * 64
    shutil.copy(output_tmp, output)
  os.remove(output_tmp)

def usage():
  """
  print usage and die
  """
  print sys.argv[0], "[Options]"
  print "Options:"
  print " --input=, or -i"
  print " --verbose=[0-9], -v"
  print " --remove-legacy, --keep-legacy"
  exit(1)

def doargs():
  '''
  Handle common parameters from command line options
  results saved to module attributes
  '''

  global input, title, nver, rmlegacy, verbose

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hv:i:t:n:",
         ["help", "input=", "verbose=", "title=", "nver=",
           "keep-legacy", "remove-legacy"])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("--keep-legacy"):
      rmlegacy = False;
    elif o in ("--remove-legacy"):
      rmlegacy = True;
    elif o in ("-n", "--nver"):
      nver = a
    elif o in ("-t", "--title"):
      title = a
    elif o in ("-i", "--input"):
      input = a
    elif o in ("-h", "--help"):
      usage()

  if (title == "" or input == ""):
    print "missing title or input file"
    usage()


def main():
  doargs()

  for ver in range(1, nver + 1):
    genver(input, "", title, ver, rmlegacy, verbose)

if __name__ == "__main__":
  main()



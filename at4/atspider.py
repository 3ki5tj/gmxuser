#!/usr/bin/env python

'''
build a self-contained GROMACS engine (deprecated)
* make mdx.c from md.c + runner.c + mdrun.c, see mkmdx_c() (key function)
* make mdxutil.c from force.c + sim_util.c, see mkmdxutil_h()
especially useful for GROMACS 4.5
'''

import re, getopt, os, sys

tgtver = (1, 2)  # do atver from 1 to 2, that is md1 and md2
atver = 1
use_at = 1



def str2re(s):
  ''' translate a plain string an equivalent regular expression
      for functions that require regular expression as input '''

  # list of special characters
  special = "+-*.&?!|()[]{}"
  t = ""
  s = re.sub("\s+", " ", s)
  for i in range( len(s) ):
    if s[i].isspace(): # convert spaces to `\s+'
      c = "\s+"
    elif s[i] in special: # add a backslash `\'
      c = "\\" + s[i]
    else: c = s[i]
    t += c
  return t



def mkcodepretty(s):
  ''' make code better looking by calling cspacer.py
      `s' is a string array (lines) '''

  # remove trailing spaces first
  s = [ln.rstrip() + '\n' for ln in s]

  try: # to improve code by cspacer
    import cspacer as cs
    cs.verbose = 0
    # turn on advanced options
    cs.use_rule_add = True
    cs.use_rule_paren2 = True
    cs.use_rule_nocppcmt = True
    s, nchanges = cs.addspace0(s)
  except ImportError: pass
  return s



def savelines(fn, s):
  s = mkcodepretty(s)
  # update only if necessary
  news = "".join(s)
  olds = None
  if os.path.exists(fn):
    olds = open(fn).read()
  if olds != news:
    print "writing %s ..." % fn
    open(fn, "w").write(news)
  else:
    print "no need to update", fn



def proto2call(proto, indent = 2):
  ''' translate function definition to a typical call
      `proto' is a string of lines
      indent is the number of spaces in an indent '''

  call = proto[:]
  # first get function name
  par = call[0].find("(")
  name = call[0][:par]
  j = name.rfind(" ")
  call[0] = call[0][j+1:]
  n = len(call)
  for i in range(n):
    s = call[i]
    # j0 is the position of first parameter
    if i == 0:
      j0 = s.find("(") + 1
    else: j0 = s.index(s.lstrip()[0]) # first nonblank character

    ls = call[i][j0:].strip().rstrip(",)") # ls the parameter list
    als = ""
    for pa in ls.split(","):
      arr = pa.split()
      if len(arr) >= 2:
        als += arr[1].strip("*[]") + ", "
    if i == n-1: als = als.rstrip(" ,") + ");"
    call[i] = ' '*indent + s[:j0] + als.rstrip() + "\n"
  return call



def getgmxroot(root = None):
  ''' determine the GROMACS root directory '''

  if root and os.path.exists(root): # the user given one exists
    return root
  curdir = os.path.abspath( os.getcwd() )
  while not (curdir == "/" or curdir.endswith(":\\")):
    curdir = os.path.abspath( os.path.join(curdir, os.pardir) )
    # there is a file `AUTHORS' under the root
    if os.path.exists( os.path.join( curdir, "AUTHORS" ) ):
      break
  else:
    print "cannot determine the GROMACS root from " + os.getcwd()
    raise Exception
  return curdir



# @staticmethod
def getgmxver(gmxroot = None):
  ''' get GROMACS version from CMakeLists.txt or configure.ac
      works for GROMACS v4.0 and v4.5 '''
  gmxroot = getgmxroot(gmxroot)

  cfgac = os.path.join(gmxroot, "configure.ac")
  cmake = os.path.join(gmxroot, "CMakeLists.txt")

  if os.path.exists(cmake): # version 4.5 or later
    for s in open(cmake):
      # we look for a line that looks like
      # set(PROJECT_VERSION "4.5.6-dev")
      ln = s.strip()
      if ln.startswith('set(PROJECT_VERSION'):
        # we take the first 5 letters, i.e., 4.5.6
        ln = ln.split()[1].strip().replace('"', '').replace(')', '')[:5]
        # convert to int 40506
        version = int( ln.replace('.', '0') )
        break
    else:
      print "cannot determine version from %s" % cmake
      raise Exception
    isgmx4 = 0
    sopenmm = "GMX_OPENMM"

  elif os.path.exists(cfgac): # v4.0 or earlier
    for s in open(cfgac):
      if s.startswith("AC_INIT"):
        # a typical string is like
        #   AC_INIT(gromacs, 4.5.6-dev, [gmx-users@gromacs.org])
        # so we take the second number in the parentheses ()
        sver = s.split(",")[1].strip()
        if sver.startswith("4.0"):
          isgmx4 = 1
        else:
          print "bad gromacs version %s" % sver
          raise Exception
        version = int( sver[:5].replace(".", "0") )
        break
    else:
      print "cannot determine version number from %s" % cfgac
      raise Exception
    # v4 uses USE_OPENMM
    sopenmm = "USE_OPENMM"

  else:  # there is neither CMakeLists.txt nor configure.ac
    raise Exception

  return version, gmxroot, isgmx4, sopenmm



class CC:
  ''' Code-Changer: basic class of changing code '''

  def __init__(c, fn, hdrs = {}):
    ''' `fn' is the input file name
        `hdrs' gives a list of existing #include '''
    if fn == None:
      c.s = []
    elif type(fn) == str:
      c.s = open(fn).readlines()
    else: # treat it as input string
      raw_input( "bad input" + str(fn) )
      c.s = fn
      fn = None

    c.filename = fn
    c.begin = 0
    c.hdrs = hdrs


  def subs(c, find, repl, pat = None, doall = False, verbose = False):
    ''' search line by line for `pat', when found, replace `find' by `repl'
        `pat' is the by default the same as `find'
        if `doall', all pattern's are replaced, otherwise the first occurance '''
    if not pat: pat = find
    i = 0
    prog = re.compile(pat)
    while i < len(c.s):
      if prog.search(c.s[i]):
        if verbose:
          print "pattern %s, line %d, %s" % (pat, i, c.s[i])
        c.s[i] = re.sub(find, repl, c.s[i])
        c.begin = i
        if not doall: break
        else: continue
      i += 1


  def substt(c, find, repl, pat = None, doall = False, verbose = False):
    ''' same as `subs', but `find' and `pat' are strings '''
    if find: find = str2re(find)
    if pat: pat = str2re(pat)
    if verbose: print "find %s, pattern %s" % (find, pat)
    c.subs(find, repl, pat, doall, verbose)


  def findln(c, pat):
    ''' find a single line that contains a pattern `pat' '''
    i = 0
    prog = re.compile(pat)
    while i < len(c.s):
      if prog.search(c.s[i]):
        break
      i += 1
    else:
      print "findln cannot find [%s]" % pat
      return -1
    return i


  def findline(c, pat):
    ''' same as `findln', but `pat' is a plain string, not a regular expression '''
    return c.findln( str2re(pat) )


  def rmln(c, pat, doall = False, off0 = 0, off1 = 1, wcmt = False, verbose = False):
    '''
    remove a single line that contains a pattern `pat'
    more generally, for each line that matches the pattern `pat'
    remove the following lines from `off0' to `off1'
    if `doall', delete all occurences, otherwise stop after the first one
    if `wcmt', remove a prev. comment line, if any, that starts with `/*'
    '''
    i = 0
    prog = re.compile(pat)
    while i < len(c.s):
      if prog.search(c.s[i]):
        # the line block to be removed
        start = i + off0
        end = i + off1
        if ( wcmt and start > 0
            and c.s[start-1].strip().startswith("/*") ):
          start -= 1
        c.s = c.s[:start] + c.s[end:]
        c.begin = i = start
        if not doall: break
        else: continue
      else:
        i += 1
    else:
      if (not doall) and verbose:
        print "cannot find pattern [%s]" % pat


  def rmline(c, pat, doall = False, off0 = 0, off1 = 1, wcmt = False, verbose = False):
    ''' same as `rmln', but `pat' is a plain string, not a regular expression '''
    c.rmln( str2re(pat), doall, off0, off1, wcmt, verbose)


  def addln(c, index, line):
    ''' add a line or several lines '''
    if type(line) == str: line = [line]
    c.s = c.s[:index] + line + c.s[index:]


  def getblockend(c, i, starter = "if", ending = None, sindent = None, wsp = True, verbose = 0):
    ''' find a code block starting from line `i', the block should look like
          if (a > b &&
              b > c ) {
            ...
          }
        if the `ending' is `}' (default),
        check if the indent of the `}' matches the supposed one (determined below)
          1. if `sindent' is given, use that
          2. if `starter' is given, use whatever that is before `starter'
          3. otherwise, use the indent of line `i'
          4. if `starter' is `^', use "" (deprecated)
        return the pair (start, end)
    '''
    start = i

    # compute leading indent
    if sindent != None:
      indent = sindent
    elif starter == "^":
      indent = ""
    elif starter:
      try:
        indent = c.s[i][ : c.s[i].index(starter)]
      except:
        print "bad starter [%s] on line %s" % (starter, i)
        raise Exception
    else:
      ch = c.s[i].strip()[0]
      indent = c.s[i][ : c.s[i].index(ch)]
    if verbose:
      print "indent = %d, [%s]" % (len(indent), indent)

    # guess what's the default ending: ";"
    # or "}" if the line ends with "{"
    if not ending:
      ending = ";"
      while i < len(c.s):
        line = c.s[i].rstrip()
        if line.endswith("{"):
          ending = "}"
          break
        elif line.endswith(")"):
          # likely the if (...) ends, see if the next line starts with "{"
          if c.s[i+1].strip() == "{":
            ending = "}"
          break
        i += 1
      else:
        print "cannot find if starter from line %s (%s)\n%s" % (
            start, c.filename, c.s[start])
        raise Exception

    # search the block ending
    while i < len(c.s):
      line = c.s[i].strip()
      if ending == "}":
        # for "}" we want to make sure that
        # 1. it starts a line
        # 2. it has the right # of indents
        if line.startswith(ending):
          myindent = c.s[i][ : c.s[i].index(ending)]
          if myindent == indent:
            break
      elif line.endswith(ending): # regular ending
        break
      i += 1
    else:
      print "no ending `%s' from line %s (%s)" % (
          ending, start, c.filename)
      raise Exception

    # include the blank lines after the ending into the block
    i += 1
    if wsp:
      while i < len(c.s):
        if c.s[i].strip() != "":
          break
        i += 1
    end = i - 1
    return start, end


  def rmblk(c, pat, doall = False, starter = "if", ending = None,
      sindent = None, wsp = True, verbose = False):
    '''
    once a pattern `pat' is found (regular expression),
    remove a following code block,
    that starts with `starter' (normal string),
    ends with `}' or `;' (if `ending is not given)
    '''
    i = 0
    prog = re.compile(pat)
    while i < len(c.s):
      if prog.search(c.s[i]): # pattern found
        if verbose: print "pattern %s found in %s" % (pat, c.filename)

        # search a starting if
        start = i
        if starter != None:
          while start >= 0:
            if c.s[start].lstrip().startswith(starter):
              break
            start -= 1
          else:
            print "no starter [%s], pat [%s], line %d, file %s" % (
                starter, pat, i, c.filename)
            #open("dump.txt", "w").write(''.join( c.s ))
            raise Exception

        start, end = c.getblockend(start, starter = starter, ending = ending,
            sindent = sindent, wsp = wsp)
        if verbose:
          print "pat %s, start %s, end %s\n%s" % (
              pat, start, end, ''.join(c.s[start:end]))

        # try to include a preceeding comment line
        if start > 0 and c.s[start-1].strip().startswith("/*"):
          start -= 1

        c.s = c.s[:start] + c.s[end+1:]

        # save the block indices
        c.begin = i = start
        c.end = end
        if not doall: break
        else: continue
      i += 1


  def rmblock(c, pat, doall = False, starter = "if", ending = None,
      sindent = None, wsp = True, verbose = False):
    ''' pattern can be a regular expression '''
    c.rmblk( str2re(pat), doall, starter, ending, sindent, wsp, verbose)


  def rmfunc(c, functag, starter = None, regex = False):
    ''' convenient wrapper for rmblk '''
    if not regex: functag = str2re(functag)
    c.rmblk(functag, starter = starter, ending = "}", wsp = True)


  def findblk(c, pat, ending = "}", wsp = False, i0 = 0, i1 = -1,
      sindent = None, verbose = False):
    ''' find a code block that starts by a line with the pattern `pat'
        the wrapper for getblockend() '''
    prog = re.compile(pat)
    if i1 < 0: i1 = len(c.s)
    for i in range(i0, i1):
      if prog.search(c.s[i]): # pattern found
        if ending:
          c.begin, c.end = c.getblockend(i, starter = None, ending = ending,
              sindent = sindent, wsp = wsp, verbose = verbose)
          # ending is not so good
          #if c.begin == -1: continue
        else:
          c.begin = c.end = i

        return c.begin
    else:
      if verbose:
        print "cannot find pattern %s" % pat
      return -1



  def findblock(c, pat, ending = "}", wsp = False, i0 = 0, i1 = -1, sindent = None, verbose = 0):
    return c.findblk(str2re(pat), ending, wsp, i0, i1, sindent, verbose)


  def shdr(c):
    ''' remove redundant #include '''
    i = 0
    while i < len(c.s):
      if c.s[i].startswith("#include"):
        key = c.s[i].strip()
        if key in c.hdrs:
          start = i
          end = i+1
          # a preceding comment
          if i > 0 and re.search("/\*.*\*/", c.s[i-1]):
            start -= 1
          # a preprocessor block
          if (start > 0 and c.s[start-1].startswith("#if") and
              end < len(c.s) and c.s[end].startswith("#endif")):
            start -= 1
            end += 1
          c.s = c.s[:start] + c.s[end:]
          i = start
          continue
        else:  # add the `#include ...' as a new key
          c.hdrs[key] = True
      i += 1


  def rmdf(c, tag, ifp = "#ifdef", ifn = "#ifndef", verbose = 0):
    ''' remove an #ifdef or #if block
          #ifdef tag
            ...
          #endif
    '''
    i = 0
    while i < len(c.s):
      # we found a block starter
      if (c.s[i].startswith(ifp + " " + tag) or
         (ifn != None and c.s[i].startswith(ifn + " " + tag)) ):
        level = 1
        nelse = -1
        # found a block ending
        for j in range(i+1, len(c.s)):
          if c.s[j].startswith("#endif"):
            level -= 1
            if level == 0: break
          # increase the level
          elif c.s[j].startswith("#if"):
            level += 1
          # register the location of `else'
          elif c.s[j].startswith("#else") and level == 1:
            nelse = j
          elif c.s[j].startswith("#elif") and level == 1:
            print "cannot handle #elif\n%s"  % c.s[i:j+1]
            raise Exception
        else:
          print "cannot find #endif"
          raise Exception

        # now the `#if' block is [i : j+1]
        if verbose:
          print "preprocessor in file %s (%d : %d)" % (c.filename, i, j+1)
          print ''.join( c.s[i : j+1] )

        if c.s[i].startswith(ifp):
          if nelse >= 0:
            # if there is an `#else', we simply enable the `#else' part
            c.s = c.s[:i] + c.s[nelse+1:j] + c.s[j+1:]
          else:
            # otherwise, we remove the entire block
            c.s = c.s[:i] + c.s[j+1:]
        else:
          # this block starts with "#ifndef"
          if nelse >= 0:
            # if there is an else, we enable the "#ifndef" part
            c.s = c.s[:i] + c.s[i+1:nelse] + c.s[j+1:]
          else:
            # we remove the two lines of "#ifndef" and "#endif"
            c.s = c.s[:i] + c.s[i+1:j] + c.s[j+1:]
      else:
          i += 1


  def rmif0(c, tag = "0", verbose = False):
    ''' remove an impossible to reach code block '''
    c.rmdf(tag, ifp = "#if", ifn = None, verbose = verbose)


  def funcdef(c, funcpat, wsp = True):
    ''' return the entire function body (definition) '''
    if c.findblk(funcpat, ending = "}", wsp = wsp) < 0:
      print "cannot find function pattern %s" % funcpat
      raise Exception
    return c.s[c.begin : c.end+1]




class CCGMX(CC):
  ''' changes for GROMACS, it extends class CC (cc.py)
  each source code file (e.g., md.c, mdrun.c)
  should declare an instance via the constructor '''

  # these are class level ``static'' variables
  # so they can be set only once
  version = -1
  rootdir = None
  isgmx4 = 0
  sopenmm = "GMX_OPENMM"

  # use a dictionary to map function names
  # the mapped functions will further have a prefix
  fmap0 = {
    "mdrunner":           "runner",
    "do_md":              "domd",
    "do_force":           "doforce",
    "do_force_lowlevel":  "doforcelow",
    "calc_bonds":         "calcbonds",
  }

  obj = "foo" # common variable name for the object
  pfx = "gmxfoo" # function prefix
  varmode = "foomode" # this is an integer parameter
                      # passed from the command line
                      # better mechanism in the future?

  def __init__(c, fn, obj = None, hdrs = {}):
    # call the base class constructor
    CC.__init__(c, fn, hdrs)

    # compute the GROMACS version
    c.getversion()
    if obj != None:
      c.obj = obj
      c.pfx = c.obj + "gmx"
      c.varmode = c.obj + "mode"
      c.fmap = {}
      for key in c.fmap0: # add a prefix
        c.fmap[key] = c.pfx + "_" + c.fmap0[key]


  def temprepl(c, s, parse = False, d0 = None):
    ''' template replacement '''
    d = {
      "%OBJ%"   : c.obj,
      "%PFX%"   : c.pfx,
      "%MODE%"  : c.varmode
    }
    if d0:
      d = dict( list(d.items()) + list(d0.items()) )
    for key in d:
      s = s.replace( key, d[key] )
    if parse: # parse the string into lines with '\n'
      s = s.splitlines(True)
    return s


  def mutfunc(c, func, doall = True):
    ''' replace the function name `func' by fmap[func] '''

    # construct a new function name
    if func in c.fmap:
      nfunc = c.fmap[func]
    else:
      nfunc = c.pfx + "_" + func

    # (?<!\w) means that the preceeding character is not a word [a-zA-Z_]
    # (?!\w)  means that the following character is not a word
    prog = re.compile("(?<!\w)" + "(" + func + ")" + "(?!\w)")
    for i in range(len(c.s)):
      # search function name line by line
      m = prog.search(c.s[i])
      if m:
        c.s[i] = c.s[i][:m.start(1)] + nfunc + c.s[i][m.end(1):]
        if not doall:
          break


  def mutfunc2(c, func, iscall = True, newend = None,
               pat = None, verbose = False):
    ''' replace `func' by the new name `fmap[func]'
        also add an object to the last parameter/argument '''

    if not func in c.fmap:
      print "cannot mutate function %s" % func
      raise Exception

    # (?<!\w) means that the preceeding character is not a word [a-zA-Z_]
    # (?!\w)  means that the following character is not a word
    if pat == None:
      pat = "(?<!\w)" + "(" + func + ")" + "\s*\(\s*\w"

    if iscall: # function call, looks like: func(...);
      ending = ");"
      oldend = "\);$"
      if not newend:
        newend = ", %s);" % c.obj
    else: # function definition, looks like: func()\n{
      ending = ")"
      oldend = "\)$"
      if not newend:
        newend = ", %s_t *%s)" % (c.pfx, c.obj)

    if c.findblk(pat, ending = ending, verbose = verbose) < 0:
      print "cannot find [%s], ending [%s]" % (func, ending)
      raise Exception

    # add function prefix
    c.s[c.begin] = re.sub(func, c.fmap[func], c.s[c.begin])

    # append obj to the function argument list
    s_end = c.s[c.end].rstrip() + "\n"
    c.s[c.end] = re.sub(oldend, newend, s_end)



  def rmcopyrite(c):
    c.rmblock("thanx(stderr)")
    c.rmline('#include "copyrite.h"')
    c.rmblock('CopyRight(stderr')
    c.rmline('CopyRight(fplog', doall = True)
    c.rmline('please_cite', doall = True)


  def rmionize(c):
    ''' bIonize (-ionize) is about X-ray bombardment
        we don't need it '''
    c.rmblock('if (bIonize)', doall = True)
    c.rmline('"-ionize",', off1 = 2)
    c.rmline('#include "ionize.h"')
    c.rmline('bool bIonize = FALSE;', doall = True)
    c.rmline('gmx_bool bIonize=FALSE;', doall = True)
    c.rmline('bIonize = (Flags & MD_IONIZE);', doall = True)
    c.rmline('Flags = Flags | (bIonize ? MD_IONIZE : 0);')
    c.substt(',bIonize=FALSE', '')


  def rmtcr(c):
    c.rmblock('if (bTCR', doall = True)
    c.rmline('bTCR = ftp2bSet', wcmt = True)
    c.rmline('t_coupl_rec *tcr=NULL;')
    c.substt('bTCR=FALSE,', '')


  def rmffscan(c):
    c.rmblock('if (bFFscan', doall = True)
    c.substt(' || bFFscan', '')
    c.substt(' && !bFFscan', '', doall = True)
    c.rmline('bFFscan = (Flags & MD_FFSCAN);', doall = True)
    c.rmline('bool bFFscan;')
    c.substt('bFFscan ? step+1 : step', ' step')


  def rmre(c):
    ''' remove replica exchange stuff '''
    c.rmblock('if (repl_ex_nst > 0 && MASTER(cr))')
    c.subs('\(bExchanged \|\| \((.*)\)\)\s*\{', r'(\1) {',
        pat = 'if \(bExchanged \|\| \(ir\-\>nstlist')
    c.rmblock('if (repl_ex_nst != 0 && nmultisim < 2)')
    c.rmblock('if (bExchanged && PAR(cr))')
    c.rmline('bExchanged = FALSE;', wcmt = True, doall = True)
    c.rmblock('print_replica_exchange_statistics(fplog')
    c.rmblock('init_multisystem(cr')
    c.substt(' || bExchanged', '', doall = True)
    c.rmblock('if ((repl_ex_nst > 0) && (step > 0)')
    c.rmline('etINT, {&repl_ex_nst}', off1 = 2)
    c.rmline('etINT, {&repl_ex_seed}', off1 = 2)
    c.rmline('etINT,{&nmultisim}', off1 = 2)
    c.rmline("gmx_repl_ex_t repl_ex=NULL;")
    c.substt(',bExchanged', '', 'bSumEkinhOld,bExchanged;')
    if c.isgmx4:
      c.rmline('int nmultisim=0;')
      c.rmline('int repl_ex_nst=0;')
      c.rmline('int repl_ex_seed=-1;')


  def rmopenmm(c):
    # remove sopenmm
    c.rmline('#include "md_openmm.h"')
    ''' remove code blocks in

          #ifdef GMX_OPENMM
          ...
          #endif
    '''
    c.rmdf(c.sopenmm)


  def rmedv4(c):
    c.rmline('#include "edsam.h"')
    c.rmline('gmx_edsam_t ed;')

    # remove the if block in the runner.c (v4.5) or mdrun.c (v4.0)
    c.rmblock('if (opt2bSet("-ei",NFILE,fnm)) {', ending = 'ed=NULL;')

    # change the argument `ed' to `NULL', such as in calling do_force()
    c.substt(',ed,', ",NULL,", pat = "constr = init_constraints(fplog")


  def rmetc(c):
    c.substt("extern", "extern volatile", pat = "bGotTermSignal, bGotUsr1Signal")
    c.rmblock('"Getting Loaded')
    c.rmblock('"Loaded with')


  def rmglasv4(c):
    ''' bGlas are in md.c and mdrun.c of GROMACS 4
        but they are removed in v4.5 '''
    c.rmblock("if (bGlas)")
    c.rmline('etBOOL,{&bGlas}', off1 = 2)
    c.substt(',bGlas=FALSE', "")
    c.rmline('static bool bGlas = FALSE;')
    c.rmline('bGlas = (Flags & MD_GLAS);')
    c.rmline('Flags = Flags | (bGlas ? MD_GLAS')


  def rmetcv4(c):
    ''' clean up GROMACS v4.0 '''

    c.rmblock('if (ir->efep != efepNO)', doall = True)
    c.rmline('etINT, {&nthreads}', off1 = 2)  # threads don't work in v4.0
    c.substt(',mu_aver=0', '')
    c.substt('boxcopy,lastbox', 'lastbox')
    c.substt(',gnx=0', '', 'a0,a1,gnx=0,ii')
    c.rmline('*xcopy=NULL,*vcopy=NULL')
    c.rmblock('if(fr->bQMMM)', doall = True)
    c.substt('nsteps_done;', 'nsteps_done = 0;', 'int nsteps_done;')
    c.rmblock('if (bSimAnn)')
    c.rmline('atom_id *grpindex=NULL;')
    c.rmline('int nthreads=1;')
    c.substt("global_state", "global_stat")
    c.rmline("GROMACS compiled without threads support", off0 = -2, off1 = 2)


  def rmdbg(c):
    ''' remove debug code '''
    c.rmblock("if (TAKETIME)", doall = True)
    c.rmline("where();", doall = True)
    c.rmline("debug_gmx();", doall = True)
    c.rmblock("if (debug)", doall = True)
    c.rmline("#define TAKETIME FALSE", off0 = -2, off1 = 2)


  def mdcom(c):
    c.rmif0()
    c.rmblock("G R O M A C S", doall = False, starter = "/*", ending = "*/")
    c.shdr()
    c.rmdf("GMX_FAHCORE")
    c.rmcopyrite()
    c.rmionize()
    c.rmtcr()
    c.rmffscan()
    c.rmline("debug_gmx()", doall = True)
    c.rmre()
    c.rmopenmm()
    c.rmetc()

    c.rmline('"pullx", ffOPTWR', doall = True)
    c.rmline('"pullf", ffOPTWR', doall = True)
    if c.isgmx4:
      c.rmedv4()
      c.rmglasv4()
      c.rmetcv4()


  def addifdef(c, pat, tag, extra = ""):
    ''' add #ifdef `tag' ... #endif around a declaration that contains `pat'
        Note that `pat' is a regular expression '''
    ii = c.findln(pat)
    if ii >= 0:
      c.s = c.s[:ii] + ["#ifdef %s\n" % tag,
        c.s[ii] + extra, "#endif\n", ] + c.s[ii+1 : ]


  def getversion(c):
    ''' same as the static function getgmxver()
        but it avoids repeated calling that '''

    # this function has been called before
    # no need to repeat it for every instance
    if CCGMX.version > 0: return

    # we assign the class-level, instead of instance-level,
    # variables, such that, when the next time an instance
    # is create, the version needs not to be computed again
    (CCGMX.version, CCGMX.rootdir, CCGMX.isgmx4, CCGMX.sopenmm
        ) = getgmxver()

    # verify the instance has the values
    print "GROMACS version: %d; v4.0? %s; OPENMM string: %s" % (
        c.version, bool(c.isgmx4), c.sopenmm)



class CCAT(CCGMX):
  ''' here we extend CCGMX (ccgmx.py), which, in turn, is an extension
      of CC (cc.py) '''

  ''' the constructor is that of CCGMX
      c = CCAT(filename, hdrs) '''

  def __init__(c, fn, atver, hdrs = {}):
    CCGMX.__init__(c, fn, str(atver), hdrs)
    c.fmap["do_md"] = "md"
    c.fmap["mdrunner"] = "runner"

  def mdcommon(c):
    c.mdcom()
    if use_at:
      c.addpremode()

  def rmjunkvarsv4(c, begin, end):
    ''' remove variable declarations for replica exchange and edsam '''
    if c.isgmx4:
      for i in range(c.begin, c.end+1):
        c.s[i] = re.sub(r"gmx_edsam_t\s*ed,", "", c.s[i])
        c.s[i] = re.sub(r"int\s*repl_ex_nst,", "", c.s[i])
        c.s[i] = re.sub(r"int\s*repl_ex_seed,", "", c.s[i])


  def rmjunkargsv4(c, begin, end):
    ''' remove arguments for replica exchange and edsam '''
    if c.isgmx4:
      for i in range(c.begin, c.end+1):
        c.s[i] = re.sub(r"(?<!\w)ed\,", "", c.s[i])
        c.s[i] = re.sub(r"(?<!\w)repl_ex_nst\,", "", c.s[i])
        c.s[i] = re.sub(r"(?<!\w)repl_ex_seed\,", "", c.s[i])


  def addpremode(c):
    if c.findblk("struct\s*mdrunner_arglist") >= 0:
      c.addln(c.end-1, "    int atpremode;\n")
    if c.findblk("t_commrec\s*\*mdrunner_start_threads", ending = ")") >= 0:
      c.s[c.end] = c.s[c.end].rstrip()[:-1] + ", int atpremode)\n"
    if c.findblk('\w\s*\=\s*mdrunner_start_threads\(', ending = ");") >= 0:
      c.s[c.end] = re.sub(r"Flags\)", "Flags, atpremode)", c.s[c.end])
    if c.findblk('mda\-\>Flags=Flags;', ending = None) >= 0:
      c.addln(c.end+1, '    mda->atpremode = atpremode;\n')

    # in v4.5 runner() is called in runner.c and mdrun.c
    if c.isgmx4: pat = '^\s*runner\('
    else: pat = '^\s*[\w\-\>]+\s*\=\s*runner\('
    if c.findblk(pat, ending = ";") >= 0:
      c.rmjunkargsv4(c.begin, c.end)
      if c.s[c.end].find("mc.") >= 0: pfx = "mc."
      else: pfx = ""
      c.s[c.end] = c.s[c.end].rstrip()[:-2] + ", %satpremode);\n" % pfx

  def addatopt(c):
    if c.findblk('t_filenm\s*fnm\[\]', ending = "};") >= 0:
      if not c.s[c.end - 1].endswith(","):
        c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
      c.addln(c.end, '    { efMDP, "-at",     NULL,       ffOPTRD } /* at.cfg */\n')

    if c.findblk('t_pargs\s*pa\[\]', ending = "};") >= 0:
      if not c.s[c.end - 1].endswith(","):
        c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
      c.addln(c.end, [ '    { "-atpre", FALSE, etINT, {&atpremode},\n',
            '      "a preparation run" }\n' ])
      if c.s[c.begin].strip().startswith("static"): pfx = "static "
      else: pfx = ""
      c.addln(c.begin, '  %sint atpremode = 0;\n' % pfx)


  def atinclude(c):
    i = len(c.s) - 1
    while i >= 0:
      if c.s[i].startswith("#include"):
        while i < len(c.s):
          if not c.s[i].startswith("#"): break
          i += 1
        break
      i -= 1
    c.addln(i, ['\n', '#define GMXVERSION %s\n' % c.version,
        '#include "md%dutil.h"\n' % atver])


  def atrunner(c):
    ''' modify the function runner() '''

    # search the declaration of the function runner()
    if c.findblk("^\s*\w+\s+" + "mdrunner" + "\(", ending = ")") < 0:
      print "cannot find runner() in %s" % c.filename
      open("dump.txt", "w").writelines(c.s)
      raise Exception
    runnerbegin = c.begin

    c.rmjunkvarsv4(c.begin, c.end)
    if use_at:
      c.s[c.begin] = re.sub("mdrunner", "runner", c.s[c.begin])
      c.s[c.end] = re.sub("\)$", ", int atpremode)", c.s[c.end])
      # add at declaration
      for j in range(c.end, len(c.s)):
        if c.s[j].strip() == "":
          break
      else: raise Exception
      c.addln(j, "    at_t *at;\n")
    proto = c.s[c.begin : c.end + 1]
    last = len(proto) - 1
    proto[last] = proto[last].rstrip() + ";\n"
    proto = ["\n", "/* declare runner() before mdrunner_start_fn() uses it */\n",
        "static\n"] + proto

    # add at_close() at the end
    k1, k2 = c.getblockend(c.end, sindent = "", ending = "}", wsp = False)
    if not c.isgmx4: k2 -= 1
    c.addln(k2, ["    if(at != NULL)\n",
        "      at_close(at); /* release memory */\n"])

    # add atgmx_init()
    for i in range(k1, k2):
      if c.s[i].strip().startswith("snew(nrnb,"):
        break
    else: raise Exception
    atinit = ['    ' + s + '\n' for s in
r'''/* initialize aT, cr->duty PP/PME has been assigned */
at = atgmx_init(atgmx_opt2fn("-at", nfile, fnm),
    Flags & MD_STARTFROMCPT,
    mtop, inputrec, cr, atpremode);
if ((cr->duty & DUTY_PP) && at == NULL) {
  fprintf(stderr, "Failed to initialize aT\n");
  exit(1);
}'''.splitlines()] + ['\n']
    c.addln(i, atinit)

    # VI. change the call do_md()
    c.substt("integrator[inputrec->eI].func", "do_md", doall = True)
    if c.findblk("(?<!\w)do_md\(", ending = ");", i0 = k1) >= 0:
      c.s[c.begin] = re.sub("do_md", c.fmap["do_md"], c.s[c.begin])
      c.rmjunkargsv4(c.begin, c.end)
      c.s[c.end] = re.sub("\);$", ", at);", c.s[c.end])

    if not c.isgmx4:
      c.mutfunc2("mdrunner", iscall = True, newend = ", mc.atpremode);")
    c.mutfunc("mdrunner")

    # insert prototype after a bunch of #include
    if not c.isgmx4:
      i = len(c.s) - 1
      while i >= 0:
        if c.s[i].startswith("#include"):
          break
        i -= 1
      c.addln(i+1, proto)
    else:
      # for v4.0, both mdrunner() and do_md() are in md.c
      # move mdrunner() to the end of file, so it can call do_md()
      k1, k2 = c.getblockend(runnerbegin, sindent = "", ending = "}", wsp = True)
      c.s = c.s[:k1] + c.s[k2+1:] + c.s[k1:k2+1]

    # these calls will mess up with the index,
    # so we add them here
    c.rmline('const gmx_intp_t', wcmt = True, doall = True)
    c.rmline('gmx_integrator_t *func;', off0 = -1, off1 = 3)

  def atmd(c):
    if c.findblk("\w+\s+" + "do_md" + "\(", ending = ")") < 0:
      open("dump.txt", "w").writelines(c.s)
      raise Exception
    c.s[c.begin] = re.sub("do_md", c.fmap["do_md"], c.s[c.begin])
    c.s[c.end] = re.sub("\)$", ", at_t *at)", c.s[c.end])
    if c.isgmx4:
      c.rmjunkvarsv4(c.begin, c.end)

    # insert atmove()
    if c.isgmx4:
      sig = "Time for performance"
      offset = -1
    else:
      sig = "END PREPARING EDR OUTPUT"
      offset = 1
    if c.findblk("/\*.*" + sig, ending = None) < 0: raise Exception
    i = c.begin + offset
    if c.isgmx4: bxtc = "bXTC"
    else: bxtc = "mdof_flags & MDOF_XTC"
    atmove = r'''
/* update data, temperature and change at->scale */
if (0 != atgmx_move(at, enerd,
     step, bFirstStep, bLastStep, bGStat,
     %s, bNS, cr) ) {
  exit(1);
}''' % bxtc
    atmove = [' '*8 + s + '\n' for s in atmove.splitlines()] + ['\n']
    c.addln(i, atmove)

    if c.findblk("(?<!\w)do_force\(", ending = ");") < 0: raise Exception
    c.s[c.begin] = re.sub("do_force", "atgmx_doforce", c.s[c.begin])
    c.s[c.end] = re.sub("\);$", ", at);", c.s[c.end])
    c.subs(",ed,", ",NULL,", "fp_field,ed,")


def do_md(fn):
  ''' read md.c, and remove junks '''
  c = CCAT(fn, atver, hdrs = {})
  c.mdcommon()
  if c.isgmx4:
    c.atrunner()
    c.substt("if (do_md == do_md) {", "{")
  if use_at:
    c.atinclude()
    c.atmd()
  return c


def do_runner(fn, hdrs):
  ''' read runner.c, and remove junks '''
  c = CCAT(fn, atver, hdrs)
  c.mdcommon()
  if not c.isgmx4:
    c.atrunner()
    # clear the ``md == md'' mess
    c.rmline("if (do_md == do_md", off1 = 2)
  return c


def do_mdrun(src, hdrs):
  ''' read mdrun.c, and remove junks '''
  c = CCAT(src, atver, hdrs)
  c.mdcommon()
  if c.findblk("static const char \*desc", ending = None) >= 0:
    pfx = "static "
  else:
    pfx = ""
  c.rmblock("char *desc[]", starter = None, ending = "};", verbose = False)

  c.addln(c.begin, '  %sconst char *desc[] = {"self-contained GROMACS"};\n' % pfx)
  if use_at: c.addatopt()

  c.mutfunc2("mdrunner", iscall = True, newend = ", atpremode);")
  c.mutfunc("mdrunner")

  if c.isgmx4:
    c.substt(",ed,repl_ex_nst,repl_ex_seed,pforce", ",pforce")
  return c


def mkmdx_c(atver):
  ''' combine source code from md.c, runner.c and mdrun.c '''

  md_c = do_md("md.c")
  if not md_c.isgmx4: runner_c = "runner.c"
  else: runner_c = None
  runner_c = do_runner( runner_c, md_c.hdrs )
  mdrun_c = do_mdrun("mdrun.c", runner_c.hdrs )

  # combine md.c, runner.c and mdrun.c
  s = md_c.s + runner_c.s + mdrun_c.s

  if md_c.isgmx4: fnout = "md%da.c" % atver
  else: fnout = "md%d.c" % atver

  savelines(fnout, s)


class CCutil(CC):
  ''' functions about creating mdxutil.c '''

  def funcdef(c, funcpat, wsp = True):
    ''' return the entire function body (definition) '''
    if c.findblk(funcpat, ending = "}", wsp = wsp) < 0:
      print "cannot find function pattern %s" % funcpat
      raise Exception
    return c.s[c.begin : c.end+1]


  def common(c):
    c.rmblock("if (TAKETIME)", doall = True)
    c.rmline("where();", doall = True)
    c.rmline("debug_gmx();", doall = True)
    c.rmblock("if (debug)", doall = True)
    c.rmline("#define TAKETIME FALSE", off0 = -2, off1 = 2)


  def rmdoforce(c):
    ''' remove existing atgmx_doforce() definition '''
    c.rmblk("^\w+\s*atgmx_doforce\(", starter = None, ending = "}", wsp = False)


  def rmforcelow(c):
    ''' remove existing atgmx_forcelow() definition '''
    c.rmblk("^\w+\s*atgmx_forcelow\(", starter = None, ending = "}", wsp = False)


  def doforcedef(c, forcescl = None):
    ''' modify do_force() definition, c = force.c '''
    c.common()
    nm = "do_force"
    if c.findblk("^\w+\s*" + nm + "\(", ending = ")", wsp = False) < 0:
      print "cannot find %s in %s" % (nm, c.filename)
      raise Exception
    fbegin = c.begin
    proto = c.s[fbegin : c.end + 1]
    last = len(proto) - 1
    proto0 = proto[0: last+1] # make a copy
    proto[0] = re.sub(nm, "atgmx_doforce", proto[0])
    proto[last] = re.sub("\)\s*$", ", at_t *at)\n", proto[last]) # add at
    if atver == 2:
      # 1. add force scaling code
      if c.findblk("\/\* Sum the potential energy terms from group contributions \*\/",
          ending = None, wsp = False) < 0: raise Exception
      if not forcescl: raise Exception
      c.addln(c.begin, forcescl)

      # 2. do_force_lowlevel --> atgmx_forcelow
      if c.findblk("^\s*do_force_lowlevel\(", ending = ");", wsp = False) < 0:
        raise Exception
      c.s[c.begin] = re.sub("do_force_lowlevel", "atgmx_forcelow", c.s[c.begin])
      c.s[c.end] = re.sub("\);$", ", at);", c.s[c.end])
      #print ''.join(c.s[c.begin:c.end+1]),; raw_input()
      c.addln(c.begin - 1, "    at_clear_energy(at);\n")

    b0, b1 = c.getblockend(fbegin, sindent = "", ending = "}", wsp = False)
    body = proto + c.s[b0 + last + 1: b1 + 1]
    return proto0, proto, body


  def forcelowdef(c):
    ''' existing forcelow() definition '''
    c.common()
    nm = "do_force_lowlevel"
    if c.findblk("^\w+\s*%s\(" % nm, ending = ")", wsp = False) < 0:
      print "cannot find %s in %s" % (nm, c.filename)
      raise Exception
    fbegin = c.begin
    proto = c.s[fbegin : c.end + 1]
    last = len(proto) - 1
    proto0 = proto[0: last+1] # make a copy
    proto[0] = re.sub(nm, "atgmx_forcelow", proto[0])
    proto[last] = re.sub("\)\s*$", ", at_t *at)\n", proto[last]) # add at

    if atver == 2:
      if c.findblk("^\s*calc_bonds\(", ending = ");", wsp = False) < 0:
        raise Exception
      c.s[c.begin] = re.sub("calc_bonds\(", "atgmx_calcbonds(cr, ", c.s[c.begin])
      c.s[c.end] = re.sub("\);$", ", at);", c.s[c.end])
      for l in range(c.begin, c.end+1):
        c.s[l] = re.sub("\s*atype,\s*born,", "", c.s[l])

    b0, b1 = c.getblockend(fbegin, sindent = "", ending = "}", wsp = False)
    body = proto + c.s[b0 + last + 1: b1 + 1]
    return proto0, proto, body

  def get_forcescl(c):
    ''' get force scaling code '''
    if c.findblk("\w+\s+atgmx_doforce\(", ending = "}") < 0:
      print "cannot find atgmx_doforce"
      raise Exception
    df0 = c.begin
    df1 = c.end
    if c.findblk("if\s*\(at\s*\!=\s*NULL\)\s*\{", wsp = True, i0 = df0, i1 = df1) < 0:
      print "cannot find force scaling code"
      raise Exception
    return c.s[c.begin : c.end]



def do_force1(cu, cf):
  '''
  cu: md1util.h
  cf: sim_util.h, which contains force()
  '''

  # create atgmx_doforce that calls force()
  call, proto, body = cf.doforcedef()
  #print ''.join(call), ''.join(atproto); raw_input()
  # call do_force()
  call = proto2call(call)
  decl = ["  int k;\n", "  real scl;\n", "\n"]
  forcescl = [s + '\n' for s in
  '''
  /* scale the force */
  scl = (real) at->scale;
  for (k = mdatoms->start; k < mdatoms->start + mdatoms->homenr; k++) {
    f[k][0] *= scl;
    f[k][1] *= scl;
    f[k][2] *= scl;
  }'''.splitlines()]
  f = proto + ["{\n"] + decl + call + forcescl + ["}\n"]

  cu.rmdoforce()

  # add the force to the declaration
  cu.addln(cu.begin, f)



def do_force2(cu, cf, cl):
  ''' cu: md2util.h,  cf: sim_util.c, cl: force.c '''
  call, proto, f = cl.forcelowdef()
  cu.rmforcelow()
  cu.addln(cu.begin, f)

  forcescl = cu.get_forcescl()
  #print "force scaling code:\n", ''.join(forcescl); raw_input()
  call, proto, f = cf.doforcedef(forcescl)

  # collect auxillary static functions from sim_util.c
  aux = []
  aux += cf.funcdef("^static void sum_forces\(")
  aux += cf.funcdef("^static void calc_f_el\(")
  aux += cf.funcdef("^static void calc_virial\(")
  aux += cf.funcdef("^static void print_large_forces\(")

  # remove atgmx_doforce() from md2util.h
  # add the new force
  cu.rmdoforce()
  cu.addln(cu.begin, aux +  f)


def do_force(cu, cf, cl):
  ''' cu: mdxutil.c
      cf: sim_util.c which contains force()
      cl: forcelow.c
  '''
  if atver == 1:
    do_force1(cu, cf)
  elif atver == 2:
    do_force2(cu, cf, cl)


def mkmdxutil_h(atver):
  pfx = "md%s" % atver
  fnin = pfx + "uti.h"
  fnout = pfx + "util.h"
  if not os.path.exists(fnin):  # v4.0
    fnin = pfx + "util.h"
    fnout = pfx + "utila.h"

  util_h = CCutil(fnin)
  force_c = CCutil("../mdlib/sim_util.c")
  if atver == 2:
    forcelow_c = CCutil("../mdlib/force.c")
  else: forcelow_c = None

  do_force(util_h, force_c, forcelow_c)

  savelines(fnout, util_h.s)


def usage():
  ''' print usage and die '''
  print ''' The program grabs several core GROMACS files,
  kernel/runner.c, kernel/mdrun.c, kernel/md.c, mdlib/sim_util.c, mdlib/force.c
  and manages to construct a self-contained MD engine.
  Code from these files are copied to
  mdx.c and mdxutil.h, where X = 1 or 2
  md1.c and md2.c are mostly the same,
  md2util.h embeds the content of do_force(), while md1util.h only make the function call
  for GROMACS 4.5, mdxutil.h are constructed from modifying the corresponding GROMACS 4.0 verion
  which are linked to here in different names as mdxuti.h
  '''
  exit(0)


def doargs():
  ''' Handle args '''
  global tgtver

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hv:",
         ["help", "version=", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  for o, a in opts:
    if o in ("-v", "--version"):
      tgtver = range(1, int(a)+1)
    elif o in ("-h", "--help"):
      usage()


def main():
  ''' main function '''
  global atver
  doargs()
  for atver in tgtver:
    mkmdx_c(atver)
    mkmdxutil_h(atver)


if __name__ == "__main__":
  # try to accelerate the process
  try:
    import psyco
    psyco.full()
  except ImportError: pass

  main()


#!/usr/bin/env python

'''
generic code changer
'''

import re, getopt, sys

class CC:
  ''' Code-Changer: basic class of changing code '''
  def __init__(c, fn, hdrs = {}):
    ''' read the input file `fn', parses into lines `c.s'
        hdrs gives a list of existing #include
        `fn' can alternatively be a list of strings '''
    if type(fn) == list: # input is a list of strings
      c.s = fn
    elif fn:  # input is a file name
      c.s = open(fn).readlines()
    else: c.s = []
    c.s = CC.tab2sp0( c.s ) # tabs to spaces
    c.fname = fn
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
    if find: find = CC.tore0(find)
    if pat: pat = CC.tore0(pat)
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
    return c.findln(CC.tore0(pat))


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
    c.rmln(CC.tore0(pat), doall, off0, off1, wcmt, verbose)


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
            start, c.fname, c.s[start])
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
          ending, start, c.fname)
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
        if verbose: print "pattern %s found in %s" % (pat, c.fname)

        # search a starting if
        start = i
        if starter != None:
          while start >= 0:
            if c.s[start].lstrip().startswith(starter):
              break
            start -= 1
          else:
            print "no starter [%s], pat [%s], line %d, file %s" % (
                starter, pat, i, c.fname)
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
    c.rmblk(CC.tore0(pat), doall, starter, ending, sindent, wsp, verbose)


  def rmfunc(c, functag, starter = None, regex = False):
    ''' convenient wrapper for rmblk '''
    if not regex: functag = CC.tore0(functag)
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
    return c.findblk(CC.tore0(pat), ending, wsp, i0, i1, sindent, verbose)


  def shdr(c):
    ''' remove redundant #include '''
    i = 0
    while i < len(c.s):
      if c.s[i].startswith("#include"):
        key = c.s[i].strip()
        if key in c.hdrs:
          start = i
          end = i+1
          # a preceeding comment
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
          print "preprocessor in file %s (%d : %d)" % (c.fname, i, j+1)
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


  """  deprecated
  def addfuncpfx(c, func, pfx, newfunc = None, doall = True):
    ''' change the function name `func' to `pfx_func'
    or if `newfunc' is given, change `func' to `pfx_newfunc' '''

    if newfunc == None: newfunc = func

    for i in range(len(c.s)):
      # search function name line by line
      pos0 = c.s[i].find(func)
      if pos0 >= 0:
        pos1 = pos0 + len(func)
        c.s[i] = c.s[i][:pos0] + pfx + "_" + newfunc + c.s[i][pos1:]
        if not doall:
          break
  """


  @staticmethod
  def beautify0(s):
    ''' make code better looking by calling cspacer.py '''
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


  @staticmethod
  def tore0(s):
    ''' many functions below require the pattern to be regular expression
    so this function translate a plain string to a regular expression '''

    # list of special characters
    special = "+-*.&?!|()[]{}"
    t = ""
    s = re.sub("\s+", " ", s)
    for i in range(len(s)):
      if s[i].isspace(): c = "\s+"
      elif s[i] in special: c = "\\" + s[i]
      else: c = s[i]
      t += c
    return t


  @staticmethod
  def save0(fn, s):
    # remove trailing spaces
    s = CC.beautify0(s)
    print "writing %s ..." % fn
    open(fn, "w").writelines(s)


  @staticmethod
  def temprepl0(s, d):
    ''' for each key in `d' in `s', replace it with `d[key]' '''
    snew = s
    for key in d:
      snew = snew.replace(key, d[key])
    return snew


  @staticmethod
  def tagrepl0(templ, d):
    '''
    in the template `templ' (string array)
    find a line that matches any key in `d'
    replace it with `d[key]' (string array)
    This is used in the main function of spider0v45.py
    '''
    src = templ[:]
    for key in d:
      spat = key.strip()
      # strip away `/*' and `*/'
      if spat.startswith("/*"): spat = spat[2:].lstrip()
      if spat.endswith("*/"): spat = spat[:-2].rstrip()

      # `srep' is the actual code, as a string
      srep = ''.join(d[key]) + '\n'

      for i in range(len(src)):
        ln = src[i].strip()
        if (ln.find(spat) >= 0 and ln.startswith("/*")):
          src[i] = srep   # replace it with `srep'

      '''
      # the regular expression version is slow and troublesome
      spat = CC.tore0(key)
      prog = re.compile(spat)
      srep = ''.join(d[key]) + '\n'
      for i in range(len(src)):
        if prog.search(src[i]):
          src[i] = prog.sub(srep, src[i])
      '''
      src = ''.join(src).splitlines(True)
    return src


  @staticmethod
  def gettabsize0(s, deftab = 4):
    ''' try to guess the tab size '''
    for i in range(len(s)):
      if len(s[i].strip()): break
    else:
      return deftab

    ln = s[i]
    if not ln.startswith("/*"): return deftab
    m = re.search("tab-width:\s*([0-9]+);", ln)
    if not m: return deftab
    return int(m.group(1))


  @staticmethod
  def tab2sp0(s, ntab = 0):
    ''' convert tabs to spaces '''
    if ntab <= 0:
      ntab = CC.gettabsize0(s)
      #print "tab size is %d" % ntab

    for i in range(len(s)):
      while True: # don't stop until every tab is killed
        pos = s[i].find("\t")
        if pos < 0: break
        nsp = ntab - (pos % ntab)
        s[i] = s[i][:pos] + " " * nsp + s[i][pos+1:]
    return s

  @staticmethod
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


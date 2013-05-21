#!/usr/bin/env python

''' a generic code changer class '''

import re, getopt, os, sys
from ccutil import str2re, tab2sp

class CC:
  ''' Code-Changer: basic class of changing code '''

  def __init__(c, fn, src = [], hdrs = {}):
    ''' `fn' is the input file name
        `src' is the input text, if not given explicitly
              `src' is read from the input file `fn'
        `hdrs' gives a list of existing #include '''
    c.s = src
    if not src and fn: # src is not explicitly specified
      c.s = open(fn).readlines()

    if type(c.s) == str: # input is a string
      c.s = c.s.splitlines(True) # parse it to lines

    c.s = tab2sp( c.s ) # tabs to spaces
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



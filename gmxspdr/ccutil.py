#!/usr/bin/env python

''' standalone and reusable functions

    export function list:
    * str2re(s)               translate a plain string to R.E.
    * tab2sp(s)               convert leading tabs to spaces
    * savelines(fn, s)        save `s' to file `fn'
    * tmphastag(templ, tag)   if `templ' has `tag' or its variants
    * tmptagrep(templ, dic)   replace tags by the dictionary values
    * getgmx()                get GROMACS version, rootdir, ...
'''

import re, os, sys


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



def gettabsize(s, deftab = 4):
  ''' try to guess the tab size '''
  if not s: return deftab
  for i in range(len(s)):
    if len(s[i].strip()): break
  else:
    return deftab

  ln = s[i]
  if not ln.startswith("/*"): return deftab
  m = re.search("tab-width:\s*([0-9]+);", ln)
  if not m: return deftab
  return int(m.group(1))


def tab2sp(s, ntab = 0):
  ''' convert tabs to spaces '''

  if not s: return s
  if ntab <= 0:
    ntab = gettabsize(s)
    #print "tab size is %d" % ntab

  for i in range(len(s)):
    while True: # don't stop until every tab is killed
      pos = s[i].find("\t")
      if pos < 0: break
      nsp = ntab - (pos % ntab)
      s[i] = s[i][:pos] + " " * nsp + s[i][pos+1:]
  return s



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





def mktagvars(tag):
  ''' make a list of variants of `tag '''

  if tag[0] != '%': tag = '%' + tag + '%'
  ls = [tag, tag.lower(), tag.upper()]
  ls += [ d.replace('.', '_') for d in ls ]
  ls += [ d.replace('_', '.') for d in ls ]
  return list(set(ls)) # remove duplicates


def tmphastag(templ, tag):
  ''' detect if the template `templ' has `tag' or its variants '''

  tags = mktagvars(tag)
  src = templ
  if type(src) == list: src = ''.join(templ)
  for k in tags:
    if src.find(k) >= 0: return 1
  return 0



def addcasekeys(d):
  ''' add lower-case keys to the dictionary '''
  newd = {}
  for key in d:
    newd[key] = d[key]
    newd[key.lower()] = d[key]
    newd[key.upper()] = d[key]
  return newd


def addvarkeys(d):
  ''' add keys with `_' replaced by `.' to the dictionary '''
  newd = {}
  for key in d:
    newd[key] = d[key]
    newd[key.replace('_', '.')] = d[key]
    newd[key.replace('.', '_')] = d[key]
  return newd


def tmptagrep(templ, d0):
  ''' in the template `templ' (string array)
      replace any line that matches any key in `d'
      by `d[key]' (string array) '''

  # add variance of the keys
  d = addvarkeys( addcasekeys(d0) )
  #raw_input( str(list(d)) )

  # make a copy of the template
  if type(templ) == str: src = templ.splitlines(True)
  else: src = templ[:]

  for i in range(len(src)):
    ln = src[i].strip()
    # find keywords wrapped in comments `/* */'
    if ln.startswith("/*") and ln.endswith("*/"):
      for key in d:
        if ln.find(key) >= 0:
          src[i] = ''.join(d[key]) + '\n'
  return ''.join(src).splitlines(True)



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



def getgmx():
  ''' get GROMACS version from CMakeLists.txt or configure.ac
      should work for GROMACS v4.0+ '''

  # only do the computation for the first time
  if not hasattr(getgmx, "root"):
    getgmx.root = getgmxroot()
    cfgac = os.path.join(getgmx.root, "configure.ac")
    cmake = os.path.join(getgmx.root, "CMakeLists.txt")

    if os.path.exists(cmake): # version 4.5 or later
      for s in open(cmake):
        # we look for a line that looks like
        # set(PROJECT_VERSION "4.5.6-dev")
        # here, `*?' in `(.*?)' means a non-greedy search
        #       `(-dev)?' means `-dev' may or may not present
        m = re.search(r'PROJECT_VERSION\s+"(.*?)(-dev)?"', s.strip())
        if m:
          # convert 4.5.6 to an integer 40506
          sver = m.group(1).replace('.', '0')
          if len(sver) < 5: sver += '0' * (5 - len(sver))
          getgmx.ver = int(sver)
          break
      else:
        print "cannot determine version from %s" % cmake
        raise Exception
      getgmx.isv4 = 0
      getgmx.sopenmm = "GMX_OPENMM"

    elif os.path.exists(cfgac): # v4.0 or earlier
      for s in open(cfgac):
        if s.startswith("AC_INIT"):
          # a typical string is like
          #   AC_INIT(gromacs, 4.5.6-dev, [gmx-users@gromacs.org])
          # so we take the second number in the parentheses ()
          sver = s.split(",")[1].strip()
          if sver.startswith("4.0"):
            getgmx.isv4 = 1
          else:
            print "bad gromacs version %s" % sver
            raise Exception
          getgmx.ver = int( sver[:5].replace(".", "0") )
          break
      else:
        print "cannot determine version number from %s" % cfgac
        raise Exception
      # v4 uses USE_OPENMM
      getgmx.sopenmm = "USE_OPENMM"

    else:  # there is neither CMakeLists.txt nor configure.ac
      raise Exception

    print "GROMACS version: %d; root %s, v4.0? %s; OPENMM string: %s" % (
        getgmx.ver, getgmx.root, bool(getgmx.isv4), getgmx.sopenmm)

  return getgmx.ver, getgmx.root, getgmx.isv4, getgmx.sopenmm



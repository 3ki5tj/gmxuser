#!/usr/bin/env python

''' standalone and reusable functions

    export function list:
    * str2re(s)               translate a plain string to R.E.
    * tab2sp(s)               convert leading tabs to spaces
    * savelines(fn, s)        save `s' to file `fn'
    * tmphastag(templ, tag)   if `templ' has `tag' or its variants
    * tmptagrep(templ, dic)   replace tags by the dictionary values
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


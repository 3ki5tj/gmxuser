#!/usr/bin/env python

''' standalone and reusable functions

    export function list:
    * str2re(s)               translate a plain string to R.E.
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



def mkcodepretty(s, verbose = 0):
  ''' make code better looking by calling cspacer.py
      `s' is a string array (lines) '''

  # remove trailing spaces first
  s = [ln.rstrip() + '\n' for ln in s]

  try: # reindent the file
    import cindent
    s = cindent.reindent(s, verbose = verbose)
  except ImportError:
    print "cannot find cindent.py"
    pass

  try: # improve code by cspacer
    import cspacer
    # turn on advanced options
    cspacer.use_rule_add = True
    cspacer.use_rule_paren2 = True
    cspacer.use_rule_nocppcmt = True
    s, nchanges = cspacer.addspace(s, verbose = verbose)
  except ImportError:
    print "cannot find cspacer.py"
    pass

  return s



def savelines(fn, s, verbose = 0):
  s = mkcodepretty(s, verbose)
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
    if k in src: return 1
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
        if key in ln:
          src[i] = ''.join(d[key])
  return ''.join(src).splitlines(True)




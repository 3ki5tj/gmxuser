#!/usr/bin/env python

'''
common GROMACS related modifications,
such as retrieving version and removing useless code
arranged in a class CCGMX, which extends CC
'''

import re, getopt, os, sys
import gmxcom
from cc import CC


class CCGMX(CC):
  ''' changes for GROMACS, it extends class CC (cc.py)
  each source code file (e.g., md.c, mdrun.c)
  should declare an instance via the constructor '''

  # use a dictionary to map function names
  # the mapped functions will further have a prefix
  fmap0 = {
    "mdrunner":           "runner",
    "do_md":              "domd",
    "do_force":           "doforce",
    "do_force_lowlevel":  "doforcelow",
    "do_force_cutsVERLET":"doforcecv",
    "do_force_cutsGROUP": "doforcecg",
    "calc_bonds":         "calcbonds",
  }

  obj = "foo" # common variable name for the object
  pfx = "gmxfoo" # function prefix
  varmode = "gmxfoo_mode" # this is an integer parameter
                          # passed from the command line
                          # better mechanism in the future?

  def __init__(c, fn, src = [],
               obj = None, pfx = None, hdrs = {}):
    # call the base class constructor
    CC.__init__(c, fn, src, hdrs)

    if obj != None:
      c.obj = obj
      c.pfx = pfx
      if not c.pfx:
        c.pfx = "gmx" + c.obj
      c.varmode = c.pfx + "_mode"
      c.fmap = {}
      for key in c.fmap0: # add a prefix
        c.fmap[key] = c.pfx + "_" + c.fmap0[key]


  def temprepl(c, s, parse = False, d0 = None):
    ''' template replacement '''
    d = {
      "%obj%"   : c.obj,
      "%pfx%"   : c.pfx,
      "%mode%"  : c.varmode
    }
    if d0:
      d = dict( list(d.items()) + list(d0.items()) )

    for key in d:
      s = s.replace( key, d[key] )
      s = s.replace( key.upper(), d[key] )
      s = s.replace( key.lower(), d[key] )

    if parse: # parse the string into lines with '\n'
      s = s.splitlines(True)
    return s


  def mutfunc(c, func, doall = True, pfx = None):
    ''' replace the function name `func' by fmap[func] '''

    if pfx == None: pfx = c.pfx

    # construct a new function name
    if func in c.fmap:
      nfunc = c.fmap[func]
    else:
      nfunc = pfx + "_" + func

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
      pat = "(?<!\w)" + "(" + func + ")" + "\s*\(\s*[\w&*]"

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
      print "cc.mutfunc2: no [%s|%s], ending [%s], %s" % (
          func, pat, ending, c.fn)
      return
      #raise Exception

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
    c.subs('bTCR\s*=\s*FALSE,', '')


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
    c.rmblock('if (repl_ex_nst != 0 && nmultisim < 2)')
    c.rmblock('print_replica_exchange_statistics(fplog')
    c.rmblock('init_multisystem(cr')
    c.rmblock('if ((repl_ex_nst > 0) && (step > 0)')
    c.rmline('etINT, {&repl_ex_nst}', off1 = 2)
    c.rmline('etINT, {&repl_ex_seed}', off1 = 2)
    c.rmline('etINT,{&nmultisim}', off1 = 2)
    c.rmline("gmx_repl_ex_t repl_ex=NULL;")

    if gmxcom.isver4():
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
    c.rmdf(gmxcom.sopenmm())


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
    if gmxcom.isver4():
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


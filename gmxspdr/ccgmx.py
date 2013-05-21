#!/usr/bin/env python

'''
GROMACS related code changer, a class CCGMX
also has a function getgmxver
'''

import re, getopt, os, sys
from cc import CC

class CCGMX(CC):
  ''' changes for GROMACS, it extends class CC (cc.py)
  each source code file (e.g., md.c, mdrun.c)
  should declare an instance via the constructor '''

  # these are class level ``static'' variables
  # so they can be set only once
  isgmx4 = 0
  version = -1
  sopenmm = "GMX_OPENMM"

  # use a dictionary to map function names
  fmap0 = {
    "mdrunner": "runner",
    "do_md": "domd",
    "do_force": "doforce",
    "do_force_lowlevel": "doforcelow",
    "calc_bonds": "calcbonds",
  }

  obj = "foo" # project name
  pfx = "foogmx" # function prefix
  varmode = "foomode" # this is an integer parameter to be passed from the command line

  def __init__(c, fn, obj = None, hdrs = {}):
    # call the base class constructor
    CC.__init__(c, fn, hdrs)
    # compute GROMACS version
    c.getgmxver()
    if obj != None:
      c.obj = obj
      c.pfx = c.obj + "gmx"
      c.varmode = c.obj + "mode"
      c.fmap = {}
      for key in c.fmap0: # add a prefix
        c.fmap[key] = c.pfx + "_" + c.fmap0[key]


  def temprepl(c, s, parse = False):
    ''' template replacement '''
    d = {
      "$OBJ" : c.obj,
      "$PFX" : c.pfx,
      "$MODE" : c.varmode
    }
    s = CC.temprepl0(s, d)
    if parse:
      s = [ln + '\n' for ln in s.splitlines()] + ['\n']
    return s

  def mutfunc(c, func, doall = True):
    ''' replace the function name `func' by fmap[func] '''
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


  def mutfunc2(c, func, iscall = True, newend = None, pat = None, verbose = False):
    ''' replace `func' by fmap[func] '''
    if not func in c.fmap:
      print "cannot mutate function %s" % func

    # (?<!\w) means that the preceeding character is not a word [a-zA-Z_]
    # (?!\w)  means that the following character is not a word
    if pat == None:
      pat = "(?<!\w)" + "(" + func + ")" + "\s*\(\s*\w"

    if iscall:
      ending = ");"
      oldend = "\);$"
      if not newend:
        newend = ", %s);" % c.obj
    else:
      ending = ")"
      oldend = "\)$"
      if not newend:
        newend = ", %s_t *%s)" % (c.obj, c.obj)

    if c.findblk(pat, ending = ending, verbose = verbose) < 0:
      print "cannot find [%s], ending [%s]" % (func, ending)
      raise Exception

    #if func.startswith("mdrunner"):
    #  raw_input("New function:\n" + ''.join(c.s[c.begin:c.end+1]) + "\nfor %s\n" % func)

    # add function prefix
    c.s[c.begin] = re.sub(func, c.fmap[func], c.s[c.begin])

    # append obj to the function argument list
    s_end = c.s[c.end].rstrip() + "\n"
    c.s[c.end] = re.sub(oldend, newend, s_end)

    #raw_input("New function:\n" + ''.join(c.s[c.begin:c.end+1]) + "\nfor %s\n" % func)


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
    c.rmcopyrite()
    c.rmionize()
    c.rmtcr()
    c.rmffscan()
    c.rmline("debug_gmx()", doall = True)
    c.rmre()
    c.rmopenmm()
    c.rmetc()

    # remove the declarations of the variable for the integrator
    c.rmline('const gmx_intp_t', wcmt = True, doall = True)
    c.rmline('gmx_integrator_t *func;', off0 = -1, off1 = 3)

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


  @staticmethod
  def getgmxver0():
    ''' get GROMACS version from CMakeLists.txt or configure.ac
        works for up to v4.5 '''

    cfgac = "../../configure.ac"
    cmake = "../../CMakeLists.txt"

    if os.path.exists(cmake): # version 4.5 or later
      for s in open(cmake):
        # we look for something like that
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

    return version, isgmx4, sopenmm


  def getgmxver(c):
    ''' get GROMACS version from configure.ac
        works for up to v4.5 '''

    # this function has been called before
    # no need to repeat it for every instance
    if CCGMX.version > 0: return

    (CCGMX.version, CCGMX.isgmx4, CCGMX.sopenmm) = CCGMX.getgmxver0()

    # verify the instance has the values
    print "gromacs version: %d; v4.0? %s; openmm string: %s" % (
        c.version, bool(c.isgmx4), c.sopenmm)




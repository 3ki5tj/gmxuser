#!/usr/bin/env python

'''
build a self-contained GROMACS engine
* make mdx.c from md.c + runner.c + mdrun.c, see mkmdx_c() (key function)
* make mdxutil.c from force.c + sim_util.c, see mkmdxutil_h()
especially useful for GROMACS 4.5
'''

import re, getopt, os, sys
from cc import CC, savelines, proto2call
from ccgmx import CCGMX

tgtver = (1, 2)  # do atver from 1 to 2, that is md1 and md2
atver = 1
use_at = 1

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
      print "cannot find runner() in %s" % c.fname
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
  c.rmblock("char *desc[]", starter = None, ending = "};", verbose = True)

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
      print "cannot find %s in %s" % (nm, c.fname)
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
      print "cannot find %s in %s" % (nm, c.fname)
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


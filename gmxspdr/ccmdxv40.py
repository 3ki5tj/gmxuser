#!/usr/bin/env python

''' similar to ccmdxv45.py, but for GROMACS 4.0 '''


import re, getopt, os, sys
from ccgmx import CCGMX
from ccutil import tmphastag

def get_bondfree_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  if not tmphastag(txtinp, "bondfree.c"): return ""
  fn = "../gmxlib/bondfree.c"
  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()

  c.rmfunc("int glatnr(int *global_atom_index,int i)")
  c.rmfunc("static int pbc_rvec_sub(const t_pbc *pbc,const rvec xi,")
  c.rmfunc("real morse_bonds(int nbonds,", starter = "/*")
  c.rmfunc("real cubic_bonds(int nbonds,")
  c.rmfunc("real FENE_bonds(int nbonds,")
  c.rmfunc("real harmonic(real kA,real kB,real xA,real xB,real x,")
  c.rmfunc("real bonds(int nbonds,")
  c.rmfunc("real restraint_bonds(int nbonds,")
  c.rmfunc("real polarize(int nbonds,")
  c.rmfunc("real water_pol(int nbonds,")
  c.rmfunc("static real do_1_thole(const rvec xi,const rvec xj,")
  c.rmfunc("real thole_pol(int nbonds,")
  c.rmfunc("real bond_angle(const rvec xi,const rvec xj,const rvec xk,")
  c.rmfunc("real angles(int nbonds,")
  c.rmfunc("real urey_bradley(int nbonds,")
  c.rmfunc("real quartic_angles(int nbonds,")
  c.rmfunc("real dih_angle(const rvec xi,const rvec xj,")
  c.rmfunc("void do_dih_fup(int i,int j,int k,int l,real ddphi,")
  c.rmfunc("real dopdihs(real cpA,real cpB,real phiA,real phiB,int mult,")
  c.rmfunc("static real dopdihs_min(real cpA,real cpB,real phiA,real phiB,int mult,")
  c.rmfunc("real pdihs(int nbonds,")
  c.rmfunc("real idihs(int nbonds,")
  c.rmfunc("real posres(int nbonds,")
  c.rmfunc("static real low_angres(int nbonds,")
  c.rmfunc("real angres(int nbonds,")
  c.rmfunc("real angresz(int nbonds,")
  c.rmfunc("real unimplemented(int nbonds,")
  c.rmfunc("real rbdihs(int nbonds,")

  c.rmfunc("real g96harmonic(real kA,real kB,real xA,real xB,",
      starter = "/*")
  c.rmfunc("real g96bonds(int nbonds,")
  c.rmfunc("real g96bond_angle(const rvec xi,const rvec xj,const rvec xk,")
  c.rmfunc("real g96angles(int nbonds,")
  c.rmfunc("real cross_bond_bond(int nbonds,")
  c.rmfunc("real cross_bond_angle(int nbonds,")
  c.rmfunc("static real bonded_tab(")
  c.rmfunc("real tab_bonds(int nbonds,")
  c.rmfunc("real tab_angles(int nbonds,")
  c.rmfunc("real tab_dihs(int nbonds,")

  # change function name
  c.mutfunc2("calc_bonds", iscall = False)

  return c.s


def get_force_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  if not tmphastag(txtinp, "force.c"): return ""
  if not fn: fn = "../mdlib/force.c"

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()

  c.rmfunc("t_forcerec *mk_forcerec(void)")
  c.rmdf("DEBUG")
  c.rmfunc("static real *mk_nbfp(const gmx_ffparams_t *idef,bool bBHAM)")
  c.rmblock("/* This routine sets fr->solvent_opt", starter = "/*",
      ending = "} solvent_parameters_t;")
  c.rmfunc("check_solvent_cg(const", starter = "static void")
  c.rmfunc("check_solvent(FILE * fp", starter = "static void")
  c.rmfunc("static int *init_cginfo(FILE *fplog")
  c.rmfunc("void set_chargesum(FILE *")
  c.rmfunc("void update_forcerec(FILE *")
  c.rmfunc("void set_avcsixtwelve(FILE *")
  c.rmfunc("static void set_bham_b_max(FILE *")
  c.rmfunc("static void make_nbf_tables(FILE *")
  c.rmfunc("static void count_tables(int")
  c.rmfunc("static bondedtable_t *make_bonded_tables(FILE *fplog,")
  c.rmfunc("void init_forcerec(FILE *fp,")
  c.rmline("#define pr_real(")
  c.rmline("#define pr_int(")
  c.rmline("#define pr_bool(")
  c.rmfunc("void pr_forcerec(FILE *fp,")
  c.rmfunc("void ns(FILE *fp,")
  c.rmfunc("void init_enerdata(FILE *")

  c.mutfunc2("do_force_lowlevel", iscall = False)

  c.rmline("static int timesteps=0;")
  c.rmline("static double t_fnbf=0.0, t_wait=0.0;")
  c.substt("bDoEpot,bSepDVDL", "bSepDVDL")
  c.substt("i,nit,status", "i,status")
  return c.s


def get_simutil_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  if not tmphastag(txtinp, "sim_util.c"): return ""
  if not fn: fn = "../mdlib/sim_util.c"

  c = CCGMX(fn, None, obj, pfx, hdrs = hdrs)
  c.mdcom()

  c.rmline("#define difftime(end,start)")
  c.rmfunc("void print_time(FILE *out,")
  c.rmfunc("time_t print_date_and_time(FILE *fplog,")
  c.rmline("static double runtime=0",
      off0 = -3)
  # the #include statement
  c.rmblock("#ifdef GMX_CRAY_XT3",
      starter = "#ifdef", ending = "#endif")
  # the function start_time()
  c.rmblock("#ifdef GMX_CRAY_XT3",
      starter = "#ifdef", ending = "#endif")
  c.rmfunc("double node_time(void)")
  c.rmfunc("void do_shakefirst(FILE *fplog,")
  c.rmfunc("static void calc_enervirdiff(FILE *fplog,")
  c.rmfunc("void calc_dispcorr(FILE *fplog,")
  c.rmfunc("void do_pbc_first(FILE")
  c.rmfunc("static void low_do_pbc_mtop(FILE")
  c.rmfunc("void do_pbc_first_mtop(FILE")
  c.rmfunc("void do_pbc_mtop(FILE")
  c.rmfunc("void finish_run(FILE")
  c.rmfunc("void init_md(FILE")

  c.mutfunc("sum_forces")
  c.mutfunc("reset_energies")
  c.mutfunc("calc_f_el")
  c.mutfunc("calc_virial")
  c.mutfunc("sum_v")
  c.mutfunc("sum_epot")
  c.mutfunc("print_large_forces")

  # change function names
  c.mutfunc2("do_force_lowlevel", iscall = True)
  c.mutfunc2("do_force", iscall = False)

  c.rmline("tensor virtest;")
  c.subs("int i,j;", "int i;")
  c.rmline("static rvec box_size")
  c.subs(",\s*bBS", "", "bCalcCGCM,\s*bBS")
  c.subs(",\s*cycles_pme", "", "cycles_pme,\s*cycles_force")
  c.subs("e,\s*v,", "v,")
  c.addifdef("matrix\s+boxs;", "GMX_MPI",
      extra = "  bool bBS;\n  float cycles_pme;\n")
  return c.s


def mdrun_c_addopt(c):
  ''' add options -cfg and -mode '''

  # add `-cfg'
  if c.findblk('t_filenm\s*fnm\[\]', ending = "};") >= 0:
    if not c.s[c.end - 1].endswith(","):
      c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
    c.addln(c.end, '    { efMDP, "-cfg",     NULL,       ffOPTRD } /* .cfg file */\n')

  # add `-mode'
  if c.findblk('t_pargs\s*pa\[\]', ending = "};") >= 0:
    if not c.s[c.end - 1].endswith(","):
      c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
    c.addln(c.end, [ '    { "-mode", FALSE, etINT, {&%s},\n' % c.varmode,
          '      "a preparation run" }\n' ])
    if c.s[c.begin].strip().startswith("static"): pfx = "static "
    else: pfx = ""
    c.addln(c.begin, '  %sint %s = 0;\n' % (pfx, c.varmode))


def md_c_changemdrunner(c):
  ''' modify the function mdrunner()
      this function is located in md.c
  '''

  # I. search the declaration of the function mdrunner()
  c.mutfunc2("mdrunner", iscall = False, newend = ", int %s)" % c.varmode,
      pat = "void\s+mdrunner\(")
  runnerbegin = c.begin

  # II. add %obj% declaration
  for j in range(c.end, len(c.s)): # look for a blank line
    if c.s[j].strip() == "":
      break
  else:
    print "cannot find the variable list of the function %s" % c.fmap["mdrunner"]
    raise Exception
  c.addln(j, "    %s_t *%s;\n" % (c.pfx, c.obj) )

  # III. add %pfx%_done() at the very end of the function
  k1, k2 = c.getblockend(c.end, sindent = "", ending = "}", wsp = False)
  c.addln(k2, c.temprepl(
    ' if (%obj% != NULL) %pfx%_done(%obj%, cr);\n',
    True) )

  # IV. add a call to the %pfx%_init() function
  # after the signals are installed
  for i in range(k1, k2):
    # cf. v4.0, md.c, line 1414
    if c.s[i].strip().startswith("signal(SIGUSR1,signal_handler);"):
      l, i = c.getblockend(i, ending = "}", sindent = "  ")
      break
  else:
    print "cannot find the insertion point for %s_init" % c.pfx
    raise Exception

  prjinit = c.temprepl( r'''
  /* initialize project %obj%, cr->duty PP/PME has been assigned */
  %obj% = %pfx%_init(%pfx%_opt2fn("-cfg", nfile, fnm),
      Flags & MD_STARTFROMCPT, state, mtop, inputrec, cr, %MODE%);
  if ((cr->duty & DUTY_PP) && %obj% == NULL) {
    fprintf(stderr, "Failed to initialize %obj%\n");
    exit(1);
  }''', parse = True)
  c.addln(i, prjinit)

  # V. add an handle `obj' to the call to do_md()
  c.substt("integrator[inputrec->eI].func", "do_md", doall = True)
  c.mutfunc2("do_md", iscall = True)

  # VI. place mdrunner() (v4.0, md.c, line 123) after do_md() (v4.0, md.c, line 484)
  k1, k2 = c.getblockend(runnerbegin, sindent = "", ending = "}", wsp = True)
  c.s = c.s[:k1] + c.s[k2+1:] + c.s[k1:k2+1]

  # remove unused variables
  c.rmline("char *gro;")
  c.subs(",\s*nalloc;", ";", "status,\s*nalloc;")
  c.subs("i,\s*m,\s*", "", "i,\s*m,\s*nChargePerturbed")
  c.rmln("real\s+tmpr1,\s*tmpr2;")


def md_c_changedomd(c):
  '''
  change the function md()
  * add a function %pfx%_move()
  * change do_force() to %pfx%_doforce()
  '''

  c.mutfunc2("do_md", iscall = False)

  # add a call %pfx%_move()
  #offset = -1
  #if c.findblk("/\*.*Time for performance", ending = None) < 0:
  #  raise Exception
  #i = c.begin + offset
  i = c.findline("bFirstStep = FALSE;", verbose = True)
  if i < 0: raise Exception
  callmove = c.temprepl(r'''
    /* update %obj% at the end of an MD step */
    %pfx%_move(%obj%, fplog, step,
         bFirstStep, bLastStep, bGStat, bXTC, bNS, enerd,
         state_global, state, &f, top_global, top,
         ir, cr, mdatoms, fr, vsite, shellfc, constr,
         nrnb, wcycle);''' + "\n\n", parse = True)
  c.addln(i, callmove)

  # change the call do_force to %pfx%_doforce
  c.mutfunc2("do_force", iscall = True)
  c.subs(",ed,", ",NULL,", "fp_field,ed,")

  c.rmline("char *grpname;")
  c.rmline("real timestep=0;")


def get_md_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read md.c, and remove junks '''

  if not tmphastag(txtinp, "md.c"): return ""
  if not fn: fn = "md.c"

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()
  md_c_changemdrunner(c)
  md_c_changedomd(c)
  # remove unused variables
  c.substt("gmx_edsam_t ed,", "", doall = True)
  c.substt("mdatoms,nrnb,wcycle,ed,fr,", "mdatoms,nrnb,wcycle,fr,")
  c.rmline("int repl_ex_nst,int repl_ex_seed,")
  c.substt("int repl_ex_nst,int repl_ex_seed,", "", doall = True)
  c.substt("repl_ex_nst,repl_ex_seed,", "", doall = True)
  c.substt("if (do_md == do_md) {", "{")
  return c.s


def get_runner_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  # v4.0 doesn't have a functional runner.c
  # the functions are included in md.c
  return ""


def get_mdrun_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read mdrun.c, and remove junks '''

  if not tmphastag(txtinp, "mdrun.c"): return ""
  if not fn: fn = "mdrun.c"

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()

  # add a parameter `mode' to mdrunner() call
  c.mutfunc2("mdrunner", iscall = True, newend = ", %s);" % c.varmode)
  c.mutfunc("mdrunner")

  # get the static qualifier
  if c.findblk("static const char \*desc", ending = None) >= 0: pfx = "static "
  else: pfx = ""
  c.rmblk("char\s*\*desc", starter = None, ending = "};")
  # simplify the desc
  c.addln(c.begin, '  %sconst char *desc[] = {"self-contained GROMACS"};\n' % pfx)
  mdrun_c_addopt(c)
  c.substt("ed,repl_ex_nst,repl_ex_seed,", "", doall = True)

  c.rmline("bool bPPPME = FALSE;")
  c.rmline("bool bCart = FALSE;")
  c.rmline("bool HaveCheckpoint;")
  c.substt(",*fptest;", ";")
  return c.s



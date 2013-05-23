#!/usr/bin/env python

'''
  specific changes and tweaks to various GROMACS v4.5 files
  no class is defined in this file
  * get_bondfree_c()
  * get_force_c()
  * get_simutil_c()
  * get_md_c()
  * get_runner_c()
  * get_mdrun_c()
'''

import re, getopt, os, sys
from ccgmx import CCGMX
from ccutil import tmphastag, getgmx

def get_bondfree_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read bondfree.c and remove junks '''

  if not tmphastag(txtinp, "bondfree.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "gmxlib")
    else:
      fn = os.path.join(getgmx.root, "src", "gromacs", "gmxlib")
    fn = os.path.join(fn, "bondfree.c")

  c = CCGMX(fn, None, obj, pfx, hdrs = hdrs)
  c.mdcom()

  # I. remove unnecessary functions
  c.rmdf("SIMD_BONDEDS")
  c.rmfunc("const int cmap_coeff_matrix[] = {")
  c.rmfunc("int glatnr(int *global_atom_index,")
  c.rmfunc("static int pbc_rvec_sub(const t_pbc *pbc,")
  c.rmfunc("real morse_bonds(int nbonds,", starter = "/*")
  c.rmfunc("real cubic_bonds(int nbonds,")
  c.rmfunc("real FENE_bonds(int nbonds,")
  c.rmfunc("real harmonic(real kA,")
  c.rmfunc("real anharm_polarize(int nbonds,")
  c.rmfunc("real bonds(int nbonds,")
  c.rmfunc("real restraint_bonds(int nbonds,")
  c.rmfunc("real polarize(int nbonds,")
  c.rmfunc("real water_pol(int nbonds,")
  c.rmfunc("static real do_1_thole(const rvec xi,")
  c.rmfunc("real thole_pol(int nbonds,")
  c.rmfunc("real bond_angle(const rvec xi,")
  c.rmfunc("real angles(int nbonds,")
  c.rmfunc("real linear_angles(int nbonds,")
  c.rmfunc("real urey_bradley(int nbonds,")
  c.rmfunc("real quartic_angles(int nbonds,")
  c.rmfunc("real dih_angle(const rvec xi,")
  c.rmfunc("void do_dih_fup(int i,")
  c.rmfunc("do_dih_fup_noshiftf(int i,", starter = "static void")
  c.rmfunc("do_dih_fup_noshiftf_precalc(int i,", starter = "static")
  c.rmfunc("real dopdihs(real cpA,")
  c.rmfunc("dopdihs_noener(real cpA,", starter = "static void")
  c.rmfunc("static real dopdihs_min(real cpA,")
  c.rmfunc("dopdihs_mdphi(real cpA,", starter = "static void")
  c.rmfunc("real pdihs(int nbonds,")
  c.rmfunc("void make_dp_periodic(real *dp)")
  c.rmfunc("pdihs_noener(int nbonds,", starter = "static void")
  c.rmfunc("real idihs(int nbonds,")
  c.rmfunc("real posres(int nbonds,")
  c.rmfunc("static void posres_dx(const rvec x,")
  c.rmfunc("real fbposres(int nbonds,")
  c.rmfunc("static real low_angres(int nbonds,")
  c.rmfunc("real angres(int nbonds,")
  c.rmfunc("real angresz(int nbonds,")
  c.rmfunc("real dihres(int nbonds,")
  c.rmfunc("real unimplemented(int nbonds,")
  c.rmfunc("real rbdihs(int nbonds,")
  c.rmfunc("int cmap_setup_grid_index(int ip,")

  # we need a prototype of cmap_dihs
  if c.findblock("real cmap_dihs(int nbonds", ending = ")", wsp = False) < 0:
    print "no %s in %s" % (cmap_dihs, c.fn)
    raise Exception
  proto = c.s[c.begin : c.end + 1]
  proto[-1] = proto[-1].rstrip() + ";\n\n"
  c.rmfunc("real cmap_dihs(int nbonds,")
  c.addln(c.begin, proto)

  c.rmfunc("real g96harmonic(real kA,",
      starter = "/*")
  c.rmfunc("real g96bonds(int nbonds,")
  c.rmfunc("real g96bond_angle(const rvec xi,")
  c.rmfunc("real g96angles(int nbonds,")
  c.rmfunc("real cross_bond_bond(int nbonds,")
  c.rmfunc("real cross_bond_angle(int nbonds,")
  c.rmfunc("static real bonded_tab(")
  c.rmfunc("real tab_bonds(int nbonds,")
  c.rmfunc("real tab_angles(int nbonds,")
  c.rmfunc("real tab_dihs(int nbonds,")
  c.rmfunc("calc_bonded_reduction_mask(const t_idef *idef,",
      starter = "static")
  c.rmfunc("void init_bonded_thread_force_reduction(t_forcerec")
  c.rmfunc("static real calc_one_bond_foreign(FILE *fplog,")
  c.rmfunc("void calc_bonds_lambda(FILE *fplog,")

  c.mutfunc("reduce_thread_forces", pfx = c.pfx + "cb")

  # change function name
  c.mutfunc2("calc_bonds", iscall = False)
  return c.s


def get_force_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' real force.c and remove junks '''

  if not tmphastag(txtinp, "force.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "mdlib")
    else:
      fn = os.path.join(getgmx.root, "src", "gromacs", "mdlib")
    fn = os.path.join(fn, "force.c")

  c = CCGMX(fn, None, obj, pfx, hdrs = hdrs)
  c.mdcom()

  # I. remove unnecessary functions
  # cf. v4.5, force.c, lines 68-115
  c.rmfunc("void ns(FILE *fp,")
  # cf. v4.5, force.c, lines 585-615
  c.rmfunc("void init_enerdata(")
  # cf. v4.5, force.c lines 616-630
  c.rmfunc("void destroy_enerdata(gmx_enerdata_t *enerd)")
  # cf. v4.5, force.c lines 631-642
  c.rmfunc("static real sum_v(int n,real v[])")
  # cf. v4.5, force.c lines 643-673
  c.rmfunc("void sum_epot(t_grpopts *opts,gmx_enerdata_t *enerd")
  # cf. v4.5, force.c lines 674-711
  c.rmfunc("void sum_dhdl(gmx_enerdata_t *enerd,double lambda,")
  # cf. v4.5, force.c lines 712-749
  c.rmfunc("void reset_enerdata(t_grpopts *opts,")

  # II. remove unused variables
  c.rmline("real dvdgb")
  c.substt("bDoEpot,", "")

  c.mutfunc("reduce_thread_forces", pfx = c.pfx + "fl")

  # III. change function name
  c.mutfunc2("do_force_lowlevel", iscall = False)

  # if we have to change bondfree.c, then
  # change the call calc_bonds()
  if tmphastag(txtinp, "bondfree.c"):
    c.mutfunc2("calc_bonds", iscall = True)
  return c.s



def get_simutil_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read sim_util.c, and remove junks '''

  if not tmphastag(txtinp, "sim_util.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "mdlib")
    else:
      fn = os.path.join(getgmx.root, "src", "gromacs", "mdlib")
    fn = os.path.join(fn, "sim_util.c")

  c = CCGMX(fn, None, obj, pfx, hdrs = hdrs)
  c.mdcom()

  # I. remove unecessary functions
  # cf. v4.5, sim_util.c, lines 105-109
  c.rmfunc("gmx_ctime_r(const time_t *clock,char *buf, int n);", starter = "/*")
  # cf. v4.5, sim_util.c, lines 110-131
  c.rmfunc("gmx_gettime()", starter = "double")
  # cf. v4.5, sim_util.c, lines 132-133
  c.rmline("#define difftime(end,start)")
  # cf. v4.5, sim_util.c, lines 134-188
  c.rmfunc("void print_time(FILE *out,gmx_runtime_t *runtime,gmx_large_int_t step,")
  # cf. v4.5, sim_util.c, lines 189-219
  c.rmfunc("static double set_proctime(gmx_runtime_t *runtime)",
      starter = "#ifdef NO_CLOCK")
  # cf. v4.5, sim_util.c, lines 220-230
  c.rmfunc("void runtime_start(gmx_runtime_t")
  # cf. v4.5, sim_util.c, lines 231-241
  c.rmfunc("void runtime_end(gmx_runtime_t")
  # cf. v4.5, sim_util.c, lines 242-246
  c.rmfunc("void runtime_upd_proc(gmx_runtime_t")
  # cf. v4.5, sim_util.c, lines 247-276
  c.rmfunc("void print_date_and_time(FILE")
  # cf. v4.5, sim_util.c, lines 916-1002
  c.rmfunc("void do_constrain_first(FILE")
  # cf. v4.5, sim_util.c, lines 1003-1132
  c.rmfunc("void calc_enervirdiff(FILE")
  # cf. v4.5, sim_util.c, lines 1133-1243
  c.rmfunc("void calc_dispcorr(FILE")
  # cf. v4.5, sim_util.c, lines 1243-1267
  c.rmfunc("void do_pbc_first(FILE")
  # cf. v4.5, sim_util.c, lines 1268-1308
  c.rmfunc("static void low_do_pbc_mtop(FILE")
  # cf. v4.5, sim_util.c, lines 1309-1314
  c.rmfunc("void do_pbc_first_mtop(FILE")
  # cf. v4.5, sim_util.c, lines 1315-1320
  c.rmfunc("void do_pbc_mtop(FILE")
  # cf. v4.5, sim_util.c, lines 1321-1423
  c.rmfunc("void finish_run(FILE")
  # cf. v4.5, sim_util.c, lines 1424-1503
  c.rmfunc("void init_md(FILE")

  # II. add funtion prefix
  # some static functions have to be kept here
  c.mutfunc("sum_forces")
  c.mutfunc("calc_f_el")
  c.mutfunc("calc_virial")
  c.mutfunc("print_large_forces")

  # change function names
  c.mutfunc2("do_force", iscall = False)
  if getgmx()[0] >= 40600:
    c.mutfunc2("do_force_cutsVERLET", iscall = False)
    c.mutfunc2("do_force_cutsVERLET", iscall = True)
    c.mutfunc2("do_force_cutsGROUP",  iscall = False)
    c.mutfunc2("do_force_cutsGROUP",  iscall = True)

  # if we have to change force.c, then
  # change the call do_force_lowlevel()
  if tmphastag(txtinp, "force.c"):
    c.mutfunc2("do_force_lowlevel", iscall = True)

  # III. remove uncessary blocks
  c.rmblk("do_flood\(fplog", starter = "if (ed)", ending = "}")

  # IV. unused variables
  c.rmline("tensor virtest;")
  c.substt("int i,j;", "int i;");
  c.substt(",bBS;", ";", "bCalcCGCM,bBS");
  c.addifdef("matrix\s+boxs;", "GMX_MPI", extra = "    gmx_bool bBS;\n")
  return c.s



def get_md_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read md.c, and remove junks '''

  if not tmphastag(txtinp, "md.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "kernel")
    else:
      fn = os.path.join(getgmx.root, "src", "programs", "mdrun")
    fn = os.path.join(fn, "md.c")

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()

  # change the name of do_md()
  c.mutfunc2("do_md", iscall = False)

  # add a call %PFX%_move()
  #offset = 1
  #if c.findblk("/\*.*END PREPARING EDR OUTPUT", ending = None) < 0:
  #  raise Exception
  i = c.findline("bFirstStep = FALSE;", verbose = True)
  if i < 0: raise Exception
  callmove = c.temprepl(r'''
        /* update %OBJ% */
        if (0 != %PFX%_move(%OBJ%, enerd,
             step, bFirstStep, bLastStep, bGStat,
             mdof_flags & MDOF_XTC, bNS, cr) ) {
          exit(1);
        }

''', parse = True)
  c.addln(i, callmove)

  # if we also need to modify sim_util.c, then
  # change the call do_force()
  if tmphastag(txtinp, "sim_util.c"):
    c.mutfunc2("do_force", iscall = True)

  # remove useless blocks
  # v4.5, md.c, lines 526-534, 795-804, 934-946
  #c.rmblk("if\s+(bFFscan)", ending = "}", doall = True)

  # remove unused variable
  c.rmline("int iter_i;")
  c.rmline("real timestep=0;")
  c.rmline("tensor tmpvir;")
  c.rmline("char *grpname;")
  c.rmline("atom_id *grpindex=NULL;")
  c.substt(",gnx=0,ii;", ";", "gnx=0,ii;")
  c.substt("*scale_tot,pcoupl_mu,M,ebox;", "pcoupl_mu,M;", "matrix *scale_tot,pcoupl_mu,M,ebox;")
  # temp0 has been fixed in v4.6
  c.substt("temp0,mu_aver=0,", "", "temp0,mu_aver=0,dvdl;")
  c.rmline("temp0 = enerd->term[F_TEMP];")
  c.substt("fom,oldfom,", "", "fom,oldfom,veta_save")
  c.substt("boxcopy={{0}},", "")
  c.substt("*xcopy=NULL,*vcopy=NULL,", "")
  c.rmline("gmx_bool bAppend;")
  c.rmline("bAppend = (Flags & MD_APPENDFILES);")
  return c.s



def get_runner_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read runner.c, and remove junks '''

  if not tmphastag(txtinp, "runner.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "kernel")
    else:
      fn = os.path.join(getgmx.root, "src", "programs", "mdrun")
    fn = os.path.join(fn, "runner.c")

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()

  # add `varmode' to the declaration of mdrunner_arglist, cf. v4.5, runner.c, line 113
  if c.findblk("struct\s*mdrunner_arglist") >= 0:
    c.addln(c.end-1, "    int %s;\n" % c.varmode)
  # add `varmode' to the definition mdrunner_start_threads, cf. v4.5, runner.c, line 182
  if c.findblk("t_commrec\s*\*mdrunner_start_threads", ending = ")") >= 0:
    c.s[c.end] = c.s[c.end].rstrip()[:-1] + ", int %s)\n" % c.varmode
  # add `varmode' to the call of mdrunner_start_threads(), cf. v4.5, runner.c, line 398
  if c.findblk('\w\s*\=\s*mdrunner_start_threads\(', ending = ");") >= 0:
    c.s[c.end] = re.sub(r"Flags\)", "Flags, %s)" % c.varmode, c.s[c.end])
  # append `varmode' to the mda assignment, cf. v4.5, runner.c, line 235
  if c.findblk('mda\-\>Flags=Flags;', ending = None) >= 0:
    c.addln(c.end+1, '    mda->%s = %s;\n' % (c.varmode, c.varmode))

  # I. search the declaration of the function mdrunner()
  c.mutfunc2("mdrunner", iscall = False, newend = ", int %s)" % c.varmode,
      pat = "int\s+mdrunner\(")

  # II. add %OBJ% declaration
  for j in range(c.end, len(c.s)): # look for a blank line
    if c.s[j].strip() == "":
      break
  else:
    print "cannot find the variable list of the function %s" % c.fmap["mdrunner"]
    raise Exception
  c.addln(j, "    %s_t *%s;\n" % (c.pfx, c.obj) )

  # III. save a new prototype of the function
  proto = c.s[c.begin : c.end + 1]
  #raw_input("The prototype\n" + ''.join(proto))
  last = len(proto) - 1
  proto[last] = proto[last].rstrip() + ";\n"
  # add a comment to the prototype
  proto = ["\n", "/* declare runner() before mdrunner_start_fn() uses it */\n",
      "static\n"] + proto

  # IV. add %PFX%_done() at the very end of the function
  k1, k2 = c.getblockend(c.end, sindent = "", ending = "}", wsp = False)
  k2 -= 1 # skip the statement "return rc;" cf. v4.5, runner.c, line 892
  c.addln(k2, c.temprepl('''    if (%OBJ% != NULL)
      %PFX%_done(%OBJ%);
''', True) )

  # V. add a call to the %PFX%_init() function
  # after the signals are installed
  for i in range(k1, k2):
    # cf. v4.5, runner.c, line 797
    if c.s[i].strip().startswith("signal_handler_install();"):
      # search for the signal finishing
      l, i = c.getblockend(i, ending = "}", sindent = "    ")
      break
  else:
    print "cannot find insersion point for %s_init" % c.pfx
    raise Exception

  callinit = c.temprepl(r'''
    /* initialize project %OBJ%, cr->duty PP/PME has been assigned */
    %OBJ% = %PFX%_init(%PFX%_opt2fn("-cfg", nfile, fnm),
        Flags & MD_STARTFROMCPT, mtop, inputrec, cr, %MODE%);
    if ((cr->duty & DUTY_PP) && %OBJ% == NULL) {
      fprintf(stderr, "Failed to initialize %OBJ%\n");
      exit(1);
    }''', parse = True)
  c.addln(i, callinit)

  # VI. if md.c has to be changed, then
  # change the call do_md()
  # NOTE: must compile with md.c, if %md.c% is not included
  if tmphastag(txtinp, "md.c"):
    c.substt("integrator[inputrec->eI].func", "do_md", doall = True)
    c.mutfunc2("do_md", iscall = True)
    # remove the declarations of the variable for the integrator
    c.rmline('const gmx_intp_t', wcmt = True, doall = True)
    c.rmline('gmx_integrator_t *func;', off0 = -1, off1 = 3)


  # VII. change mdrunner()
  c.mutfunc2("mdrunner", iscall = True, newend = ", mc.%s);" % c.varmode)
  c.mutfunc("mdrunner")

  # VIII. insert the prototype of mdrunner() after a bunch of #include
  # at the beginning of the file
  i = len(c.s) - 1
  while i >= 1:
    if c.s[i].startswith("#include") and c.s[i-1].startswith("#include"):
      break
    i -= 1
  c.addln(i+1, proto)

  # clear the ``if (do_md == do_md\n)'' mess, after replacing integrator
  # cf. v4.5, runner.c, line 785
  c.rmline("if (do_md == do_md", off1 = 2)

  # remove unused variables
  c.addifdef("t_commrec\s+\*cr_old\s*=\s*cr;", "GMX_THREADS")
  c.addifdef("int\s+nthreads=1;", "GMX_THREADS")
  c.rmline("int list;")
  c.rmline("char *gro;")
  c.substt(",nalloc;", ";", "status,nalloc;")
  c.substt("i,m,", "", "i,m,nChargePerturbed")
  c.rmline("real tmpr1,tmpr2;")
  c.rmline("double nodetime=0,realtime;")
  return c.s



def mdrun_c_addopt(c):
  ''' add options -cfg and -mode '''

  # add `-cfg'
  if c.findblk('t_filenm\s*fnm\[\]', ending = "};") >= 0:
    if not c.s[c.end - 1].endswith(","):
      c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
    c.addln(c.end, '    { efMDP, "-cfg",     NULL,       ffOPTRD } /* .cfg file */\n')
  else:
    print "cannot find fnm section!"
    open("dump.txt", "w").writelines(c.s)
    raise Exception

  # add `-mode'
  if c.findblk('t_pargs\s*pa\[\]', ending = "};") >= 0:
    if not c.s[c.end - 1].endswith(","):
      c.s[c.end - 1] = c.s[c.end - 1].rstrip() + ",\n"
    c.addln(c.end, [ '    { "-mode", FALSE, etINT, {&%s},\n' % c.varmode,
          '      "a preparation run" }\n' ])
    if c.s[c.begin].strip().startswith("static"): pfx = "static "
    else: pfx = ""
    c.addln(c.begin, '  %sint %s = 0;\n' % (pfx, c.varmode))



def get_mdrun_c(txtinp, obj, pfx, fn = None, hdrs = {}):
  ''' read mdrun.c, and remove junks '''

  if not tmphastag(txtinp, "mdrun.c"): return ""
  if not fn:
    if getgmx()[0] < 40600:
      fn = os.path.join(getgmx.root, "src", "kernel")
    else:
      fn = os.path.join(getgmx.root, "src", "programs", "mdrun")
    fn = os.path.join(fn, "mdrun.c")

  c = CCGMX(fn, None, obj, pfx, hdrs)
  c.mdcom()
  if c.findblk("static const char \*desc", ending = None) >= 0: pfx = "static "
  else: pfx = ""
  c.rmblk("char\s*\*desc", starter = None, ending = "};")
  # simplify the desc
  c.addln(c.begin, '  %sconst char *desc[] = {"self-contained GROMACS"};\n' % pfx)
  mdrun_c_addopt(c)

  c.mutfunc2("mdrunner", iscall = True, newend = ", %s);" % c.varmode)
  c.mutfunc("mdrunner")

  # remove unused variables
  c.rmline("gmx_bool bPPPME = FALSE;")
  c.rmline("gmx_bool bCart = FALSE;")
  c.rmline("gmx_edsam_t ed;")
  c.rmln("      int i;")
  c.substt(",*fptest;", ";")

  c.rmblock("Comment this in to do fexist",
      starter = "/*", ending = "*/")
  c.rmblock("PCA_Flags |= PCA_NOT_READ_NODE;",
      starter = "/*", ending = "*/")
  return c.s



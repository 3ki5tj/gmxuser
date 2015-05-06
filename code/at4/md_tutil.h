/* =========================================================================================

         GROMACS integration

 * ========================================================================================= */

#include "md_tcore.h"

static __inline const char * yesno(int x) { return (x) ? "yes" : "no"; }

/* update the current scale for all nodes */
static void atgmx_updscl(at_t *at, t_commrec *cr)
{
  int ismaster = SIMMASTER(cr);

  /* call the single processor version */
  if (ismaster) at_updscl(at);

  /* here we tell all PP nodes the new scale(s) */
  if (PAR(cr)) {
#if AT_VER == 2
    gmx_bcast(AT_ETOT * sizeof(at->scale[0]), at->scale, cr);
    /* currently PME-only nodes are too stupid to handle at->scale */

    /* calculate for non-masters the current temperature */
    if (!ismaster) {
      at->beta = at->scale[AT_EBG] * at_beta_T(at, at->T0);
      at->th_f = at->scale[AT_EFG]/at->scale[AT_EBG] - 1.0;
    }
#elif AT_VER == 1
    gmx_bcast(sizeof(at->scale), &at->scale, cr);
#endif
  }
}

#if AT_VER == 2

typedef int (*fexcl_t)(void *, const int *, int);

/* Add an atom set (ids[0], ids[1], ids[2], ...) to the exclusion list,
 * we use this list to avoid unwanted atom set (e.g., those of PRO/GLY)
 *
 * The implementation of exclusion list is based on atoms' global indices,
 * and therefore is not as neat as relying on a particular itype.
 * Alternative implementation would require we add a new topology entry
 * at runtime to distinguish an unwanted set from a normal set of atoms.
 * (runtime addition makes the program more flexible)
 * However, this can be more messy.
 * As long as the exclusion list is short, the current approach should be okay.
 * */
static int atgmx_excl(void *spb, fexcl_t fexcl,
    gmx_mtop_t *mtop, int imt, int *ids, int n)
{
  int im, imb, mbsize, ag, ret, j;
  gmx_molblock_t *mb;

  for (ag = 0, imb = 0; imb < mtop->nmolblock; imb++) {
    mb = mtop->molblock + imb;
    mbsize = mb->natoms_mol;
    /* found a match for moltype id */
    if (mb->type == imt) {
      /* remove all instances of the moltype */
      for (im = 0; im < mb->nmol; im++) {
        ret = (*fexcl)(spb, ids, ag + im * mbsize);
        if (ret != 0) {
          fprintf(stderr, "cannot exclude ");
          for (j = 0; j < n; j++) fprintf(stderr, "%d, ",  ids[j]);
          fprintf(stderr, "\n");
          exit(1);
        }
      }
    }
    ag += mbsize * mb->nmol;
  }
  return 0;
}

/* For a single spb, given by the names of the involved atoms,
 * return the global index, itype, for a particular interaction,
 * such as "C-N-CA-C".
 * itype is used to identify all instances of the interaction type.
 * We also put unwanted exceptions (such as dihedrals involving a GLY
 * or PRO) to the exclusion list associated with the spb
 * if FSPB_NOCATER is specified, the phi angle with ACE (N-terminal)
 * and the psi angle with NH2 (C-terminal) are excluded
 *
 *  mtop      :   gromacs global topology
 *  anames    :   an array of the involved atoms, e.g., {"C", "N", "CA", "C"}
 *  nratoms   :   # of involved atoms
 *  functype  :   functional type of the interaction, e.g., F_RBDIHS
 */
#define FSPB_VERBOSE  0x1
#define FSPB_NOCAGLY  0x2
#define FSPB_NOCAPRO  0x4
#define FSPB_NOCATER  0x8

#if defined(GMXVERSION) && (GMXVERSION >= 40099)
#define ATOM2RES(atoms, aid) atoms.atom[aid].resind
#define RES2NAME(atoms, rid) atoms.resinfo[rid].name[0]
#else
#define ATOM2RES(atoms, aid) atoms.atom[aid].resnr
#define RES2NAME(atoms, rid) atoms.resname[rid][0]
#endif

/* search interaction list in `mtop'
 * for atom names in `anames' and `functype' (F_RBDIHS)
 * `flags' can be combination of FSPB_NOCAGLY, FSPB_NOCAPRO, FSPB_NOCATER ...
 * if `onlyres' > 0, all residues except (onlyres-1) are excluded */
static int atgmx_findspb(void *spb, fexcl_t fexcl, gmx_mtop_t *mtop,
    const char **anames, int nratoms, int functype, unsigned flags,
    int onlyres)
{
  t_ilist *il;
  gmx_moltype_t *mt;
  int i, j, imt, *ia, natoms;
  int functype2, itype, itype_prev = -1, skipthis;
  char *allatoms, *atlist = NULL, strtmp[32];
  int resid, resid_key;
  char *resnm, *resnm_key;
  int hasace;
  int verbose = flags & FSPB_VERBOSE;
  int nocagly = flags & FSPB_NOCAGLY;
  int nocapro = flags & FSPB_NOCAPRO;
  int nocater = flags & FSPB_NOCATER;

  fprintf(stderr, "Searching for special %s... GLY: %s, PRO: %s, TER: %s\n",
      interaction_function[functype].longname, yesno(!nocagly), yesno(!nocapro), yesno(!nocater));

  atlist = ssdup("");

  natoms = interaction_function[functype].nratoms;
  die_if (natoms != nratoms,
    "# of atoms for %s is %d, not %d\n",
      interaction_function[functype].longname, natoms, nratoms);

  /* join names of all the involved atoms to `allatoms' */
  allatoms = ssdup("|");
  for (i = 0; i < natoms; i++) {
    sscat(allatoms, anames[i]);
    sscat(allatoms, "|");
  }

  /* loop over moltypes to search a interaction that matches `allatoms' */
  for (imt = 0; imt < mtop->nmoltype; imt++) {
    mt = mtop->moltype + imt;

    /* grap the interaction list */
    il = &mt->ilist[functype];

    /* check N-cap ACE */
    resnm = RES2NAME(mt->atoms, 0);
    hasace = (strcmp(resnm, "ACE") == 0);

    for (i = 0; i < il->nr; ) {
      /* The interaction list is arranged as
       * {typeA, atomA0, atomA1,
       *  typeB, atomB0, atomB1, atomB2, atomB3,
       *  ... }
       * the leading number, `itype', is an interaction index, like 305,
       * that corresponds to a particular set of atoms such as N-C-CA-N.
       * This `itype' is different from the general functype, which just
       * something like F_RBDIHS. */
      itype = il->iatoms[i++];  /* the leading interaction id */
      ia = il->iatoms + i; /* the following members are atom indices */
      i += natoms;

      /* create a combined atom-list */
      for (atlist[0] = '\0', j = 0; j < natoms; j++) {
        sprintf(strtmp, "%4d", ia[j] + 1);
        if (j > 0) sscat(atlist, "|");
        sscat(atlist, strtmp);
      }

      /* double check if the corresponding functype is right or not */
      functype2 = mtop->ffparams.functype[itype];
      die_if (functype2 != functype,
        "atgmx_findspb: bad functype %d vs. %d\n", functype2, functype);

      /* compare atom names, we do require they share the same case */
      for (j = 0; j < natoms; j++)
        if (strcmp(mt->atoms.atomname[ia[j]][0], anames[j]) != 0)
          break;

      /* atom names do not match, goto the next interaction */
      if (j != natoms) continue;

      /* determine the residue name of this interaction */
      /* a. use the residue from the first atom by default */
      resid_key = ATOM2RES(mt->atoms, ia[0]);
      resnm_key = RES2NAME(mt->atoms, resid_key);

      /* b. replace the residue name by that of an involved CA atom, if any,
       *    also decide if to add the interaction to the exclusion list */
      skipthis = 0;
      for (j = 0; j < natoms; j++) {
        /* trace it back to the residue ID */
        resid = ATOM2RES(mt->atoms, ia[j]);
        resnm = RES2NAME(mt->atoms, resid);

        if (nocater) {
          if (strcmp(resnm, "ACE") == 0
           || strcmp(resnm, "NH2") == 0) {
            resid_key = resid;
            resnm_key = resnm;
            skipthis = 1;
            atgmx_excl(spb, fexcl, mtop, imt, ia, natoms);
            break;
          }
        }

        if (strcmp(anames[j], "CA") == 0) { /* check a non-CA atom */
          resid_key = resid;
          resnm_key = resnm;
          if ( (nocagly && strstr(resnm, "GLY") != NULL)
           ||  (nocapro && strstr(resnm, "PRO") != NULL) ) {
            skipthis = 1;
            /* search all molblocks, and add all instances of
             * this moltype to the exclusion list */
            atgmx_excl(spb, fexcl, mtop, imt, ia, natoms);
            break;
          } else if (onlyres > 0 && resid - hasace != onlyres -1) {
            skipthis = 1;
            atgmx_excl(spb, fexcl, mtop, imt, ia, natoms);
            break;
          }
        }
      }

      if (itype_prev < 0 || verbose)
        fprintf(stderr, "%3d:%4s %s %s for %s (%d)\n",
            resid_key, resnm_key, (skipthis ? "! SKIP ": "  found"), atlist, allatoms, itype);

      if (!skipthis) {
        if (itype_prev < 0) { /* no existing itype */
          itype_prev = itype;
        } else {
          die_if (itype != itype_prev, /* interaction itype conflict */
            "double interaction types=%d[new(%s)] %d(old) for %s\n",
             itype, atlist, itype_prev, allatoms);
        }
      }
    } /* loop over interactions */
  } /* loop over moltypes */

  ssdel(allatoms);
  ssdel(atlist);
  return itype_prev;
}

/* for each spb, search its corresponding interaction type
 * every node calls this */
static int atgmx_initspb(at_t *at, gmx_mtop_t *mtop, t_commrec *cr)
{
  int i, itype, err = 0;
  spb_t *spb;

  fflush(stderr); /* to improve the output */
  if (PAR(cr)) gmx_barrier(cr);

  if (SIMMASTER(cr)) {
    unsigned flags = FSPB_VERBOSE;

    /* assign appropriate flags for skipping GLY and/or PRO */
    if (!at->spb_dogly)
      flags |= FSPB_NOCAGLY;
    if (!at->spb_dopro)
      flags |= FSPB_NOCAPRO;
    if (!at->spb_doter)
      flags |= FSPB_NOCATER;

    /* loop over all spbs */
    for (i = 0; !err && i < at->spbs->cnt; i++) {
      spb = at->spbs->arr + i;
      /* search topology to extract for interaction type for
       * a particular R.B. dihedral, quit if the corresponding
       * spb is not found */
      itype = atgmx_findspb(spb, (fexcl_t) spb_excl, mtop,
          (const char **) spb->atoms, spb->natoms,
          F_RBDIHS, flags, at->spb_onlyres);
      if (itype < 0) {
        fprintf(stderr, "failed to find an itype for spb %d\n", i);
        err = 1;
        break;
      }
      spb->type = itype;
      fprintf(stderr, "spb %d, moltype itype: %d, %d in exclusion list\n",
        i, itype, spb->exclcnt);
    }
  }

  fflush(stderr); /* allow spb info to finish */
  /* if there is an error on the master, every node knows and quits */
  if (PAR(cr)) {
    gmx_bcast(sizeof(err), &err, cr);
  }
  if (err) {
    fprintf(stderr, "spb: error occurred in the master branch\n");
    return -1;
  }

  if (PAR(cr)) {
    if ( !(cr->duty & DUTY_PP) ) {
      at->spbs = NULL; /* make sure it is NULL for a PME-only node */
    } else {
      /* we assume that SIMMASTER(cr) is the MASTER(cr)
       * Note, we will allocate spaces for spbs, even in case that spbs->cnt
       * is zero, just to make the interface of accessing spbs consistent.
       * The memory space for the master is already allocated */
      if ( !MASTER(cr) ) { /* allocate space for nonmaster */
        xnew(at->spbs, 1);
      }

#ifdef GMX_MPI
      /* now that the space for every node is allocated,
       * we share information saved in the master.
       * spbs->rank and spbs->comm will be initialized here and forever */
      if (0 != spbs_initmpi(at->spbs, cr->mpi_comm_mygroup)) {
        fprintf(stderr, "NODE %d: error during spreading spb info\n",
          cr->sim_nodeid);
        return -1;
      }
#endif
      gmx_barrier(cr);
    }
  }
  return 0;
}

/* search backbone pairs
 * every node calls this */
static int atgmx_initbb(at_t *at, gmx_mtop_t *mtop, t_commrec *cr)
{
  int i, itype, err = 0;
  bb_t *bb;

  fflush(stderr); /* to improve the output */
  if (PAR(cr)) gmx_barrier(cr);

  if (SIMMASTER(cr)) {
    unsigned flags = FSPB_VERBOSE;
    const char *atoms[5] = {"C", "N", "CA", "C", "N"};

    /* assign appropriate flags for skipping GLY and/or PRO */
    if (!at->bb_dogly)
      flags |= FSPB_NOCAGLY;
    if (!at->bb_dopro)
      flags |= FSPB_NOCAPRO;

    /* loop over all bbs */
    for (i = 0; !err && i < at->bbs->cnt; i++) {
      bb = at->bbs->arr + i;
#if GMXVERSION >= 40099
      itype = atgmx_findspb(bb, (fexcl_t) bb_excl, mtop, atoms, BB_ATOMS,
          F_CMAP, flags, at->bb_onlyres);
#else
      (void) atoms;
      itype = -(mtop != NULL) - 1; /* raise an error in version 4 */
#endif
      if (itype < 0) {
        fprintf(stderr, "bb: failed to find an itype for bb %d\n", i);
        err = 1;
        break;
      }
      bb->type = itype;
      fprintf(stderr, "bb %d: moltype itype: %d, %d in exclusion list\n",
        i, itype, bb->exclcnt);
    }
  }
  fflush(stderr);
  /* if there is an error on the master, every node knows and quits */
  if (PAR(cr)) {
    gmx_bcast(sizeof(err), &err, cr);
  }
  if (err) {
    fprintf(stderr, "bb: error occurred in the master branch\n");
    return -1;
  }

  if (PAR(cr)) {
    if ( !(cr->duty & DUTY_PP) ) {
      at->bbs = NULL; /* make sure it is NULL for a PME-only node */
    } else {
      /* we assume that SIMMASTER(cr) is the MASTER(cr)
       * Note, we will allocate spaces for spbs, even in case that bbs->cnt
       * is zero, just to make the interface of accessing spbs consistent.
       * The memory space for the master is already allocated */
      if ( !MASTER(cr) ) { /* allocate space for nonmaster */
        if ((at->bbs = calloc(1, sizeof(*at->bbs))) == NULL) {
          fprintf(stderr, "node %d: no memory for struct bbs_t.\n",
            cr->sim_nodeid);
          return -1;
        }
      }

#ifdef GMX_MPI
      /* now that the space for every node is allocated,
       * we share information saved in the master.
       * spbs->rank and spbs->comm will be initialized here and forever */
      if (0 != bbs_initmpi(at->bbs, cr->mpi_comm_mygroup)) {
        fprintf(stderr, "NODE %d: error during spreading spb info\n",
          cr->sim_nodeid);
        return -1;
      }
#endif
      gmx_barrier(cr);
    }
  }
  return 0;
}
#endif

/* Classify energies to H0 and H1
 * and calculate Ea, H0, and H1
 *
 * To simiplify the logic, I now assume that global_stat() has been called,
 * so everything is correct now. */
static void atgmx_sum(at_t *at, real *eterm, t_commrec *cr,
    llong_t step, int dirty)
{
  int ismaster = MASTER(cr);
#if AT_VER == 2
  int verbose;
  static int once;

  die_if(!(cr->duty & DUTY_PP), /* Only a node does MD calls this */
      "pme node %d calls atgmx_sum\n", cr->nodeid);

  verbose = ismaster && (once < 3 || once % 100000 == 0);
  if (verbose) fprintf(stderr, "\n");

  die_if (dirty,
    "node = %d, step = " llong_pfmt ": must have energy to do tempering!\n",
        cr->nodeid, step);

  at->H0 = eterm[F_EPOT];

  /* special bond */
  if (PAR(cr)) {
    /* communicate special energy */
    gmx_sum(AT_VTOT, at->vacc, cr);
  }

  if (at->spbs->cnt > 0) {
    at->H1 = at->vacc[AT_VSPB];
  } else if (at->bbs->cnt > 0) {
    at->H1 = at->vacc[AT_VSPB];
  } else {
    at->H1 = 0.0;
  }
  at->Ea = at->H0 + at->H1 * at->th_fd;

  if (verbose) {
    fprintf(stderr, "Sum: Ea=%g, H0=%g, H1=%g, th_f=%g, th_fd=%g\n",
        at->Ea, at->H0, at->H1, at->th_f, at->th_fd);
  }

  once++;
#else
  double epot = eterm[F_EPOT];

  if (dirty) { /* summarize local potential energy */
    die_if (!PAR(cr),
      "node: %d, step: " llong_pfmt ": no energy available\n",
        cr->nodeid, step);
    gmx_sumd(1, &epot, cr);
  }
  if (ismaster)
    at->Ea = epot;
#endif
}

/* update thermostat temperature */
static void atgmx_updTstat(at_t *at, t_inputrec *ir)
{
  int i;
  for (i = 0; i < ir->opts.ngtc; i++)
    ir->opts.ref_t[i] = (real) at->T0;
}

/* initialize aT, every node calls this */
at_t *atgmx_init(const char *atcfg, unsigned fromcpt,
    gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr, int premode)
{
  at_t *at = NULL;
  int err = 0;

#if AT_VER < 2
  die_if (mtop == NULL, "no topology\n");
#endif

  if (!(cr->duty & DUTY_PP)) /* return NULL for a PME-only node */
    return NULL;

  if (SIMMASTER(cr)) {
    /* call gromacs independent routine to initalize */
    if ((at = at_open(atcfg, fromcpt, ir->delta_t, premode)) == NULL) {
      fprintf(stderr, "error occured during initialization\n");
      err = -1; goto PAR1;
    }

    atgmx_updTstat(at, ir);
  } else { /* for non-master, simply allocate space */
    die_if ( (at = calloc(1, sizeof(*at)) ) == NULL,
      "no memory for at_t.\n");
  }

PAR1: /* check error in the master branch */
  if (PAR(cr)) gmx_bcast(sizeof(err), &err, cr);
  if (err) return NULL;

#ifdef GMX_MPI
  /* tell everyone settings on the master
   * valid only for PP only node, maybe we need to
   * consider using mpi_comm_mysim for more advanced versions
   * we pass MPI_COMM_NULL to avoid the case of one-node-mpi */
  at_initmpi(at, PAR(cr) ? cr->mpi_comm_mygroup : MPI_COMM_NULL);
#endif
  atgmx_updTstat(at, ir);

#if AT_VER == 2
  /* initialize special dihedrals, needs to called even if at->spbs->cnt = 0
   * on the master, for other nodes are uninitialized */
  if (0 != atgmx_initspb(at, mtop, cr)) {
    fprintf(stderr, "%2d: spb failed to initialize\n", cr->nodeid);
    return NULL;
  }

  if (0 != atgmx_initbb(at, mtop, cr)) {
    fprintf(stderr, "%2d: bb failed to initialize\n", cr->nodeid);
    return NULL;
  }
  if (PAR(cr)) gmx_barrier(cr);

#endif
  /* set initial scales (esp. useful for a prerun) */
  atgmx_updscl(at, cr);

  return at;
}

/* only a PP-processor calls this function
 * also assume global_stat() has been called
 * SIMMASTER(cr) should be equivalent to MASTER(cr) */
int atgmx_move(at_t *at, gmx_enerdata_t *enerd,
    llong_t step, int bFirstStep, int bLastStep,
    int bGStat, int bXTC, int bNS,
    t_commrec *cr)
{
  int do_temp, dirty;
  int ismaster = SIMMASTER(cr);

  /* nsttemp < 0 means do tempering at an NS step */
  do_temp = (at->nsttemp > 0) || bNS || bLastStep;
  if (at->nsttemp > 0 && (step % at->nsttemp) != 0 && !bLastStep)
    do_temp = 0; /* if nsttemp is set, do tempering at a regular interval */
  if (!do_temp) return 0;

  /* no tempering during prerun, temperature is fixed */
  if (at->premode == 0) {
    dirty = PAR(cr) && !bGStat;
    /* calculate H0, H1 and Ea */
    atgmx_sum(at, enerd->term, cr, step, dirty);
    /* change temperature, and regularly write output files */
    if (ismaster) {
      die_if (0 != at_move(at, step, bFirstStep, bLastStep, bXTC),
        "#%d, step = " llong_pfmt ", error during moving master\n", cr->nodeid, step);
    }
    atgmx_updscl(at, cr); /* change scale */
  }

#if AT_VER == 2
  /* synchronize local data (hist_l, sbf_l, etc., not cache,which is updated each step)
   * we then calculate the new potential of mean force */
  /* since only a node does PP calls this, at->spbs should not be NULL */
  if (at->spbs->cnt > 0 && at->spb_collect
     && doevery(step, at->spb_nstdata, bFirstStep, bLastStep) ) {
    die_if (0 != spbs_syncdata(at->spbs),
        "spbs_syncdata failed on %d\n", cr->nodeid);
    if (ismaster && doevery(step, at->mb->av_nstsave, bFirstStep, bLastStep)) {
      spbs_writebin(at->spbs, at->spbs->bin_file, /* ver = */ 2);
      spbs_write(at->spbs, at->spbs->txt_file, /* ver */ 2, /* write as is */ 1);
    }
  }

  if (at->spbs->cnt <= 0 && at->bbs->cnt > 0 && at->bb_collect) {
    if ( doevery(step, at->bb_nstflush, bFirstStep, bLastStep) ) {
      die_if (0 != bbs_syncdata(at->bbs),
          "bbs_syncdata failed on %d\n", cr->nodeid);
    }
    if (ismaster && doevery(step, at->mb->av_nstsave, bFirstStep, bLastStep)) {
      bbs_writebin(at->bbs, at->bbs->bin_file, /* ver = */ 2);
      bbs_write(at->bbs, at->bbs->txt_file, /* ver = */ 0);
    }
  }
#endif
  return 0;
}

/* to be used as a replacement of opt2fn(),
 * it will replace the file extension from .mdp to .cfg */
char *atgmx_opt2fn(const char *opt, int nfile, const t_filenm fnm[])
{
  int i;
  for (i = 0; i < nfile; i++) {
    if (strcmp(opt, fnm[i].opt) == 0) { /* a match is found */
      char *fname = fnm[i].fns[0], *p;
      if (fnm[i].ftp == efMDP) {
        /* modify the extension from .mdp to .cfg */
        if (strcmp(fname, "grompp.mdp") == 0) { /* replace the default name by NULL */
          return NULL; /* we do not have default name for .cfg files */
        } else if ((p = strstr(fname, ".cfg.mdp")) != NULL) {
          p[4] = '\0';
        } else if ((p = strstr(fname, ".mdp")) != NULL) {
          strcpy(p, ".cfg");
        }
      }
      return fname;
    }
  }
  return NULL;
}

#if AT_VER == 2
#include "txtdump.h"
#include "bondf.h"
#include "qmmm.h"

/* compute difference of two points, with periodic boundary condition,
 * cf. gmxlib/bondfree.c/pbc_rvec_sub()
 * this is usually reduced to dx = xi-xj, only for bMolPBC it calls pbc_dx_aiuc,
 * defined in gmxlib/pbc.c/pbc_dx_aiuc
 * NOTE: different from GROMACS, the first parameter is dx, instead of pbc */
int atgmx_pbcdiff(real *dx, const real *xi, const real *xj, const void *pbc)
{
  if (pbc) {
    return pbc_dx_aiuc((const t_pbc*) pbc, xi, xj, dx);
  } else {
    rv3_diff(dx, xi, xj);
    return CENTRAL;  /* = (5*3*3)/2 = 22  pbc box id */
  }
}

/* add force from biased dihedral potential
 * cf. gmxlib/bondfree.c/do_dih_fup() */
static int atgmx_updfdih(dihcalc_t *dih, int id[],
    double mfph, double lambda,
    rvec f[], rvec fshift[], const t_pbc *pbc, const t_graph *g,
    const rvec x[])
{
  real sc = (real)(mfph*lambda);
  rvec f_i, f_j, f_k, f_l;
  int i, j, k, l;

  i = id[0];
  j = id[1];
  k = id[2];
  l = id[3];
  rv3_smul2(f_i, dih->g[0], sc);
  rv3_smul2(f_j, dih->g[1], sc);
  rv3_smul2(f_k, dih->g[2], sc);
  rv3_smul2(f_l, dih->g[3], sc);
  rv3_inc(f[i], f_i);
  rv3_inc(f[j], f_j); /* note, GROMACS uses rvec_dec(), b/c its f_j is -f_j */
  rv3_inc(f[k], f_k);
  rv3_inc(f[l], f_l);

  if (g) {
    ivec jt, dt_ij, dt_kj, dt_lj;
    copy_ivec(SHIFT_IVEC(g, j), jt);
    ivec_sub (SHIFT_IVEC(g, i), jt, dt_ij);
    ivec_sub (SHIFT_IVEC(g, k), jt, dt_kj);
    ivec_sub (SHIFT_IVEC(g, l), jt, dt_lj);
    dih->t1 = IVEC2IS(dt_ij);
    dih->t2 = IVEC2IS(dt_kj);
    dih->t3 = IVEC2IS(dt_lj);
  } else if (pbc) {
    rvec dx_lj;
    dih->t3 = atgmx_pbcdiff(dx_lj, x[l], x[j], pbc);
  } else {
    dih->t3 = CENTRAL;
  }
  rv3_inc(fshift[dih->t1], f_i);
  rv3_inc(fshift[CENTRAL], f_j); /* see prev. for rvec_dec */
  rv3_inc(fshift[dih->t2], f_k);
  rv3_inc(fshift[dih->t3], f_l);
  return 0;
}

/* specially modified version for R.B. dihedrals */
real atgmx_rbdihs(int nbonds,
    const t_iatom forceatoms[], const t_iparams forceparams[],
    const rvec x[], rvec f[], rvec fshift[],
    const t_pbc *pbc, const t_graph *g,
    const t_mdatoms *md, t_fcdata *fcd,
    int *gindex,
    at_t *at, llong_t step, double *vsp)
{
  const real c0 = 0.0f, c1 = 1.0f, c2 = 2.0f, c3 = 3.0f, c4 = 4.0f, c5 = 5.0f;
  int  type, ai, aj, ak, al, i, j;
  int  t1, t2, t3;
  rvec r_ij, r_kj, r_kl, m, n;
  real parm[NR_RBDIHS];
  real phi, cos_phi, rbp;
  real v, sgn, ddphi, sin_phi;
  real cosfac, vtot;

  spbonds_t *spbs = at->spbs;
  spb_t *spb;
  int ispb, iex, isspec;
  int idx[4]; /* local indices of the four involved atoms */
  int gidx[4]; /* global indices */

  (void) md; (void) fcd;

  vtot = 0.0;
  *vsp = 0.0;
  for (i = 0; (i < nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];

    /* We always first apply the default GROMACS code first */
#if defined(GMXVERSION) && (GMXVERSION >= 40099)
    phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                  &sgn, &t1, &t2, &t3); /*  84 */
    cos_phi = (real) cos(phi);
#else
    phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                  &cos_phi, &sgn, &t1, &t2, &t3); /*  84  */
#endif

    /* Change to polymer convention */
    if (phi < c0)
      phi += M_PI;
    else
      phi -= M_PI;       /* 1 */
    cos_phi = -cos_phi;  /* 1 */

    sin_phi = (real) sin(phi);

    for (j = 0; (j < NR_RBDIHS); j++) {
      parm[j] = forceparams[type].rbdihs.rbcA[j];
    }
    /* Calculate cosine powers */
    /* Calculate the energy */
    /* Calculate the derivative */
    v       = parm[0];
    ddphi   = c0;
    cosfac  = c1;
    rbp     = parm[1];
    ddphi  += rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[2];
    ddphi  += c2*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[3];
    ddphi  += c3*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[4];
    ddphi  += c4*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[5];
    ddphi  += c5*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;

    ddphi = -ddphi*sin_phi;                     /*  11      */

    do_dih_fup(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
               f, fshift, pbc, g, x, t1, t2, t3);           /* 112      */
    vtot += v;

    /* we first use the global interaction index `type'
     * to determine if it is a special one */
    isspec = 0;
    for (ispb = 0; ispb < spbs->cnt; ispb++) {
      if (spbs->arr[ispb].type == type)
        break;
    }
    if (ispb < spbs->cnt) { /* found a type match */
      spb = spbs->arr + ispb;

      /* set local indices */
      idx[0] = ai;
      idx[1] = aj;
      idx[2] = ak;
      idx[3] = al;

      /* map the local indices to the global ones */
      for (j = 0; j < 4; j++)
        gidx[j] = (gindex != NULL) ? gindex[ idx[j] ] : idx[j];

      if (spb->exclcnt <= 0) {
        isspec = 1; /* a shortcut in case we have no exclusion list */
      } else {
        /* search the exclusion list */
        for (iex = 0; iex < spb->exclcnt; iex++) {
          for (j = 0; j < 4; j++)
            if (gidx[j] != spb->excl[iex * 4 + j])
              break;
          if (j >= 4) break; /* a match in the exclusion list is found */
        }
        isspec = (iex >= spb->exclcnt); /* no match in exclusion list */
      }
    }

    if (isspec) { /* use our own code, we assume spb and idx are set correctly */
      dihcalc_t dih;
      unsigned int dcflags = (at->spb_dihends ? DIH_ENDS : DIH_ALL);
      double mfph = 0., dblphi;

      memset(&dih, '\0', sizeof(dih));
      dih.szreal = sizeof(real);
      dih.pbcdiff = &atgmx_pbcdiff;
      dih.pbcdata = (const void *) pbc;
      dblphi = rv3_calcdihv(&dih, x, idx, dcflags);
      if (dblphi < -M_PI) dblphi = -M_PI;
      else if (dblphi > (M_PI - 1e-14)) dblphi = M_PI - 1e-14;
      dblphi += (dblphi < 0.0) ? M_PI : (-M_PI);
      if (at->spb_biased) {
        v = (real) spb_getpmf(spb, dblphi, &mfph);
      } else {
        v = 0.f;
        mfph = 0.;
      }

      /* buffer the dihedral information, wait for the force to use spb_add */
      if (at->spb_collect && step > at->nsteql
          && doevery(step, at->spb_nstdata, 0, 0)) {
        spb_buf(spb, idx, dih.g, dblphi, dih.div, dih.g2, mfph);
      }

      /* add a force component to the force profile */
      if (at->spb_biased) {
        atgmx_updfdih(&dih, idx,
          mfph, at->scale[AT_LAMBDA],
          f, fshift, pbc, g, x);
      }
      *vsp += v; /* this part of the energy is not included in vtot */
    }
  }
  return vtot;
}

#if GMXVERSION >= 40099
/* note we totally override the original cmap, if any */
real atgmx_cmap(int nbonds,
               const t_iatom forceatoms[],const t_iparams forceparams[],
               const gmx_cmap_t *cmap_grid,
               const rvec x[],rvec f[],rvec fshift[],
               const t_pbc *pbc,const t_graph *g,
               const t_mdatoms *md,t_fcdata *fcd,
               int *gindex,
               at_t *at, llong_t step, double *vsp)
{
  int n, ai, aj, ak, al, ah, type;
  int i, j, isspec, ibb, iex, idx[BB_ATOMS], gidx[BB_ATOMS];
  bb_t *bb;
  bbs_t *bbs = at->bbs;

  *vsp = 0.0;
  for (n = 0; n < nbonds; )
  {
    /* Five atoms are involved in the two torsions */
    type   = forceatoms[n++];
    /* assign local indices */
    for (i = 0; i < BB_ATOMS; i++)
      idx[i] = forceatoms[n++];

    /* we first use the global interaction index `type'
     * to determine if it is a special one */
    isspec = 0;
    for (ibb = 0; ibb < bbs->cnt; ibb++) {
      if (bbs->arr[ibb].type == type)
        break;
    }
    if (ibb < bbs->cnt) { /* found a type match */
      bb = bbs->arr + ibb;

      /* map the local indices to the global ones */
      for (j = 0; j < BB_ATOMS; j++)
        gidx[j] = (gindex != NULL) ? gindex[ idx[j] ] : idx[j];

      if (bb->exclcnt <= 0) {
        isspec = 1; /* a shortcut in case we have no exclusion list */
      } else {
        /* search the exclusion list */
        for (iex = 0; iex < bb->exclcnt; iex++) {
          for (j = 0; j < BB_ATOMS; j++)
            if (gidx[j] != bb->excl[iex * BB_ATOMS + j])
              break;
          if (j >= BB_ATOMS) break; /* a match in the exclusion list is found */
        }
        isspec = (iex >= bb->exclcnt); /* no match in exclusion list */
      }
    }

    if (isspec) {
      dihcalc_t dcphi, dcpsi;
      unsigned int dcflags = DIH_ENDS;
      double phi, psi, mfphi = 0., mfpsi = 0., v;

      memset(&dcphi, '\0', sizeof(dcphi));
      dcphi.szreal = sizeof(real);
      dcphi.pbcdiff = &atgmx_pbcdiff;
      dcphi.pbcdata = (const void *) pbc;
      phi = rv3_calcdihv(&dcphi, x, idx, dcflags);
      phi += (phi < 0.) ? M_PI : (-M_PI);

      memset(&dcpsi, '\0', sizeof(dcpsi));
      dcpsi.szreal = sizeof(real);
      dcpsi.pbcdiff = &atgmx_pbcdiff;
      dcpsi.pbcdata = (const void *) pbc;
      psi = rv3_calcdihv(&dcpsi, x, idx+1, dcflags);
      psi += (psi < 0.) ? M_PI : (-M_PI);

      if (at->bb_biased) {
        v = bb_getpmf(bb, phi, &mfphi, psi, &mfpsi);
      } else {
        v = 0.;
        mfphi = mfpsi = 0.;
      }

      /* buffer the dihedral information, wait for the force to use bb_add */
      if (at->bb_collect && step > at->nsteql
          && doevery(step, at->bb_nstdata, 0, 0)) {
        bb_buf(bb, idx,
            phi, mfphi, dcphi.g2, dcphi.g[0], dcphi.g[3], /* only 1-4 atoms */
            psi, mfpsi, dcpsi.g2, dcpsi.g[0], dcpsi.g[3]);
      }

      /* add a force component to the force profile */
      if (at->bb_biased) {
        atgmx_updfdih(&dcphi, idx,
          mfphi, at->scale[AT_LAMBDA],
          f, fshift, pbc, g, x);
        atgmx_updfdih(&dcpsi, idx + 1,
          mfpsi, at->scale[AT_LAMBDA],
          f, fshift, pbc, g, x);
      }
      *vsp += v; /* this part of the energy is not included in vtot */
    }
  }
  return 0.0f;
}
#endif

#include "nonbonded.h"

/* bonded interaction
 * corresponds to gmxlib/bondfree.c/calc_bonds(), update to 4.5.2 */
static void atgmx_calcbonds(t_commrec *cr, FILE *fplog, const gmx_multisim_t *ms,
  const t_idef *idef,
  rvec x[], history_t *hist,
  rvec f[], t_forcerec *fr,
  const t_pbc *pbc, const t_graph *g,
  gmx_enerdata_t *enerd, t_nrnb *nrnb,
  real lambda,
  const t_mdatoms *md,
  t_fcdata *fcd, int *global_atom_index,
  int bPrintSepPot, llong_t step, at_t*at)
{
  int    ftype, nbonds, ind, nat1;
  real   *epot, v, dvdl;
  const  t_pbc *pbc_null;
  static int once[F_NRE], calls, spbfreq = 500000;

  /* only for inifinite nanotube, NULL for us */
  pbc_null = (fr->bMolPBC) ? pbc : NULL;

#ifdef DEBUG
  if (g && debug)
    p_graph(debug,"Bondage is fun",g);
#endif

  epot = enerd->term;

  if (!calls) {
    const char *nfreq = getenv("SPBFREQ");
    if (nfreq != NULL) spbfreq = atoi(nfreq);
  }

  /* Do pre force calculation stuff which might require communication */
  if (idef->il[F_ORIRES].nr) {
    epot[F_ORIRESDEV] = calc_orires_dev(ms, idef->il[F_ORIRES].nr,
                                        idef->il[F_ORIRES].iatoms,
                                        idef->iparams, md, (const rvec*) x,
                                        pbc_null, fcd, hist);
  }
  if (idef->il[F_DISRES].nr) {
    calc_disres_R_6(ms, idef->il[F_DISRES].nr,
                    idef->il[F_DISRES].iatoms,
                    idef->iparams, (const rvec*) x, pbc_null,
                    fcd, hist);
  }

  /* Loop over all bonded force types to calculate the bonded forces */
  for (ftype = 0; (ftype < F_NRE); ftype++) {
#if GMXVERSION >= 40099
    if ( !(ftype < F_GB12 || ftype > F_GB14) ) continue;
#endif
    if (interaction_function[ftype].flags & IF_BOND &&
        !(ftype == F_CONNBONDS || ftype == F_POSRES)) {
      nbonds = idef->il[ftype].nr;
      if (nbonds > 0) {
        ind = interaction_function[ftype].nrnb_ind;
        nat1 = interaction_function[ftype].nratoms + 1;
        dvdl = 0;
        if (ftype < F_LJ14 || ftype > F_LJC_PAIRS_NB) {
          /* determine if a modified function is used or not */
          int accspb = (ftype == F_RBDIHS && at->spbs->cnt > 0);
          if (accspb) {
            double vsp = 0.0;
            v = atgmx_rbdihs(nbonds, idef->il[ftype].iatoms,
                idef->iparams,
                (const rvec*) x, f, fr->fshift,
                pbc_null, g, md, fcd,
                global_atom_index,
                at, step, &vsp);
            at->vacc[AT_VSPB] += (real) vsp;

            if (spbfreq > 0 && calls % spbfreq == 0)
              fprintf(stderr, "#%3d %.1f/%.1fM %s(%d) b=%.4f l=%.4f v=%+8.2f,%+8.2f upd:%d\n",
                cr->nodeid, calls*1e-6, step*1e-6, interaction_function[ftype].longname,
                ftype, at->beta, at->scale[AT_LAMBDA], v, vsp, at->spb_collect);
          }
#if GMXVERSION >= 40099
          else if (ftype == F_CMAP && at->bbs->cnt > 0) {
            double vsp = 0.0;
            v = atgmx_cmap(nbonds, idef->il[ftype].iatoms,
                idef->iparams, &idef->cmap_grid,
                (const rvec *) x, f, fr->fshift,
                pbc_null, g, md, fcd,
                global_atom_index,
                at, step, &vsp);
            at->vacc[AT_VSPB] += (real) vsp;

            if (spbfreq > 0 && calls % spbfreq == 0)
              fprintf(stderr, "#%3d %.1f/%.1fM %s(%d) b=%.4f l=%.4f v=%+8.2f,%+8.2f upd:%d\n",
                cr->nodeid, calls*1e-6, step*1e-6, interaction_function[ftype].longname,
                ftype, at->beta, at->scale[AT_LAMBDA], v, vsp, at->bb_collect);
          }
          else if (ftype == F_CMAP)
          {
            v = cmap_dihs(nbonds,idef->il[ftype].iatoms,
                          idef->iparams,&idef->cmap_grid,
                          (const rvec*)x,f,fr->fshift,
                          pbc_null,g,lambda,&dvdl,md,fcd,
                          global_atom_index);
          }
#endif
          else
          {
            v =
            interaction_function[ftype].ifunc(nbonds, idef->il[ftype].iatoms,
                                              idef->iparams,
                                              (const rvec*) x, f, fr->fshift,
                                              pbc_null, g, lambda, &dvdl, md, fcd,
                                              global_atom_index);
          }

          if (bPrintSepPot) {
            fprintf(fplog,"  %-23s #%4d  V %12.5e  dVdl %12.5e\n",
                    interaction_function[ftype].longname, nbonds/nat1, v, dvdl);
          }
          once[ftype]++;
        } else {
          v = do_listed_vdw_q(ftype, nbonds, idef->il[ftype].iatoms,
                              idef->iparams,
                              (const rvec*) x, f, fr->fshift,
                              pbc_null, g,
                              lambda, &dvdl,
                              md, fr, &enerd->grpp, global_atom_index);

          if (bPrintSepPot) {
            fprintf(fplog,"  %-5s + %-15s #%4d                  dVdl %12.5e\n",
                    interaction_function[ftype].longname,
                    interaction_function[F_COUL14].longname, nbonds/nat1, dvdl);
          }
        }
        if (ind != -1)
          inc_nrnb(nrnb, ind, nbonds/nat1);
        epot[ftype]  += v;
#if GMXVERSION >= 40099
        enerd->dvdl_nonlin += dvdl;
#else
        epot[F_DVDL] += dvdl;
#endif
      }
    }
  }
  /* Copy the sum of violations for the distance restraints from fcd */
  if (fcd)
    epot[F_DISRESVIOL] = fcd->disres.sumviol;
  calls++;
}

/* compute local force for a PP node
 * correpsonds to mdlib/force.c/do_force_lowlevel()
 * 1. Wall force
 * 2. Bonded
 * 3. nonbonded
 * 4. electrostatic correction (PP-PME node)
 * 5. PME (PP-PME node)
 *
 * + In 4. and 5. forces are saved to fr->f_novirsum
 * + 4. and 5. can be replaced by their RF counterparts
 * + We assume flags=GMX_FORCE_ALLFORCES,
 *
 * scale[] is equal to at->scale (no offset)
 */
void atgmx_forcelow(FILE *fplog, llong_t step, t_forcerec *fr, t_inputrec *ir, t_idef *idef,
  t_commrec *cr, t_nrnb *nrnb, gmx_wallcycle_t wcycle, t_mdatoms *md, t_grpopts  *opts,
  rvec x[], history_t  *hist, rvec f[], gmx_enerdata_t *enerd, t_fcdata *fcd, matrix box,
  real lambda, t_graph *graph, t_blocka *excl, rvec mu_tot[],
  int flags,
  float *cycles_force, at_t *at)
{
    int     i, status;
    int    bSepDVDL, bSB;
    matrix  boxs;
    rvec    box_size;
    t_pbc   pbc;
    real    dvdlambda, Vsr, Vlr, Vcorr[3] = { 0 };
    real    vdip, vcharge;
    real   **earr = enerd->grpp.ener;

#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,sepdvdlformat,s,v,dvdl);

    (void) opts;
    GMX_MPE_LOG(ev_force_start);

    /* Reset box */
    for (i = 0; (i < DIM); i++)
    {
        box_size[i] = box[i][i];
    }

    bSepDVDL = (fr->bSepDVDL && do_per_step(step, ir->nstlog));

    /* do QMMM first if requested */
    if (fr->bQMMM)
    {
        enerd->term[F_EQM] = calculate_QMMM(cr, x, f, fr, md);
    }

    /* Call the short range functions all in one go. */
    GMX_MPE_LOG(ev_do_fnbf_start);

    dvdlambda = 0;


    /* 1. Wall force */
    if (ir->nwall)
    {
        dvdlambda = do_walls(ir, fr, box, md, x, f, lambda,
                             earr[egLJSR], nrnb);
        PRINT_SEPDVDL("Walls", 0.0, enerd->term[F_DVDL]);
        enerd->term[F_DVDL] += dvdlambda;
    }

    /* 2. Nonbonded force */
    do_nonbonded(cr, fr, x, f, md, fr->bBHAM ? earr[egBHAMSR] : earr[egLJSR],
        earr[egCOULSR], box_size, nrnb, lambda, &dvdlambda,
        FALSE/*short range*/, -1, -1, flags & GMX_FORCE_FORCES);

    enerd->term[F_DVDL] += dvdlambda;
    Vsr = 0;
    for (i = 0; i < enerd->grpp.nener; i++)
      Vsr += (fr->bBHAM ? earr[egBHAMSR][i] : earr[egLJSR][i])
           + earr[egCOULSR][i];
    enerd->term[F_DVDL] += Vsr;
    PRINT_SEPDVDL("VdW and Coulomb SR particle-p.", Vsr, enerd->term[F_DVDL]);

    GMX_MPE_LOG(ev_do_fnbf_finish);

    if (debug)
    {
        pr_rvecs(debug, 0,"fshift after SR", fr->fshift, SHIFTS);
    }

    /* Shift the coordinates. Must be done before bonded forces and PPPM,
     * but is also necessary for SHAKE and update, therefore it can NOT
     * go when no bonded forces have to be evaluated.
     */

    /* Here sometimes we would not need to shift with NBFonly,
     * but we do so anyhow for consistency of the returned coordinates.
     */
    if (graph)
    {
        shift_self(graph, box, x);
        if (TRICLINIC(box))
        {
            inc_nrnb(nrnb, eNR_SHIFTX, 2*graph->nnodes);
        }
        else
        {
            inc_nrnb(nrnb, eNR_SHIFTX, graph->nnodes);
        }
    }
    /* Check whether we need to do bondeds or correct for exclusions */
    if (fr->bMolPBC)
    {
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc, fr->ePBC, cr->dd, TRUE, box);
    }

    /* 3. bonds, LJ14, etc */
    GMX_MPE_LOG(ev_calc_bonds_start);
    atgmx_calcbonds(cr, fplog, cr->ms, idef, x, hist, f, fr, &pbc, graph, enerd, nrnb, lambda, md, fcd,
          DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL, fr->bSepDVDL && do_per_step(step, ir->nstlog), step,
          at);
    GMX_MPE_LOG(ev_calc_bonds_finish);

    if (EEL_FULL(fr->eeltype))
    {
        bSB = (ir->nwall == 2);
        if (bSB)
        {
            copy_mat(box, boxs);
            svmul(ir->wall_ewald_zfac, boxs[ZZ], boxs[ZZ]);
            box_size[ZZ] *= ir->wall_ewald_zfac;
        }

        clear_mat(fr->vir_el_recip);

        if (fr->bEwald)
        {
            dvdlambda = 0;
            Vcorr[0] = ewald_LRcorrection(fplog, md->start, md->start + md->homenr, cr, fr,
              md->chargeA, md->nChargePerturbed ? md->chargeB : NULL, excl, x, bSB ? boxs : box,
              mu_tot, ir->ewald_geometry, ir->epsilon_surface, lambda, &dvdlambda, &vdip, &vcharge);
            PRINT_SEPDVDL("Ewald excl./charge/dip. corr.", Vcorr[0], dvdlambda);
            enerd->term[F_DVDL] += dvdlambda;
        }
        else
        {
            Vcorr[0] = shift_LRcorrection(fplog, md->start, md->homenr, cr, fr, md->chargeA, excl, x, TRUE, box, fr->vir_el_recip);
        }

        *cycles_force = (float) wallcycle_stop(wcycle, ewcFORCE);
        /* Now we can do communication again */

        dvdlambda = 0;
        status = 0;
        switch (fr->eeltype)
        {
        case eelPPPM:
            status = gmx_pppm_do(fplog, fr->pmedata, FALSE, x, fr->f_novirsum, md->chargeA,
              box_size, fr->phi, cr, md->start, md->homenr, nrnb, ir->pme_order, &Vlr);
            break;
        case eelPME:
        case eelPMESWITCH:
        case eelPMEUSER:
            if (cr->duty & DUTY_PME)
            {
                wallcycle_start(wcycle, ewcPMEMESH);
                /* energy is saved to Vlr */
                status = gmx_pme_do(fr->pmedata, md->start, md->homenr, x, fr->f_novirsum, md->chargeA, md->chargeB,
                  bSB ? boxs : box, cr, DOMAINDECOMP(cr) ? dd_pme_maxshift(cr->dd) : 0,
                  nrnb, fr->vir_el_recip, fr->ewaldcoeff, &Vlr, lambda, &dvdlambda, FALSE);
                PRINT_SEPDVDL("PME mesh", Vlr, dvdlambda);
                wallcycle_stop(wcycle, ewcPMEMESH);
            }
            else
            {
                /* Energies and virial are obtained later from the PME nodes */
                /* but values have to be zeroed out here */
                Vlr = 0.0;
            }
            break;
        case eelEWALD:
            Vlr = do_ewald(fplog, FALSE, ir, x, fr->f_novirsum, md->chargeA, md->chargeB,
              box_size, cr, md->homenr, fr->vir_el_recip, fr->ewaldcoeff, lambda, &dvdlambda);
            PRINT_SEPDVDL("Ewald long-range", Vlr, dvdlambda);
            break;
        default:
            Vlr = 0;
            gmx_fatal(FARGS,"No such electrostatics method implemented %s",
                      eel_names[fr->eeltype]);
        }
        if (status != 0)
        {
            gmx_fatal(FARGS,"Error %d in long range electrostatics routine %s",
                      status, EELTYPE(fr->eeltype));
        }
        enerd->term[F_DVDL] += dvdlambda;
        enerd->term[F_COUL_RECIP] = Vlr + Vcorr[0];
    }
    else
    {
        if (EEL_RF(fr->eeltype))
        {
            dvdlambda = 0;

            if (fr->eeltype != eelRF_NEC)
            {
              enerd->term[F_RF_EXCL] = RF_excl_correction(fplog, fr, graph, md, excl, x, f, fr->fshift, &pbc, lambda, &dvdlambda);
            }
            enerd->term[F_DVDL] += dvdlambda;
            PRINT_SEPDVDL("RF exclusion correction", enerd->term[F_RF_EXCL], dvdlambda);
        }
        *cycles_force = (float) wallcycle_stop(wcycle, ewcFORCE);
        /* Now we can do communication again */
    }

    GMX_MPE_LOG(ev_force_finish);
}

static real sumarr(int n, real v[])
{
  real t;
  int i;
  for (t = 0.0, i = 0; i < n; i++) t += v[i];
  return t;
}
#endif

#ifdef PRINT_SEPDVDL
#undef PRINT_SEPDVDL
#endif
#define PRINT_SEPDVDL(yn,msg,V,dVdl) if(fplog && yn) fprintf(fplog,"%s: dVdl=%g(this) %g(total)\n",msg,dVdl,V);

/* calculate global force
 * the energies, however, are still local
 * always assumes bDoForces = TRUE,bStateChanged=TRUE,
 * as it is called in md(), bFillGrid==bNS
 * */
void atgmx_doforce(FILE *fplog, t_commrec *cr,
    t_inputrec *inputrec,
    llong_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
    gmx_localtop_t *top,
    gmx_groups_t *groups,
    matrix box, rvec x[], history_t *hist,
    rvec f[], rvec buf[],
    tensor vir_force,
    t_mdatoms *mdatoms,
    gmx_enerdata_t *enerd, t_fcdata *fcd,
    real lambda, t_graph *graph,
    t_forcerec *fr, gmx_vsite_t *vsite, rvec mu_tot,
    real t, FILE *field, gmx_edsam_t ed,
    int flags, at_t *at)
{
#if AT_VER == 2
  int    cg0, cg1, i, j;
  int    start, homenr;
  static double mu[2*DIM];
  rvec   mu_tot_AB[2];
  int    bSepDVDL, bStateChanged, bNS, bFillGrid, bCalcCGCM, bBS = 0, bDoForces;
  matrix boxs = {{ 0 }};
  real   e = 0.0f, v, dvdl;
  t_pbc  pbc;
  float  cycles_ppdpme, cycles_pme = 0.0f, cycles_force;
  int istart, iend;
  real *eterm = enerd->term, **earr = enerd->grpp.ener;

  (void) bBS; (void) boxs; (void) cycles_pme; (void) e; /* appease the single-processor case */
  start  = mdatoms->start;
  homenr = mdatoms->homenr;

  bSepDVDL = (fr->bSepDVDL && do_per_step(step, inputrec->nstlog));

  clear_mat(vir_force);

  if (PARTDECOMP(cr)) {
    pd_cg_range(cr, &cg0, &cg1);
  } else {
    cg0 = 0;
    if (DOMAINDECOMP(cr))
      cg1 = cr->dd->ncg_tot;
    else
      cg1 = top->cgs.nr;
    if (fr->n_tpi > 0)
      cg1--;
  }

  bStateChanged = (flags & GMX_FORCE_STATECHANGED);
  bNS           = (flags & GMX_FORCE_NS);
  bFillGrid     = (bNS && bStateChanged);
  bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
  bDoForces     = (flags & GMX_FORCE_FORCES);

  if (bStateChanged) {
    update_forcerec(fplog, fr, box);

    /* Calculate total (local) dipole moment in a temporary common array.
     * This makes it possible to sum them over nodes faster.
     */
    calc_mu(start, homenr,
        x, mdatoms->chargeA, mdatoms->chargeB, mdatoms->nChargePerturbed,
        mu, mu + DIM);
  }

  if (fr->ePBC != epbcNONE) {
    /* Compute shift vectors every step,
     * because of pressure coupling or box deformation!
     */
    if (DYNAMIC_BOX(*inputrec) && bStateChanged)
      calc_shifts(box, fr->shift_vec);

    if (bCalcCGCM) {
      put_charge_groups_in_box(fplog, cg0, cg1, fr->ePBC, box,
       &(top->cgs), x, fr->cg_cm);
      inc_nrnb(nrnb, eNR_CGCM, homenr);
      inc_nrnb(nrnb, eNR_RESETX, cg1 - cg0);
    }
    else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph) {
      unshift_self(graph, box, x);
    }
  }
  else if (bCalcCGCM) {
    calc_cgcm(fplog, cg0, cg1, &(top->cgs), x, fr->cg_cm);
    inc_nrnb(nrnb, eNR_CGCM, homenr);
  }

  if (bCalcCGCM) {
    if (PAR(cr)) {
      move_cgcm(fplog, cr, fr->cg_cm);
    }
    if (gmx_debug_at)
      pr_rvecs(debug, 0,"cgcm", fr->cg_cm, top->cgs.nr);
  }

#ifdef GMX_MPI
  if (!(cr->duty & DUTY_PME)) {
    /* Send particle coordinates to the pme nodes.
     * Since this is only implemented for domain decomposition
     * and domain decomposition does not use the graph,
     * we do not need to worry about shifting.
     */

    wallcycle_start(wcycle, ewcPP_PMESENDX);
    GMX_MPE_LOG(ev_send_coordinates_start);

    bBS = (inputrec->nwall == 2);
    if (bBS) {
      copy_mat(box, boxs);
      svmul(inputrec->wall_ewald_zfac, boxs[ZZ], boxs[ZZ]);
    }

    gmx_pme_send_x(cr, bBS ? boxs : box, x, mdatoms->nChargePerturbed, lambda);

    GMX_MPE_LOG(ev_send_coordinates_finish);
    wallcycle_stop(wcycle, ewcPP_PMESENDX);
  }
#endif /* GMX_MPI */

  /* Communicate coordinates and sum dipole if necessary */
  if (PAR(cr)) {
    wallcycle_start(wcycle, ewcMOVEX);
    if (DOMAINDECOMP(cr)) {
      dd_move_x(cr->dd, box, x, buf);
    } else {
      move_x(fplog, cr, GMX_LEFT, GMX_RIGHT, x, nrnb);
    }
    /* When we don't need the total dipole we sum it in global_stat */
    if (NEED_MUTOT(*inputrec))
      gmx_sumd(2*DIM, mu, cr);
    wallcycle_stop(wcycle, ewcMOVEX);
  }
  for (i = 0; i < 2; i++) /* for state A and state B */
    for (j = 0; j < DIM; j++)
      mu_tot_AB[i][j] = (real) mu[i*DIM + j];
  if (fr->efep == efepNO)
    copy_rvec(mu_tot_AB[0], mu_tot);
  else
    for (j = 0; j < DIM; j++)
      mu_tot[j] = (real)( (1.0 - lambda)*mu_tot_AB[0][j] + lambda*mu_tot_AB[1][j] );

  /* Reset energies */
  /* long range terms on the master at non neighbor-search steps are kept
     since these terms are already summed at the last neighbor search step. */
  for (i = 0; i < egNR; i++)
    if (bNS || (!fr->bTwinRange) || (!MASTER(cr)) || (i != egCOULLR && i != egLJLR))
      for (j = 0; j < enerd->grpp.nener; j++) enerd->grpp.ener[i][j] = 0.0;
  if (bNS || (!fr->bTwinRange)) enerd->dvdl_lr = 0.0;
  for (i = 0; i <= F_EPOT; i++) eterm[i] = 0.0;
  eterm[F_DVDL]  = enerd->dvdl_lr;
  eterm[F_DGDL_CON] = eterm[F_DKDL] = 0.0;

  if (bNS) {
    wallcycle_start(wcycle, ewcNS);

    if (graph && bStateChanged)
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(fplog, graph, fr->ePBC, box, x);

    /* Reset long range forces if necessary */
    if (fr->bTwinRange) {
      clear_rvecs(fr->f_twin_n, fr->f_twin);
      clear_rvecs(SHIFTS, fr->fshift_twin);
    }
    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl = 0;
    {
      /* Do neighbour searching and calculate long range forces and energy for twin range.  */
      static int bFirst = 1;
      GMX_MPE_LOG(ev_ns_start);
      if (bFirst) {
        /* Allocate memory for the neighbor lists */
        init_neighbor_list(fplog, fr, mdatoms->homenr);
        bFirst = FALSE;
      }
      if (fr->bTwinRange)
        fr->nlr = 0;
      i = search_neighbours(fplog, fr, x, box, top, groups, cr, nrnb, mdatoms, lambda, &dvdl, &enerd->grpp, bNS, TRUE);
      if (debug)
        fprintf(debug,"nsearch = %d\n", i);
      GMX_MPE_LOG(ev_ns_finish);
    }

    PRINT_SEPDVDL(bSepDVDL, "LR non-bonded", 0.0, dvdl);
    enerd->dvdl_lr       = dvdl;
    eterm[F_DVDL] += dvdl;
    wallcycle_stop(wcycle, ewcNS);
  }

  if (DOMAINDECOMP(cr)) {
    if (!(cr->duty & DUTY_PME)) {
      wallcycle_start(wcycle, ewcPPDURINGPME);
      dd_force_flop_start(cr->dd, nrnb);
    }
  }
  /* Start the force cycle counter.
   * This counter is stopped in do_forcelow_level.
   * No parallel communication should occur while this counter is running,
   * since that will interfere with the dynamic load balancing.
   */
  wallcycle_start(wcycle, ewcFORCE);

  if (bDoForces) {
      /* Reset PME/Ewald forces if necessary */
    if (fr->bF_NoVirSum)
    {
      GMX_BARRIER(cr->mpi_comm_mygroup);
      if (fr->bDomDec)
        clear_rvecs(fr->f_novirsum_n, fr->f_novirsum);
      else
        clear_rvecs(homenr, fr->f_novirsum + start);
      GMX_BARRIER(cr->mpi_comm_mygroup);
    }
    /* Copy long range forces into normal buffers */
    if (fr->bTwinRange) {  /* copy force/shift force from twin zone to here */
      for (i = 0; i < fr->f_twin_n; i++)
        copy_rvec(fr->f_twin[i], f[i]);
      for (i = 0; i < SHIFTS; i++)
        copy_rvec(fr->fshift_twin[i], fr->fshift[i]);
    }
    else { /* simply clear f and shift f */
      if (DOMAINDECOMP(cr))
        clear_rvecs(cr->dd->nat_tot, f);
      else
        clear_rvecs(mdatoms->nr, f);
      clear_rvecs(SHIFTS, fr->fshift);
    }
    clear_rvec(fr->vir_diag_posres);
    GMX_BARRIER(cr->mpi_comm_mygroup);
  }
  if (inputrec->ePull == epullCONSTRAINT)
    clear_pull_forces(inputrec->pull);

  /* update QMMMrec, if necessary */
  if (fr->bQMMM)
    update_QMMMrec(cr, fr, x, mdatoms, box, top);

  if ((flags & GMX_FORCE_BONDED) && top->idef.il[F_POSRES].nr > 0) {
    /* Position restraints always require full pbc */
    set_pbc(&pbc, inputrec->ePBC, box);
    v = posres(top->idef.il[F_POSRES].nr, top->idef.il[F_POSRES].iatoms,
               top->idef.iparams_posres,
               (const rvec*) x, fr->f_novirsum, fr->vir_diag_posres,
               inputrec->ePBC == epbcNONE ? NULL : &pbc, lambda, &dvdl,
               fr->rc_scaling, fr->ePBC, fr->posres_com, fr->posres_comB);
    PRINT_SEPDVDL(bSepDVDL, interaction_function[F_POSRES].longname, v, dvdl);
    enerd->term[F_POSRES] += v;
    enerd->term[F_DVDL]   += dvdl;
    inc_nrnb(nrnb, eNR_POSRES, top->idef.il[F_POSRES].nr/2);
  }

  at_clear_energy(at);

  /* Compute the bonded and non-bonded forces */
  atgmx_forcelow(fplog, step, fr, inputrec, &(top->idef),
      cr, nrnb, wcycle, mdatoms, &(inputrec->opts),
      x, hist, f, enerd, fcd, box, lambda, graph, &(top->excls), mu_tot_AB,
      flags, &cycles_force, at);
  GMX_BARRIER(cr->mpi_comm_mygroup);

  if (ed) {
    do_flood(fplog, cr, x, f, ed, box, step);
  }

  if (DOMAINDECOMP(cr)) {
    dd_force_flop_stop(cr->dd, nrnb);
    if (wcycle)
      dd_cycles_add(cr->dd, cycles_force, ddCyclF);
  }

  if (bDoForces) {
    die_if(field != NULL, "node %d: e-field must be zero\n", cr->nodeid);
    /* Compute forces due to electric field */
    /*calc_f_el(MASTER(cr) ? field : NULL, start, homenr, mdatoms->chargeA, x, f, inputrec->ex, inputrec->et, t);*/

    /* When using PME/Ewald we compute the long range virial there.
     * otherwise we do it based on long range forces from twin range
     * cut-off based calculation (or not at all).
     */

    /* Communicate the forces */
    if (PAR(cr)) {
      wallcycle_start(wcycle, ewcMOVEF);
      if (DOMAINDECOMP(cr)) {
        dd_move_f(cr->dd, f, buf, fr->fshift);
        /* Position restraint do not introduce inter-cg forces */
        if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl)
          dd_move_f(cr->dd, fr->f_novirsum, buf, NULL);
      } else { /* adding force from the two neighboring nodes */
        move_f(fplog, cr, GMX_LEFT, GMX_RIGHT, f, buf, nrnb);
      }
      wallcycle_stop(wcycle, ewcMOVEF);
    }
  }

  if (bDoForces) {
    if (vsite) {
      wallcycle_start(wcycle, ewcVSITESPREAD);
      spread_vsite_f(fplog, vsite, x, f, fr->fshift, nrnb,
                     &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
      wallcycle_stop(wcycle, ewcVSITESPREAD);
    }

    /* Calculation of the virial must be done after vsites! */
    clear_mat(vir_force);
    /* calculate virial due to shift */
    calc_vir(fplog, SHIFTS, fr->shift_vec, fr->fshift, vir_force, inputrec->ePBC == epbcSCREW, box);
    inc_nrnb(nrnb, eNR_VIRIAL, SHIFTS);

    /* Calculate partial virial, for local atoms only, based on short range.
     * Total virial is computed in global_stat(), called from md() */
    f_calc_vir(fplog, mdatoms->start, mdatoms->start + mdatoms->homenr, x, f, vir_force, graph, box);
    inc_nrnb(nrnb, eNR_VIRIAL, homenr);

    for (i = 0; i < DIM; i++) {
      vir_force[i][i] += fr->vir_diag_posres[i]; /* Add position restraint contribution */
      vir_force[i][ZZ] += (real) fr->vir_wall_z[i]; /* Add wall contribution */
    }
    if (debug)
      pr_rvecs(debug, 0,"vir_part", vir_force, DIM);
  }

  if (inputrec->ePull == epullUMBRELLA || inputrec->ePull == epullCONST_F) {
    /* Calculate the center of mass forces, this requires communication,
     * which is why pull_potential is called close to other communication.
     * The virial contribution is calculated directly,
     * which is why we call pull_potential after calc_virial.
     */
    set_pbc(&pbc, inputrec->ePBC, box);
    dvdl = 0;
    enerd->term[F_COM_PULL] =
      pull_potential(inputrec->ePull, inputrec->pull, mdatoms, &pbc,
                     cr, t, lambda, x, f, vir_force, &dvdl);
    if (bSepDVDL)
      fprintf(fplog, sepdvdlformat,"Com pull", enerd->term[F_COM_PULL], dvdl);
    enerd->term[F_DVDL] += dvdl;
  }

  if (!(cr->duty & DUTY_PME)) {
    cycles_ppdpme = (float) wallcycle_stop(wcycle, ewcPPDURINGPME);
    dd_cycles_add(cr->dd, cycles_ppdpme, ddCyclPPduringPME);
  }

#ifdef GMX_MPI
  if (PAR(cr) && !(cr->duty & DUTY_PME)) {
    /* In case of node-splitting, the PP nodes receive the long-range
     * forces, virial and energy from the PME nodes here.
     */
    wallcycle_start(wcycle, ewcPP_PMEWAITRECVF);
    dvdl = 0;
    /* the longrange force is stored in fr->f_novirsum */
    gmx_pme_receive_f(cr, fr->f_novirsum, fr->vir_el_recip, &e, &dvdl,
                      &cycles_pme);
    PRINT_SEPDVDL(bSepDVDL, "PME mesh", e, dvdl);
    enerd->term[F_COUL_RECIP] += e;
    enerd->term[F_DVDL] += dvdl;
    if (wcycle)
      dd_cycles_add(cr->dd, cycles_pme, ddCyclPME);
    wallcycle_stop(wcycle, ewcPP_PMEWAITRECVF);
  }
#endif

  if (bDoForces && fr->bF_NoVirSum) {
    if (vsite) {
      /* Spread the mesh force on virtual sites to the other particles...
       * This is parallellized. MPI communication is performed
       * if the constructing atoms aren't local.
       */
      wallcycle_start(wcycle, ewcVSITESPREAD);
      spread_vsite_f(fplog, vsite, x, fr->f_novirsum, NULL, nrnb,
                     &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
      wallcycle_stop(wcycle, ewcVSITESPREAD);
    }

    /* Now add the reciprocal space force, this is local */
    if (fr->bDomDec) {
      istart = 0;
      iend = fr->f_novirsum_n; /* domain decomposition */
    } else {
      istart = mdatoms->start;
      iend = istart + mdatoms->homenr; /* particle decomposition */
    }
    for (i = istart; i < iend; i++) {
      /* Combines the real space force and the reciprocal space force
       * rvec_inc(f[i],fr->f_novirsum[i]); */
      f[i][0] += fr->f_novirsum[i][0];
      f[i][1] += fr->f_novirsum[i][1];
      f[i][2] += fr->f_novirsum[i][2];
    }
    if (EEL_FULL(fr->eeltype)) {
      /* Add the mesh contribution to the virial */
      m_add(vir_force, fr->vir_el_recip, vir_force);
    }
    if (debug)
      pr_rvecs(debug, 0,"vir_force", vir_force, DIM);
  }

  die_if (!bDoForces,
    "at_do_force not doing force? seriously?");
#else
  /* in version 1, we just call the default do_force */
  do_force(fplog, cr, inputrec, step, nrnb, wcycle, top,
           groups,
           box, x, hist,
           f, buf,
           vir_force, mdatoms, enerd, fcd,
           lambda, graph,
           fr, vsite, mu_tot, t, field, ed,
           flags);
#endif

  if (at != NULL) {
    int k;
#if AT_VER == 2
    real scl = (real) at->scale[AT_EBG];
    double lam;

    if (at->spbs->cnt > 0
        && at->spb_collect
        && step > at->nsteql
        && doevery(step, at->spb_nstdata, 0, 0)) {
      lam = at->spb_tscal ? at->scale[AT_LAMBDA] : 1.;
      spbs_debuf(at->spbs, f, at->beta, lam);
    }
    if (at->bbs->cnt > 0
        && at->bb_collect
        && step > at->nsteql
        && doevery(step, at->bb_nstdata, 0, 0)) {
      lam = at->bb_tscal ? at->scale[AT_LAMBDA] : 1.;
      bbs_debuf(at->bbs, f, at->beta, lam);
    }
#else
    real scl = (real) at->scale;
#endif

    /* scale the force */
    for (k = mdatoms->start; k < mdatoms->start + mdatoms->homenr; k++) {
      f[k][0] *= scl;
      f[k][1] *= scl;
      f[k][2] *= scl;
    }
  }

#if AT_VER == 2
  {
    eterm[F_BHAM]     = sumarr(enerd->grpp.nener, earr[egBHAMSR]) ;
    eterm[F_COUL_SR]  = sumarr(enerd->grpp.nener, earr[egCOULSR]) ;
    eterm[F_LJ]       = sumarr(enerd->grpp.nener, earr[egLJSR])   ;
    eterm[F_COUL_LR]  = sumarr(enerd->grpp.nener, earr[egCOULLR]) ;
    eterm[F_LJ_LR]    = sumarr(enerd->grpp.nener, earr[egLJLR])   ;
    eterm[F_BHAM_LR]  = sumarr(enerd->grpp.nener, earr[egBHAMLR]) ;
    eterm[F_LJ14]     = sumarr(enerd->grpp.nener, earr[egLJ14])   ;
    eterm[F_COUL14]   = sumarr(enerd->grpp.nener, earr[egCOUL14]) ;

    for (eterm[F_EPOT] = 0, i = 0; i < F_EPOT; i++)
      if (i != F_DISRESVIOL && i != F_ORIRESDEV && i != F_DIHRESVIOL)
        eterm[F_EPOT] += eterm[i];
  }

  if (fr->print_force > 0 && MASTER(cr)) {
    real pf2 = sqr(fr->print_force), fn2;
    for (i = mdatoms->start; i < mdatoms->start + mdatoms->homenr; i++)
      if ((fn2 = norm2(f[i])) >= pf2)
        fprintf(stderr, "step " llong_pfmt "  atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
                step, ddglatnr(cr->dd, i), x[i][XX], x[i][YY], x[i][ZZ], sqrt(fn2));
  }
#endif
}


#ifndef GMXALGE_C__
#define GMXALGE_C__

/* routine GROMACS code are to be added by
 *
 *   python gmxspdr.py -b ga -p gmxalge -i gmxalge.c -o mdalge.c
 *
 * put this file in GNUmakefile */

/* %obj.decl% */

#define HAVEREAL 1
#define ZCOM_PICK
#define ZCOM_LOG
#define ZCOM_CFG
#define ZCOM_RNG
#include "zcom.h"
#include "gmxhmc.h"

typedef struct taggmxalge_t {
  int    imode;     /* command-line input mode */
  double tmstep;    /* MD integration step */
  char   *fnrng;    /* file name of random number state; */
  int    nstlog;    /* interval of writing trace file; -1: only at ns, 0: disable; $def: -1;  */
  char   *fnlog;    /* name of trace file */
  logfile_t *log;   /* logfile */

  gmxhmc_t  *hmc;   /* Hybrid-MC */
  int seglen;       /* HMC segement length */
} gmxalge_t;


/* initialize the object on the master node */
static gmxalge_t *gmxalge_open(const char *fncfg, t_state *state,
    gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr, int imode)
{
  gmxalge_t *ga;
  cfg_t *cfg;

  xnew(ga, 1);
  ga->imode = imode;

  die_if( (cfg = cfg_open(fncfg)) == NULL,
      "cannot load cfg %s\n", fncfg);

#define CFGADD(key, fmt, var, def, desc) { \
  var = def; cfg_add(cfg, key, fmt, &(var), desc); }

  CFGADD("rng", "%s", ga->fnrng, "rng.dat", "random number file");
  CFGADD("nstlog", "%d", ga->nstlog, 1000, "NST of saving log file");
  CFGADD("log", "%s", ga->fnlog, "trace.dat", "log file");
  CFGADD("seglen", "%d", ga->seglen, 100, "HMC segment length");

  die_if (0 != cfg_match(cfg, CFG_CHECKUSE|CFG_VERBOSE),
      "bad cfg %s\n", fncfg);
  cfg_close(cfg);

  ga->hmc = gmxhmc_open(state, ir, cr);

  return ga;
}


/* initialization in mdrunner() before calling `do_md'
 * applies to both PP and PME-only nodes */
static gmxalge_t *gmxalge_init(const char *fncfg, unsigned iscntn,
    t_state *state, gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr,
    int imode)
{
  gmxalge_t *ga;

  if ( SIMMASTER(cr) ) { /* initialization on the master node */
    ga = gmxalge_open(fncfg, state, mtop, ir, cr, imode);
  } else {
    xnew(ga, 1);
  }
  gmx_bcast_sim(sizeof(*ga), ga, cr);
  return ga;
}

/* release everyting at the end of mdrunner()
 * applies to both PP and PME-only nodes */
static void gmxalge_done(gmxalge_t *ga, t_commrec *cr)
{
  if ( SIMMASTER(cr) ) {
    if (ga->hmc) gmxhmc_close(ga->hmc);
  }
  if (ga != NULL) free(ga);
}

static int gmxalge_move(gmxalge_t *ga, int *dirty, FILE *fplog,
    gmx_large_int_t step, int bFirstStep, int bLastStep,
    int bGStat, int bXTC, int bNS, gmx_enerdata_t *enerd,
    t_state *state_global, t_state *state_local, rvec **f,
    gmx_mtop_t *top_global, gmx_localtop_t *top_local,
    t_inputrec *ir, t_commrec *cr, t_mdatoms *mdatoms, t_forcerec *fr,
    gmx_vsite_t *vsite, gmx_shellfc_t shellfc,
    gmx_constr_t constr, t_nrnb *nrnb, gmx_wallcycle_t wcycle)
{
  if (step % ga->seglen != 0) return 0;

  gmxhmc_push(ga->hmc, cr, state_global, state_local);
  if (rnd0() < 1e-16) {
    gmxhmc_pop(ga->hmc, fplog, step, state_global, state_local, f,
        top_global, top_local, ir, cr, mdatoms, fr,
        vsite, shellfc, constr, nrnb, wcycle);
    if (dirty) *dirty = 1;
  } else {
    if (dirty) *dirty = 0;
  }
  return 0;
}

/* % NO bondfree.c% */
/* % NO force.c% */
/* % NO sim_util.c% */
/* %md.c% */
/* %runner.c% */
/* %mdrun.c% */


#endif /* GMXALGE_C__ */


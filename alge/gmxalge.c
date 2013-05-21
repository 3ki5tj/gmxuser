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
#define USE_MPI GMX_MPI

/* BOLTZ is defined in GROMACS */
#ifndef BOLTZ
#define BOLTZ  8.314511212e-3
#endif

typedef struct taggmxalge_t {
  int    imode;     /* command-line input mode */
  double tmstep;    /* MD integration step */
  double beta;      /* current inverse temperature */
  char   *fnrng;    /* file name of random number state; */
  int    nstlog;    /* interval of writing trace file; -1: only at ns, 0: disable; $def: -1;  */
  char   *fnlog;    /* name of trace file */
  logfile_t *log;   /* logfile */
} gmxalge_t;

static gmxalge_t *gmxalge_init(const char *fncfg, unsigned fromcpt,
    gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr, int imode)
{
  gmxalge_t *ga;

  xnew(ga, 1);
  ga->imode = imode;
  return ga;
}

static void gmxalge_done(gmxalge_t *ga)
{
  free(ga);
}

static int gmxalge_move(gmxalge_t *ga, gmx_enerdata_t *enerd,
    gmx_large_int_t step, int bFirstStep, int bLastStep,
    int bGStat, int bXTC, int bNS, t_commrec *cr)
{
  return 0;
}

/* %no bondfree.c% */
/* %no force.c% */
/* %no sim_util.c% */
/* %md.c% */
/* %runner.c% */
/* %mdrun.c% */


#endif /* GMXALGE_C__ */


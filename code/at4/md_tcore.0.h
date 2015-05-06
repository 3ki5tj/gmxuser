/* Accelerated tempering core
 * Copyright (c) 2009-2010 Cheng Zhang */
#define HAVEREAL 1
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ENDN
#define ZCOM_RNG
#define ZCOM_LOG
#if AT_VER == 2
#define ZCOM_RV3
#endif
#include "zcom.h"

#include "mb_t.h"

#ifndef BOLTZ
#define BOLTZ  8.314511212e-3 /* Boltzmann constant */
#endif

#if AT_VER == 2

#include "md2spb.h" /* independent spb */
#include "md2bb.h"  /* pair of dihedrals */

#define AT_VSPB  0                /* index for special dihedrals */
#define AT_VLAST AT_VSPB          /* index for the last additional energy item */
#define AT_VTOT  (AT_VLAST + 1)   /* total number additional energy items */
#define AT_V0    2                /* starting index for additional energy items */
#define AT_LAMBDA (AT_V0 + AT_VSPB)
#define AT_ETOT  (AT_V0 + AT_VTOT)
#define AT_EFG   0
#define AT_EBG   1
#endif

/* define a long long int type */
#ifndef I32ONLY
typedef long long  llong_t;
#define llong_pfmt "%lld"
#else
typedef long llong_t;
#define llong_pfmt "%ld"
#endif

/* for object-generator */
#define USE_MPI GMX_MPI

typedef struct {
  double bmin;        /* minimal beta (highest temperature);
                         $key: beta_min;  $def: 0.20;  $must; $valid: @bmin >= 0.0; */
  double bmax;        /* maximal beta (lowest temperature);
                         $key: beta_max;  $def: 0.41;  $must; $valid: @bmax > @bmin; */
  double T0;          /* thermostat temperature; $key: T0;  $def: 300.0;  */
  double beta;        /* current temperature; $def: @bmax;  $io:none; */
  int    nsttemp;     /* frequency of tempering, 0: disable, -1: only ns; $def: -1;  */
  int    mvreps;      /* number of repeating Langevin eq.
                         $key: move_repeats; $def: 1;  */
  double tmstep;      /* MD integration step, for convenience;
                         $def: 0.002; $usr:cfg; $io:; */
  char   *rng_file;   /* file name of random number state;
                         $def: NULL; $closecall: if ($ismaster) mtsave(@@); */
  int    nsttrace;    /* interval of writing trace file; -1: only at ns, 0: disable; $def: -1;  */
  char   *trace_file; /* name of trace file; $def: "TRACE";
                         $closecall: if ($ismaster) log_close(@log); */
  logfile_t *log;     /* logfile $cfg:none; */
  int    premode;     /* preparation mode, e.g. for collecting p.m.f.,
                         mode 0: (disable, production run) freeze stat., use bin distref;
                         mode 1: update data, use cfg distref; (for a md2 prep. run)
                         mode 2: freeze data., use cfg distref;
                         mode 3: update data, use bin distref;
                         mode 4: freeze data, use bin distref;
                         $usr:cfg; */
  double boltz;       /* Boltzmann constant $usr:cfg; */

  /* tilted Hamiltonian: H = H0 + [c0 + a0 * (T/Tref - 1)] * H1 */
  int    th_mode; /* $key: boost_mode;      $def: 1;     0: disable; 1:enable;  */
  /* $prereq := @th_mode; */
  double th_Tref; /* $key: boost_Tref;      $def: 300.0; */
  double th_a0;   /* $key: boost_alpha0;    $def: 1.0;   */
  double th_amin; /* $key: boost_alpha_min; $def: -1.0;  */
  double th_amax; /* $key: boost_alpha_max; $def: 1e4;   */
  double th_f;    /* f(beta); $def: 0.0;  $io:none; */
  double th_fd;   /* (beta * f(beta))' = beta * f'(beta) + f(beta); $def: 0.0;  $io:none; */
  /* $prereq := $0; */
  mb_t *mb;       /* handle for multiple-bin estimator
                     $obj; $objfprefix: mb_; $clr;
                     $cfgargs: @bmin, @bmax, @boltz; */
#if AT_VER == 2
  int    mimic_md1;      /* mimic md1 (a fallback option); $def: 0;  */
  double H0;             /* original Hamiltonian, pot. energy; $def: 0.0;  $io:none; */
  double H1;             /* additional energy; $def: 0.0;  $io:none; */
  double Ea;             /* (beta*H)' = H0 + [beta*f'(beta) + f(beta)]*H1; $def: 0.0;  $io:none; */
  double Eeff;           /* H = H0 + f(beta)*H1, f(beta) == th_f; $def: 0.0;  $io:none; */
  real   vacc[AT_VTOT];  /* various energy items; $def: 0.0f;  */
  real   scale[AT_ETOT]; /* index 0: accelerated part, 1: unaccelerated part; $def: 1.0f;  */

  int spb_dopro;   /* include proline; $def: 0;  */
  int spb_dogly;   /* include glycine; $def: 0;  */
  int spb_doter;   /* include terminal caps; $def: 1;  */
  int spb_collect; /* update spb data, currently only in preruns;
                      true for 1, 3, 5, ...
                      $def = (@premode % 2);  $io:none; */
  int spb_nstdata; /* number of steps to write data
                      $key: spb_nststat; NOTE: old name for compatibility
                      $def = 1; $io:c; $valid: @@ > 0; */
  int spb_biased;  /* apply a force component from spb's pmf,
                      $def = (@premode == 0 || @premode <= 8); $io:none */
  int spb_binref;  /* use distref from binary file (for a continuation run)
                      true for premode = 0, 3, 4, 7
                      in regular md2 prep. run, binref = 0, so it doesn't matter
                      $def = ((@premode + 3) % 4 >= 2);  $io:none; */
  int spb_initclr; /* clear initial data (for a continuation run)
                      true for premode = 0, 5, 6, 7
                      $def = ((@premode > 0) && ((@premode + 7) % 8 >= 4)); $io:none; */
  int spb_tscal;   /* allow temperature scaling during preparation run, $def: 1; $io:c; */
  int spb_onlyres; /* if positive, index (from 1) of the only residue
                      to be included for spb; $def = 0; $io:c; */
  int spb_dihends; /* dihedral gradient from i-l atoms only $def = 0; $io:c; */
  spbonds_t *spbs;  /* handle for spbs,
                     $obj; $objfprefix: spbs_; $clr;
                     $cfgargs: (@premode > 0) ? "pre_" : NULL;
                     $rbargs: 0u; */

  /* backbone pair */
  int bb_dopro;   /* $def: 0; */
  int bb_dogly;   /* $def: 0; */
  int bb_nstdata; /* number of steps to write data
                     $def = 1; $io:c; $valid: @@ > 0; */
  int bb_nstflush;/* number of steps for MPI synchronization,
                     i.e., hist_l[] --> hist[], an expensive operation
                     $def = 100; $io:c; $valid: @@ > 0; */
  int bb_collect; /* update bb data in preruns;
                     $def = (@premode % 2); $io:none; */
  int bb_biased;  /* apply a force component from bb's pmf
                     $def = (@premode == 0 || @premode <= 8); $io:none; */
  int bb_tscal;   /* allow temperature scaling during preparation run
                     $def: 1; $io:c; */
  int bb_onlyres; /* if positive, index (starting from 1) of the only residue
                     to be included, $def = 0; $io:c */
  bbs_t *bbs; /* backbone pairs, $obj; $clr;
                 $cfgargs: (@premode > 0) ? "pre_" : NULL; */
#endif
  double nsteql; /* time for equilibration; $io:c; $def:0; */

#if AT_VER == 1
  double Ea;    /* total potential energy; $def: 0.0;  $io:none; */
  double scale; /* ratio of the thermostat temperature to the current one; $def: 1.0;  $io:none; */
#endif
} at_t; /* $cfgopen; $rb:0; $wb:0; $reduce:0; $bcast:0; */

/* convert between beta and T */
static double at_beta_T(at_t *at, double x) { return 1.0 / (at->boltz * x); }

#if  AT_VER == 2
void at_clear_energy(at_t *at)
{
  int i;
  for (i = 0; i < AT_VTOT; i++) { at->vacc[i] = 0.0f; }
}
#endif

/* update at->scale, after T is updated */
void at_updscl(at_t *at)
{
  double Tnow = at_beta_T(at, at->beta); /* current temperature */
#if AT_VER == 2
  static int updcnt;

  /* effective Hamiltonian = H0 + f(beta) H1, where f(beta) is the at->th_f */
  if (at->th_mode == 0) {
    at->th_f = 0;
    at->th_fd = 0;
  } else {
    at->th_f = at->th_a0 * (Tnow / at->th_Tref - 1);
    at->th_fd = - at->th_a0;
  }
  if (at->th_f < at->th_amin) {
    at->th_fd = at->th_f = at->th_amin;
  }
  if (at->th_f > at->th_amax) {
    at->th_fd = at->th_f = at->th_amax;
  }
  at->scale[AT_EBG] = (real) (at->T0/Tnow);
  if (at->mimic_md1) {
    at->scale[AT_EFG] = at->scale[AT_EBG];
    at->scale[AT_LAMBDA] = 0.0f;
  } else {
    /* for an energy item that is present both in H0 and H1 */
    at->scale[AT_EFG] = (real)(at->scale[AT_EBG] * (1.0 + at->th_f));
    at->scale[AT_LAMBDA]
      = (at->spbs->cnt > 0) ? (real)(at->th_f) : 0.0f;
  }

  updcnt++;
  if (updcnt < 3 || updcnt % 100000 == 0) {
    fprintf(stderr, "step: %d, tilted Hamiltonian: %g, %g, T: %g\n",
        updcnt, at->th_fd, at->th_f, at_beta_T(at, at->beta));
  }
#else
  at->scale = at->T0 / Tnow;
#endif
}


/* load previous data */
static int at_loaddata(at_t *at, int isctn)
{
  mb_t *mb = at->mb;
  int ver = 0, readmb = 1;

#if AT_VER == 2
  /* read distref from spb.bin even if it's not a continuation
   * for a cont. run, loaddata
   * for init. and spb_binref, loaddata
   * for init. and !spb_binref, don't load data */
  if (!isctn && !at->spb_binref)
    return 0;
#endif
#if AT_VER == 1
  if (!isctn) /* initial run */
    return 0;
#endif
  readmb = (at->premode > 0) ? 0 : isctn;
  if (readmb) {
    if (mb_readbin(mb, mb->av_file, &ver) != 0) { /* read stat. */
      fprintf(stderr, "cannot load mb data from %s\n", mb->av_file);
      return 1;
    }
    fprintf(stderr, "%s version: %d\n", mb->av_file, ver);
    mb_wze(mb, "ze_init");
    at->beta = at->mb->beta; /* update the current temperature */

    if (mb_eh_readbin(mb, mb->eh_file, &ver) != 0) { /* read E-histograms */
      fprintf(stderr, "cannot load energy histogram from %s\n", mb->eh_file);
      return 2;
    }
  }

#if AT_VER == 2
  if (at->spbs->cnt > 0) { /* read the spb data */
    unsigned rflags = SPB_VERIFY | (at->spb_binref ? 0 : SPB_VERREF);
    /* we do not retain the mf from pre-run for a regular run
     * this maybe is unnecessary for the current prerun/regular division */
    if (spbs_readbin(at->spbs, at->spbs->bin_file, &ver, rflags) != 0) {
      fprintf(stderr, "cannot load previous spb data\n");
      return 3;
    }
    if (!isctn && at->spb_initclr) spbs_clear(at->spbs);
    spbs_write(at->spbs, "spb_init.txt", 2, /* as is */ 1);
  }
  /* bb and spb are exclusive */
  else if (at->bbs->cnt > 0) {
    unsigned rflags = BB_VERIFY;
    if (bbs_readbin(at->bbs, at->bbs->bin_file, &ver, rflags) != 0) {
      fprintf(stderr, "cannot load previous bb data\n");
      return 4;
    }
    bbs_write(at->bbs, "bb_init.txt", 0);
  }
#endif
  fprintf(stderr, "successfully load previous data\n");
  return 0;
}

/* dump at_t to file, for an array, only the first and last
 * `arrmax' elements are printed */
int at_dump(at_t *at, const char *fname, int arrmax)
{
  FILE *fp;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  at_manifest(at, fp, arrmax);
  fclose(fp);
  return 0;
}

/* Basic initialization steps (master only)
 * 1. load configuration file
 * 2. load data
 * It should be called even during prerun */
at_t *at_open(const char *fname, int isctn, double tmstep, int premode)
{
  at_t *at;
  double beta0;

  /* this will also initialize settings for member objects such as at->mb */
  at = at_cfgopen((fname != NULL) ? fname : "at.cfg", tmstep, premode, BOLTZ);
  if (at == NULL) {
    fprintf(stderr, "failed to load configuration file.\n");
    return NULL;
  }
  at->log = log_open(at->trace_file);
  /* set an attempt value for at->beta before reading the value from mb.av */
  beta0 = at_beta_T(at, at->T0);
  at->beta = beta0;
  if (at->beta > at->bmax - 1e-5) at->beta = at->bmax - 1e-5;
  if (at->beta < at->bmin + 1e-5) at->beta = at->bmin + 1e-5;

  /* we only load previous data if it's continuation */
  if (at_loaddata(at, isctn) != 0) {
    at_close(at);
    return NULL;
  }

#if AT_VER == 2
  if (at->premode > 0) {
    die_if ( (at->spbs == NULL || at->spbs->cnt <= 0)
          && (at->bbs == NULL || at->bbs->cnt <= 0),
      "both spb and bb are disabled in premode %d\n", at->premode);
    at->beta = (at->spb_tscal ? at->bmin : beta0);
    fprintf(stderr, "premode: set T to %g\n", at_beta_T(at, at->beta));
  }
  /* ensure spb force is used during normal tempering */
  if (premode == 0) {
    if (at->spbs->cnt > 0)
      die_if(!at->spb_biased, "spb must be active in tempering\n");
    else if (at->bbs->cnt > 0)
      die_if(!at->bb_biased, "bb must be biased in tempering\n");
  }
#endif
  at_updscl(at); /* update scale since T might be changed */
  at_dump(at, "at.manifest", 3);
  return at;
}

/* don't write at the first step; write at the last step */
static int doevery(llong_t step, int nst, int bfirst, int blast)
{
  return !bfirst && (blast || (nst > 0 && step % nst == 0));
}

/* write various output files */
static void at_output(at_t *at, llong_t step,
    int ib, double invw, double t1, double t2, double Eav,
    int bfirst, int blast, int btr)
{
  int btrace;

  /* write the trace file */
  if (at->nsttrace > 0)
    btrace = (step % at->nsttrace == 0) || bfirst || blast;
  else /* tracing is disabled if at->nsttrace == 0 */
    btrace = (at->nsttrace == 0) ? btr : 0;

  if (btrace) {
    log_printf(at->log, "%10.3f %5d %10.6f %12.3f %12.3f %10.6f %8.3f %8.5f",
      step * at->tmstep, ib, t2 - t1, at->Ea, Eav, at->beta, t1, invw);
#if AT_VER == 2
    log_printf(at->log, " %12.3f %12.3f", at->H0, at->H1);
#endif
    log_printf(at->log, "\n");
  }

  if (doevery(step, at->mb->av_nstsave, bfirst, blast)) { /* save averages */
    mb_write(at->mb);
    mb_wze(at->mb, NULL);
    mtsave(at->rng_file);
  }
  if (doevery(step, at->mb->eh_nstsave, bfirst, blast)) { /* save energy histograms */
    mb_eh_writebin(at->mb, at->mb->eh_file, 0);
    mb_eh_recon(at->mb, NULL);
  }
}

/* change beta and output */
int at_move(at_t *at, llong_t step, int bfirst, int blast, int btr)
{
  double invw = 1.0, T1 = 0., T2 = 0., Eav = 0., ndlnw;
  int ib = 0, rep;
#if AT_VER == 2
  double varr[AT_VTOT + 1] = {0.0};
  varr[0] = at->H0;
  varr[1] = at->H1;
#else
  double *varr = NULL;
#endif

  /* no need to update temperature in premode */
  if (at->premode > 0) return 0;

  /* update energy data, change at->beta */
  /* repeat several times to change the temperature */
  for (rep = 0; rep < at->mvreps; rep++) {
    /* 1. deposit the current energy and temperature */
    mb_add(at->mb, at->Ea, varr, at->beta, &ib, &invw, &ndlnw);

    /* 2. use Langevin equation to update the temperature */
    T1 = at_beta_T(at, at->beta);
    at->beta = mb_move(at->mb, at->Ea, at->beta, ib, ndlnw, &Eav);
    T2 = at_beta_T(at, at->beta);
  }

  if (doevery(step, at->mb->nstrefresh, 0, blast))
    mb_refresh_et(at->mb, 1);
  at_output(at, step, ib, invw, T1, T2, Eav, bfirst, blast, btr);
  return 0;
}



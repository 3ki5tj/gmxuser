/* multiple-bin estimators based on integral identities
 * Copyright (c) 2009, 2010 Cheng Zhang */
#ifndef MB_C__
#define MB_C__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_CFG
#define ZCOM_ENDN
#define ZCOM_RNG
#include "zcom.h"

/* raw data */
typedef struct { /* $io:= bt; $clr := 1; */
  double s, se, se2, se3;
} sm_t; /* $private; $mpi:0; */

/* multiple-bin estimator parameters */
typedef struct {
  double   bmin;  /* minimal beta (highest temperature); $usr:cfg; */
  double   bmax;  /* maximal beta (lowest temperature); $usr:cfg; */
  double   bdel;  /* bin size of beta; $key: beta_del;  $def: 0.0001;  $must; $valid: @bdel > 1e-6; */
  int      n;     /* number of temperature bins;
                     $def: (int)((@bmax - @bmin)/@bdel - 1e-5) + 1;
                     $io:bt; $rbverify; */
  double   *barr; /* temperature array; $cnt: @n+1;
                     $def: @bmin + i * @bdel;
                     $bincnt: @n; $# nasty fix for binary file
                     $binprev: @lgv_tot; */
  /* $assert: @-check_barr(@) == 0; check beta array */
  /* $call: @bmax = @bmin + @bdel * @n;  fix bmax to a bin boundary */
  double   beta;  /* current value of beta; $def: 0.5 * (@bmin + @bmax);
                     $io:b; $binprev: @cnt_dbl;
                     $valid: @@ >= @bmin && @@ <= @bmax; */
  int      m;     /* maximal number of bins in a window; $def: 0;
                     $io:bt; $binprev: @n; $rbverify; */
  int      order; /* order, should be 1; $key: mbest_order; $def: 1;
                     $io:cbt; $binprev: @m; $valid: @order == 1; */
  unsigned flags; /* combination of flags; $def: 0;  $io:bt; $binprev: @order; */
  int      bwmod; /* $key: mbest_mbin_mode; 0: d(beta) 1: dT/T  2: d(kT); $def: 1;
                     $valid: @@ >= 0 && @@ <= 2; */
  double   bwdel; /* delta lnT;  $key: mbest_delta_lnT;  $def: 0.05;  $must;
                     $cfgprereq: @bwmod == 1;
                     $valid: @@ > @bdel/pow(@bmin, 1.0); */
                  /* $altvar: @bwdel;
                     delta beta; $key: mbest_delta_beta; $def: 0.02;  $must;
                     $cfgprereq: @bwmod == 0;
                     $valid: @@ > @bdel/pow(@bmin, 0.0); */
                  /* $altvar: @bwdel;
                     delta kT;   $key: mbest_delta_kT;   $def: 0.1;   $must;
                     $cfgprereq: @bwmod == 2;
                     $valid: @@ > @bdel/pow(@bmin, 2.0); */

  /* $flagprefix := MB_; $kprefix := mbest_; */
#define MB_DAMP    0x00000001 /* $key: damp;        $def: 1;  use adaptive averaging; */
#define MB_CV      0x00000002 /* $key: needcv;      $def: 1;  compute heat capacity; */
#define MB_SYMWIN  0x00000004 /* $key: sym_mbin;    $def: 1;  use symmetrical window; */
#define MB_ONEBIN  0x00000020 /* $key: single_bin;  $def: 0;  use single bin estimator; */
#define MB_VERBOSE 0x00001000 /* $key: verbose;     $def: 1;  being verbose; */
#define MB_SBCORR  0x00002000 /* $key: sbcorr;      $def: 1;  include energy fluctuation correction due to a small bin width for internal energy, etc. */

  int       *js;    /* lower  boundary of asym. beta windows for ehat (asym.); $cnt: @n+1;  $def: 0; $io:; */
  int       *jt;    /* higher boundary of asym. beta windows for ehat (asym.); $cnt: @n+1;  $def: 0; $io:; */
  int       *jset;  /* lower  boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1);
                       $cnt: @n;  $def: 0; $io:; */
  int       *jtet;  /* higher boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1);
                       $cnt: @n;  $def: 0; $io:; */
  /* $assert: mb_mkwin(@, @bwmod, @bwdel, @js, @jt, @jset, @jtet) == 0;
   * setup the temperature windows (js, jt, jset, jtet)
   * need to setup damp and symwin before calling! */
  /* $call: @m = mb_maxwinspan(@); compute m as the maximal window span */

  int      nstrefresh;  /* interval of recalculating et for all temperature;
                           $key: nstrefresh; $def: 10000; */
  int      av_nstsave;  /* interval of writing mbav and ze files;
                           $key: nstav;  $def: 10000; */
  int      av_binary;   /* use binary format in mbav file;
                           $key: mbav_binary;  $def: 1; */
  char     *av_file;    /* name of mbav file; $key: mbav_file;  $def: "mb.av";  */
  char     *ze_file;    /* name of ze file; $key: ze_file;  $def: "ZE";  */
  int      wze_reps;    /* number of iterations before writing ze file; $def: 5; */
  double   *vis;        /* number of visits; $cnt: @n;  $def: 0.0;  $io:none; */
  double   totvis;      /* total number of visits, number of tempering;
                           $def: 0.0;  $io:b; $binprev: @beta; $valid: @@ > 0; */
  double   *winstot;    /* total of sum.s over a multiple-bine temperature window;
                           $cnt: @n;  $def: 0.0; $io:none; */
  double   boltz;       /* Boltzmann constant (given by user); $usr:cfg; */

  /* langevin equation */
  double   (*grand)(void);   /* function pointer to a gaussian random number generator; $def: &grand0;  $io:none; */
  double   lgv_dt;      /* time step for the temperature Langevin eq. $key: Tdt;  $def: 1e-5;  */
  double   lgv_dTmax;   /* maximal amount of temperature change in a step; $key: dTmax;  $def: 25.0;  */
  double   lgv_rej;     /* number of attempts of langevin equation trying to change beta too drastically;
                           $def: 0.0;  $io:none; $clr; */
  double   lgv_rate;    /* $usr:bintmp; $binprev: @shk_base; $io:b;
                           $def:(@lgv_tot > 1.0) ? (@lgv_rej/@lgv_tot) : 0.0;
                           $rbdef: 0.0; */
  double   lgv_tot;     /* total number of attempts of using langevin equation;
                           $def: 0.0;  $io:b; $binprev: @lgv_rate; $clr;
                           $rbvalid: (@lgv_rej=@lgv_tot*lgv_rate) >= 0.0; */

  int      regl;        /* average within a bin first;
                           $key: @<regularize; $def: 2; */
  double   fracmin;     /* minimal allowable coefficient during left/right combination;
                           $def: 0.0;  */
  double   cvshiftmax;  /* maximal fraction for shift energy fluct.
                           if cv is monotonic, it should be 0.0,
                           for ising model, it can restrain the magnitude
                           $def: 1.0;  */
  /* $kprefix := $0; */

  /* ensemble parameters:
   *  w = [1 + amp*exp^{-0.5 ((beta - bfoc)/bdev)^2}] / beta^exp; */
  double ens_exp;   /* ensemble exponent of beta; $key: ensemble_factor;  $def: 1.0;  */
  double ens_amp;   /* ensemble focus amplitude;  $key: ensemble_amp;  $def: 0.0;  */
  double ens_bdev;  /* ensemble beta width; $key: ensemble_dbeta;  $def: 0.02;
                       $prereq: fabs(@ens_amp) > 1e-8;  */
  double ens_bfoc;  /* ensemble beta focus; $key: ensemble_focus;  $def: @bmax;
                       $prereq: fabs(@ens_amp) > 1e-8;  */
  double *ens_w;    /* array of ensemble weights at bin boundaries; $cnt: @n+1;
                       $def: 1.0/(@-ensinvw(@, @barr[i], NULL, NULL)); $io:none; */

  /* shrink parameters */
  double shk_base;    /* current generic shrink amplitude $def: 0.0; $io:b; $binprev: @totvis; */
  int    shk_winadj;  /* adjust shrink according to temperature window width;
                         $key: shrink_mbin_adjust;  $def: 1;  */
  double shk_max;     /* initial and maximal shrink (adjusted); $key: shrink0;
                         $def: 0.01;  $valid: @@ < 0.9 && @@ >= 0.0; */
  double *shk_gauge;  /* array used of modulation shrinking factors; $cnt: @n; $io:none;
                         $def: (@-ensinvw(@, 0.5*(@barr[i] + @barr[i+1]), NULL, NULL) * @m) / (@jtet[i] - @jset[i]);  */
  int    shk_mode;    /* 0: const, 1: amp/t, 2: amp/t^exp; $key: shrink_mode;  $def: 1;
                         $valid: @shk_mode >= 0 && @shk_mode <= 2; */
  double shk_min;     /* minimal value for enforcing acc. sampling; $key: shrinkmin;  $def: 0.0;  */
  int    shk_stop;    /* stop shrinking after this number of steps; $key: shrinkstop;  $def: -1;  */
  double shk_amp;     /* amp t^(-exp); $key: shrinkamp;  $def: 0.1;  $prereq: @shk_mode >= 1;  */
  double shk_exp;     /* amp t^(-exp); $key: shrinkexp;  $def: 1.0;  $prereq: @shk_mode >= 2;  */

  /* reconstructed averages $clr := 1; $io := ; */
  double   *lnz;    /* logarithm of the partition function; $cnt: @n+1;  $def: 0.0; */
  double   *ehat;   /* internal energy; $cnt: @n+1;  $def: 0.0;  */
  double   *cvhat;  /* heat capacity;   $cnt: @n+1;  $def: 0.0;  */
  double   *et;     /* bin-averaged internal energy; $cnt: @n;  $def: 0.0;
                       $io:b; $binprev: @barr; */
  double   *imbal;  /* |a+ - a-| / (a+ + a-) for left-right combination;
                       $cnt: @n+1;  $def: 0.0;  */
  int      *haset;  /* current et[i] is reasonably good;
                       $cnt: @n;    $def: 0; $io:none; */
  unsigned *qua;    /* bits represent whether estimated values are unbiased;
                       $cnt: @n+1;  $def: 0; $io:none;  */
  double   *ampf;   /* currently amplification factor for adaptive averaging;
                       $cnt: @n;    $def: 1.0;  */
  /* $io := $0; $clr := $0; */

  sm_t *sums;  /* normal data; $io: bt; $binprev: @et;
                  $obj; $cnt: @n;  $io:bt; $clr; */

#if MB_VER == 2
  /* $binprereq:= ver == 2 && @vcnt > 0; */
  int       vcnt;   /* number of additional energy items
                       $def: 2;  $io:bt; $rbverify;
                       $binprev: @sums; */
  sm_t      *vsums; /* data for additional energy items
                       $obj; $cnt: @vcnt, @n;  $clr;
                       $binprev: @vcnt; */
  double    *vb;    /* additional energy item reconstructed as function of beta
                       $cnt: @vcnt, (@n+1);  $def: 0.0;  $io:none; $clr; */
  /* $binprereq:= $0; */
#endif
  int has_xsums; /* $def: @flags & MB_DAMP; $io:b; $binprev: @flags; */
  int cnt_int; /* number of additional integer variables to be written to binary file
                  $io:b; $def=0; $binprev: @has_xsums; $rbverify; */
  int cnt_dbl; /* number of additional double variables to be written to binary file
                  $io:b; $def=5; $binprev: @cnt_int; $rbverify; */
  int       *idxcnt; /* index count; $cnt: @n;  $def: 0; $io:none; */
  /* compute idxcnt, and adjust m to the maximal value among idxcnt[j].
   * $call: @m = mb_maxidxcnt(@);
   * the call is needed before allocating midx, */
  int       *midx;    /* index look-up table;
                         $cnt: @n, @m;  $def: 0;
                         $io:none; */
  /*  $assert: mb_mkidx(@) == 0; setup indices  */
  sm_t      *xsums;  /* multiple-bin damping data;
                        $obj; $cnt: @n, @m;  $clr;
                        $wbprep: mb_normalize(@, -1);
                        $bin_jmin = @js[i]; $bin_jmax = @jt[i+1]; */

  /* energy histogram stuff, a new fold
   * $fold := eh; $fold_fprefix: mb_eh_; */
  /* $fold_bin_add: @n;
   * $fold_prereq: @eh_mode > 0;
   * $fold_valid:  @eh_mode == 1; */
  int     eh_mode;     /* 0: disable; 1: simple histogram; $key: ehist_mode;  $def: 0;
                          $valid: @@ == 0 || @@ == 1; */
  /* $prereq := @eh_mode; $# apply the test to all the following */
  int     eh_skip;     /* interval of reconstructing energy histograms;
                          $key: ehist_skip; $def: 10; $valid: @eh_skip > 0; */
  int     eh_bwmod;    /* 0: d(beta) 1: dT/T  2: d(kT)
                          $key: ehist_mbin_mode; $def = 1;
                          $valid: @@ >= 0 && @@ <= 2; */
  double  eh_bwdel;    /* delta lnT;  $key: ehist_delta_lnT;   $def: 0.05;
                          $must; $cfgprereq: @eh_mode && @eh_bwmod == 1;
                          $valid: @eh_bwdel > @bdel/pow(@bmin, 1.0); */
                       /* $altvar: eh_bwdel; delta beta
                          $key: ehist_delta_beta;  $def: 0.02;
                          $must; $cfgprereq: @eh_mode && @eh_bwmod == 0;
                          $valid: @eh_bwdel > @bdel/pow(@bmin, 0.0); */
                       /* $altvar: eh_bwdel; delta kT
                          $key: ehist_delta_kT;    $def: 0.10;
                          $must; $cfgprereq: @eh_mode && @eh_bwmod == 2;
                          $valid: @eh_bwdel > @bdel/pow(@bmin, 2.0); */
  double  eh_min;      /* minimal energy  $key: ehist_min; $def = -12.6e4;
                          $rbverify; $binprev: @eh_cnt; $io:cbt; */
  double  eh_max;      /* maximal energy $key: ehist_max; $def: -9.0e4;
                          $valid: @eh_max > @eh_min; */
  double  eh_del;      /* energy bin size  $key: ehist_del; $def = 20.0;
                          $io:cbt; $binprev: @eh_min;
                          $valid: @eh_del > 0; $rbverify; */
  int     eh_cnt;      /* number of energy bins; $io:bt; $rbverify;
                          $def: (int)((@eh_max-@eh_min)/@eh_del - 1e-5 + 1); */
  int     eh_binary;   /* binary format for ehist file;
                          $key: ehist_binary; $def: 1; */
  int     eh_nstsave;  /* interval of writing histogrm files;
                          $key: nsthist;  $def: 100000; */
  char    *eh_file;    /* name of ehist file;
                          $key: ehist_file; $def: "hist.bin";  */
  char    *eh_rfile;   /* name of reconstructed energy histogram;
                          $key: ehist_mbin_file; $def: "HMB";  */
  double  *eh_his;     /* energy histogram data;
                          $cnt: @n, @eh_cnt; $clr;
                          $bin_trim;  $# try to save space in writing binary */
  double  *eh_recon;   /* temporary space for reconstructing histogram
                          $cnt: @eh_cnt; $io:none; $clr; */
  int     *eh_is;      /* indices for temperature windows (lower);
                          $cnt: @n + 1;  $io:none; */
  int     *eh_it;      /* indices for temperature windows (higher);
                          $cnt: @n + 1;  $io:none;
  $cfgvalid: 0 == mb_mkwin(@, @eh_bwmod, @eh_bwdel, @eh_is, @eh_it, NULL, NULL);  */

  /* $flagprefix := MB_EH_; */
#define MB_EH_ADDAHALF 0x00010000  /* $key: ehist_addahalf;  $def: 1; add a half energy bin width in output; */
#define MB_EH_KEEPEDGE 0x00020000  /* $key: ehist_keepedge;  $def: 0; keep zero edge at sides; */
#define MB_EH_NOZEROES 0x00040000  /* $key: ehist_nozeroes;  $def: 0; do not output zeroes; */
} mb_t; /* $ptrname: mb; $mpi: 0; */

#define MBQ_ET    0x00000001  /*  et quality bit   */
#define MBQ_EHAT  0x00000002  /*  ehat quality bit */
#define MBQ_CV    0x00000004  /*  cv quality bit   */
#define MBQ_LNZ   0x00000008  /*  lnz quality bit  */

#define MB_LOOSE  0x00000010  /*  temporarily allow empty temperature windows */

/* multiply everything in sm by a factor `f' */
static void sm_mul(sm_t *sm, double f)
{
  sm->s *= f;
  sm->se *= f;
  sm->se2 *= f;
  sm->se3 *= f;
}

/* add a new energy to sm_t with a weight `w' */
static void sm_adde(sm_t *sm, double w, double e,
    double *de, double *dep, double *var, int cv)
{
  double s, sp, se, sep, see, seep, de2;

  sm->s = sp = (s = sm->s) + w;
  sm->se = sep = (se = sm->se) + e*w;
  if (s <= 0.0) return;
  *de = e - se/s;
  *dep = e - sep/sp;
  sm->se2 = seep = (see = sm->se2) + (de2 = (*de)*(*dep))*w;
  if (cv) sm->se3 += ((*var = de2-seep/s) - 2.0*see/s)*(*dep)*w;
}

#if MB_VER == 2
/* se   --> sv
 * se2  --> sve
 * se3  --> svee */
/* add an additional energy item */
static void sm_addv(sm_t *sm, double w, double v,
    double de, double dep, double var, int cv)
{
  double s, sp, svp, s2, dvp;

  sm->s = sp = (s = sm->s) + w;
  sm->se = svp = sm->se + v * w;
  if (s <= 0.0) return;
  dvp = v - svp/sp;
  sm->se2 = (s2 = sm->se2) + dvp*de*w;
  if (cv) sm->se3 += (dvp*var - 2.0*s2*dep/s)*w;
}
#endif

/* check if mb->barr is arranged in an ascending order */
static int mb_check_barr(mb_t *mb)
{
  int i;

  for (i = 0; i <= mb->n; i++)
    if (i > 0 && mb->barr[i] <= mb->barr[i-1]) {
      fprintf(stderr, "barr should ascend: barr[%d] = %g, barr[%d] = %g\n",
          i, mb->barr[i], i-1, mb->barr[i-1]);
      return 1;
    }
  return 0;
}

/* setup temperature windows, [ajs[i], ajt[i]) for each i in [0..n],
 * and symmetrical ones [ajset[i], ajtet[i]) for each i in [0..n) */
static int mb_mkwin(mb_t *mb, int bwmod, double bwdel,
    int ajs[], int ajt[], int ajset[], int ajtet[])
{
  int i, n, js, jt, idel, di1, di2, di;
  double bet, dbet = 0.0;

  die_if (mb == NULL || ajs == NULL || ajt == NULL,
      "null pointer mb: %p, ajs: %p, ajt: %p", mb, ajs, ajt);
  n = mb->n;
  for (i = 0; i <= n; i++) {
    bet = mb->barr[i];
    switch (bwmod) {
      case 0: dbet = bwdel; break;
      case 1: dbet = bwdel * bet; break;
      case 2: dbet = bwdel * (bet * bet); break;
      default: die_if(1, "bad bwmod=%d\n", bwmod);
    }
    idel = (int)(dbet/mb->bdel + 0.50000001);
    if (idel < 1) idel = 1; /* although idel = 0 should be fine */
    if ((js = i - idel) < 0) js = 0;
    if ((jt = i + idel) > n) jt = n;
    die_if (i < js || i > jt, "bad window %d (%d,%d)\n", i, js, jt);
    ajs[i] = js;
    ajt[i] = jt;
  }

  if (ajset == NULL || ajtet == NULL) return 0;
  /* calculate jset and jtet */
  for (i = 0; i < n; i++) {
    if (mb->flags & MB_SYMWIN) {
      di1 = i - ajs[i];
      di2 = ajt[i+1] - (i + 1);
      /* choose the smaller one from di1 and di2 */
      di  = (di1 > di2) ? di2 : di1;
      ajset[i] = i - di;
      ajtet[i] = i + 1 + di;
    } else {
      ajset[i] = ajs[ i ];
      ajtet[i] = ajt[ i + 1 ];
    }
  }
  return 0;
}

/* construct a more symmetrical bin for temperature mb->barr[i] than
 * (mb->js[i], mb->jt[i]) for quantities at a bin boundary, e.g. ehat */
static void mb_mksymwin(mb_t *mb, int i, int *js, int *jt)
{
  int j;

  j = (i < mb->n) ? i : (i - 1);
  (*js) = mb->jset[j];
  (*jt) = mb->jtet[j];
  if (i > 0 && i < mb->n) {
    if (i < mb->n/2) (*jt)--;
    else (*js)--;
  }
  die_if (*jt - *js <= 0, "empty window (%d,%d) for %d\n", *js, *jt, i);
}

/* compute the largest span of temperature windows in number of bins */
static int mb_maxwinspan(mb_t *mb)
{
  int i, cnt, max;

  /* first calculate the maximal # of bins, (attempted value) */
  for (max = 0, i = 0; i < mb->n; i++) {
    cnt = mb->jt[i + 1] - mb->js[i];
    if (cnt > max) max = cnt;
  }
  return mb->m = max;
}

/* adjust mb->m as the max. number of estimators affected by a single bin,
 * idxcnt needs to be allocated before calling */
static int mb_maxidxcnt(mb_t *mb)
{
  int i, j, js, jt, max;

  if (!(mb->flags & MB_DAMP)) return mb->m;
  /* mb->idxcnt[j]: number of estimators affected by stat. in bin j */
  for (i = 0; i < mb->n; i++) {
    js = mb->js[i];
    jt = mb->jt[i + 1];
    for (j = js; j < jt; j++)
      mb->idxcnt[j]++;
  }
  for (max = 0, j = mb->n; j; j--) /* compute max { mb->idxcnt[j] } */
    if (mb->idxcnt[j] > max) max = mb->idxcnt[j];
  if (max > mb->m) {
    if (mb->flags & MB_VERBOSE) printf("mb->m: %d => %d\n", mb->m, max);
    mb->m = max;
  }
  return mb->m;
}

/* set up matrix indices mb->midx, need mb->js, mb->jt */
static int mb_mkidx(mb_t *mb)
{
  int i, j, cnt, js, jt;

  if (!(mb->flags & MB_DAMP)) return 0;
  /* clear mb->idxcnt, whose space is to be used
   * to dynamically construct the matrix indices
   * mb->idxcnt will be restored at the end */
  for (j = 0; j < mb->n; j++) mb->idxcnt[j] = 0;
  for (i = 0; i < mb->n; i++) {
    js = mb->js[ i ];
    jt = mb->jt[ i + 1 ];
    for (j = js; j < jt; j++) {
      cnt = mb->idxcnt[j];
      die_if (cnt >= mb->m,
        "cnt: %d, m: %d, i: %d, j: %d, (js, jt) = (%d, %d).\n",
          cnt, mb->m, i, j, js, jt);
      die_if (j - js >= mb->m,
        "j: %d, js: %d, m: %d\n", j, js, mb->m);
      mb->midx[ j * mb->m + cnt ] = i * mb->m + j - js;
      mb->idxcnt[ j ] = cnt + 1;
    }
  }
  return 0;
}

/* normalize damping weight back to 1.0, averages are not affected
 * if i < 0, do for all estimators */
static void mb_normalize(mb_t *mb, int i)
{
  int i0, i1, j, js, jt;
  double fac;

  if (!(mb->flags & MB_DAMP)) return;
  if (i >= 0) { i0 = i1 = i; }
  else { i0 = 0; i1 = mb->n; } /* all estimators, if i < 0 */
  for (i = i0; i <= i1; i++) { /* loop over estimators */
    if (fabs(mb->ampf[i] - 1.0) < 1e-12) continue;
    fac = 1.0 / mb->ampf[i];
    mb->ampf[i] = 1.0;
    js = mb->js[i];
    jt = mb->jt[i+1];
    for (j = js; j < jt; j++)
      sm_mul(mb->xsums + i * mb->m + j - js, fac);
  }
}

/* return the reciprocal ensemble weight
 * wc = 1 + amp*exp{-0.5*[(beta-betfoc)/betdev]^2}, ndwc = -wc' */
double mb_ensinvw(mb_t *mb, double beta, double *wc, double *ndwc)
{
  const double eps = 1e-5;
  double x, dif, invvar, invw, wfocus, ndwfoc;
  int ifac;

  if (fabs(mb->ens_amp) > eps) {
    die_if (mb->ens_amp < 0.0, /* amplitude should be positive */
      "the ensemble amplitude %g must be positive!\n", mb->ens_amp);
    dif    = beta - mb->ens_bfoc;
    invvar = 1.0 / (mb->ens_bdev * mb->ens_bdev);
    wfocus = 1.0 + (x = mb->ens_amp * exp(-0.5 * (dif * dif) * invvar));
    ndwfoc = x * dif * invvar;
  } else {
    wfocus = 1.0;
    ndwfoc = 0.0;
  }
  if (wc   != NULL) *wc   = wfocus;
  if (ndwc != NULL) *ndwc = ndwfoc;

  beta /= mb->bmax; /* to relative beta */
  ifac = (int)(mb->ens_exp + 0.5); /* round to nearest int */
  if (fabs(mb->ens_exp - ifac) < 1e-5 && ifac >= 0) {
    for (invw = 1.0/wfocus; ifac > 0; ifac--)
      invw *= beta;
  } else  /* unable to avoid exponential */
    invw = exp(mb->ens_exp * log(beta)) / wfocus;

  die_if (invw > 1e6 || invw < 1e-6, "bad invw=%g, beta=%g\n", invw, beta);
  /* printf("calling with beta: %g, factor: %g, return: %g\n",
      beta, mb->ens_exp, invw); */
  return invw;
}

/* compute the temperature-independent shrinking factor */
static double mb_shkbase(mb_t *mb)
{
  double x, shk;

  if (mb->shk_stop >= 0 && mb->totvis > mb->shk_stop)
    return 0.0;
  else if (mb->shk_mode == 0)
    return mb->shk_min;
  shk = mb->shk_max;
  if (mb->totvis > 10 * mb->n) {
    x = mb->totvis / mb->n;
    if (mb->shk_mode == 1) {
      x = mb->shk_amp / x;
    } else if (mb->shk_mode == 2) {
      x = mb->shk_amp / pow(x, mb->shk_exp);
    } else die_if(1, "invalid shk_mode: %d\n", mb->shk_mode);
    if (x < shk) shk = x;
    if (shk < mb->shk_min) shk = mb->shk_min;
  }
  return mb->shk_base = shk;
}

/* compute the adjusted shk and invgam */
static double mb_invgam(mb_t *mb, int ib)
{
  double shk;

  shk = mb_shkbase(mb); /* compute the unadjusted */
  if (mb->shk_winadj) { /* multiply the gauge */
    die_if (mb->shk_gauge == NULL, "gauge is null\n");
    die_if (ib < 0 || ib >= mb->n, "index %d out of range\n", ib);
    shk *= mb->shk_gauge[ib];
    if (shk > mb->shk_max) shk = mb->shk_max;
  }
  return 1.0 / (1.0 - shk);
}

/* add energy and bet */
void mb_add(mb_t *mb, double e, const double v[], double bet,
    int *pib, double *pinvw, double *ndlnw)
{
  double de = 0.0, dep = 0.0, var = 0.0, ginvw, invw, wc = 1.0, ndwc = 0.0;
  int j, cv = mb->flags & MB_CV;

  *pib = j = (int)((bet - mb->bmin)/mb->bdel);
  die_if (j < 0 || j >= mb->n, "beta = %d, %g, range = (%g, %g, %g)\n",
      j, bet, mb->bmin, mb->bdel, mb->bmax);
  mb->vis[j] += 1.0;
  mb->totvis += 1.0;

  *pinvw = invw = mb_ensinvw(mb, bet, &wc, &ndwc); /* get weight */
  *ndlnw = mb->ens_exp/bet + ndwc/wc;
  sm_adde(mb->sums + j, invw, e, &de, &dep, &var, cv);

#if (MB_VER == 2)
  if (mb->vcnt > 0 && v != NULL) { /* add additional energy items */
    int iv;
    die_if (mb->vsums == NULL, "vsums is null, cnt: %d\n", mb->vcnt);
    for (iv = 0; iv < mb->vcnt; iv++)
      sm_addv(mb->vsums + iv*mb->n+j, invw, v[iv], de, dep, var, cv);
  }
#else
  die_if(v != NULL, "v[] must be null\n");
#endif

  if (mb->eh_mode > 0) { /* add to energy histogram */
    int ie = (int)((e - mb->eh_min)/mb->eh_del);
    if (ie >= 0 && ie < mb->eh_cnt) /* no invw for histogram */
      mb->eh_his[j*mb->eh_cnt+ie] += 1.0;
  }

  if (mb->flags & MB_DAMP) { /* add to damping data */
    int i, l, mid, cnt, upd;

    cnt = mb->idxcnt[j];
    for (l = 0; l < cnt; l++) { /* loop over affected estimators */
      mid = mb->midx[j * mb->m + l];
      i = mid / mb->m;
      die_if (i * mb->m + j - mb->js[i] != mid, /* check index */
        "index corruption, i=%d, m=%d, j=%d, js=%d, mid=%d(%d)\n",
            i, mb->m, j, mb->js[i], mid, i * mb->m + j - mb->js[i]);
      /* et is computed from (jset, jtet), which is a subset of (js, jt),
       * avoid weight updating when j lies outside of the former */
      if (mb->flags & MB_ONEBIN)
        upd = (j == i);
      else if (mb->flags & MB_SYMWIN)
        upd = (j >= mb->jset[i] && j < mb->jtet[i]);
      else upd = 1;
      /* apply adaptive averaging */
      if (upd) mb->ampf[i] *= mb_invgam(mb, i);
      ginvw = mb->ampf[i] * invw; /* multiply accumulated 1/gamma */
      de = 0.0;
      sm_adde(mb->xsums+mid, ginvw, e, &de, &dep, &var, cv);
      /* we call normalization when the weight starts to blow up */
      if (ginvw > 2.0) mb_normalize(mb, i);
    }
  }
}

/* compute total weighted number of visits to T. window of bin i */
static void mb_winstot(mb_t *mb)
{
  int i, j, js, jt, damp;
  double tot;
  sm_t *sm0;

  damp = (mb->flags & MB_DAMP);
  for (i = 0; i < mb->n; i++) {
    js = mb->jset[i];
    jt = mb->jtet[i];
    die_if (js < mb->js[i] || i < js || i >= jt,
      "bad window (%d, %d), js=%d, i=%d, n=%d\n", js, jt, mb->js[i], i, mb->n);
    sm0 = damp ? (mb->xsums + i*mb->m - js) : (mb->sums);
    for (tot = 0.0, j = js; j < jt; j++) tot += sm0[j].s;
    mb->winstot[i] = damp ? (tot/mb->ampf[i]) : tot;
  }
}

/* set quality bit */
static void mb_setqbit(unsigned *ptr, unsigned mask, int on)
{
  if (mask) { if (on) *ptr |= mask; else *ptr &= ~mask; }
}

/* translate quality bits into a 0-1 string */
static char *mb_qua2s(unsigned i)
{
  static char buf[64]; /* has to be static to be the return value */
  int  cnt = 0;

  buf[cnt++] = (char)((i & MBQ_ET  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_EHAT) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_CV  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_LNZ ) ? '1' : '0');
  buf[cnt] = '\0';
  return buf;
}

/* collect moments from i's left and i's right */
static void mb_lrcol(mb_t *mb, int i, int js, int jt,
    double t1[], double *tb, double s0[], double s1[],
    sm_t *sm0, int bcorr)
{
  double del, et0, el, var, x;
  int j, jx, lr, regl = mb->regl;

  die_if (0 > js || js >= jt || jt > mb->n, "bad window (%d, %d)", js, jt);
  s0[0] = s0[1] = s1[0] = s1[1] = t1[0] = t1[1] = 0.0;
  if (tb != NULL) *tb = 0.0;
  for (j = js; j < jt; j++) { /* loop over bins */
    sm_t *sm = sm0 + j;
    if (fabs(sm->s) < 1e-6) continue; /* empty bin */
    if (j < i) { jx = js; lr = 0; }
    else       { jx = jt; lr = 1; }
    /* correct energy fluctuation for small bin size:
     *   < (E - el)^2 > + (el - et0)^2,
     * el is the average energy of the bin */
    el = sm->se/sm->s;
    et0 = mb->haset[j] ? mb->et[j] : el;
    del = bcorr ? (el - et0) : 0.0;
    var = sm->se2/sm->s + del * del;
    if (regl) {
      s0[lr] += 1;
      s1[lr] += sm->se / sm->s;
      x = var;
    } else {
      s0[lr] += sm->s;
      s1[lr] += sm->se;
      x = var * sm->s;
    }
    t1[lr] += x * (j - jx + 0.5);
    if (j == i && tb != NULL) *tb = 0.5 * x;
  }
}

/* collect second order moments from i's left and i's right */
static void mb_lrcol2(mb_t *mb, int i, int js, int jt,
    double t2[], double s0[], double s2[], sm_t *sm0)
{
  double x;
  int j, jx, lr, regl = mb->regl;

  die_if (0 > js || js >= jt || jt > mb->n, "bad window (%d, %d)", js, jt);
  s0[0] = s0[1] = s2[0] = s2[1] = t2[0] = t2[1] = 0.0;
  for (j = js; j < jt; j++) { /* loop over bins */
    sm_t *sm = sm0+j;
    if (fabs(sm->s) < 1e-6) continue;
    if (j < i) { jx = js; lr = 0; }
    else       { jx = jt; lr = 1; }
    if (regl) {
      s0[lr] += 1;
      s2[lr] += sm->se2/sm->s;
      x = sm->se3/sm->s;
    } else {
      s0[lr] += sm->s;
      s2[lr] += sm->se2;
      x = sm->se3;
    }
    t2[lr] += x * (j - jx + 0.5);
  }
  t2[0] *= -mb->bdel;
  t2[1] *= -mb->bdel;
}

/* return an estimate from combining data from left and right
 * *imbal returns the difference between two coefficients
 * *good returns if an estimate is successful */
static double mb_lrbal(mb_t *mb, int ib, int js, int jt, int loose,
    const double t1[], double tb, const double s0[], const double s1[],
    double *imbal, int *good, unsigned qbit)
{
  double a, b;  /* combination coefficients, a for left, b for right */
  double num, den, del, tmp;
  int ip, im;

  a = b = 1.0;
  mb_setqbit(&mb->qua[ib], qbit, 0);

  if (fabs(s0[0] + s0[1]) < 1e-3) { /* no data in the entire window */
    if (good != NULL) *good = 0;
    die_if (!loose, "an empty window ib = %d\n", ib);
    return 0.0;
  }
  ip = jt - ib;
  im = ib - js;
  num = t1[0] + t1[1] + tb * (jt - js);
  den = t1[0] * ip - t1[1] * im; /* t1[0] >= 0 and t1[1] <= 0, so den >= 0 */
  if (fabs(den) < 1e-8) goto FALLBACK; /* visited at most once */
  del = num / den;
  if ((a = 1.0 - del*ip) < 0.0) a = 0.0;
  if ((b = 1.0 + del*im) < 0.0) b = 0.0;
  if (fabs(a + b) < 1e-8) goto FALLBACK;
  tmp = 1.0/(a + b);
  a *= tmp;
  b *= tmp;
  if (a < mb->fracmin) {
    b = 1.0 - (a = mb->fracmin);
  } else if (b < mb->fracmin) {
    a = 1.0 - (b = mb->fracmin);
  } else {
    mb_setqbit(&mb->qua[ib], qbit, 1);
  }
FALLBACK:
  /* in case denominator is 0.0, e.g., a*s0[0] + b*s0[1] = 1*0 + 0*3  */
  if (fabs(a * s0[0] + b * s0[1]) < 1e-6) {
    a = b = 1.0;
    mb_setqbit(&mb->qua[ib], qbit, 0);
  }
  if (imbal != NULL) *imbal = (a - b) / (a + b);
  if (good != NULL) *good = 1.0;
  return (a * s1[0] + b * s1[1]) / (a * s0[0] + b * s0[1]);
}

/* update the estimated energy et of the given ib
 * using a combination of data from left and right
 * note `flags' is different from mb->flags
 * if MB_LOOSE is in `flags', failure due to empty data is allowed */
static double mb_etlr(mb_t *mb, int ib, int flags)
{
  sm_t *sm0;
  int js, jt;
  double s0[2], s1[2], t1[2], tb;  /* first order */
  static int once;

  if (once++ == 0) fprintf(stderr, "etlr: fracmin = %g\n", mb->fracmin);

  js = mb->jset[ib];
  jt = mb->jtet[ib];
  die_if (js < 0 || jt <= js || jt > mb->n,
      "invalid indices %d, %d, ib = %d/%d", js, jt, ib, mb->n);

  /* choose xsums or sums */
  sm0 = (mb->flags & MB_DAMP) ? (mb->xsums + ib*mb->m-mb->js[ib]) : mb->sums;

  if (mb->flags & MB_ONEBIN || jt == js + 1) /* window reduced to a bin */
    return mb->et[ib] = (fabs(sm0[ib].s)>1e-6) ? (sm0[ib].se/sm0[ib].s) : 0.;

  /* using stat. from bins (js, jt) to form two estimates, left & right */
  mb_lrcol(mb, ib, js, jt, t1, &tb, s0, s1, sm0, 1);
  /* compute linear combination coefficients and et */
  return mb->et[ib] = mb_lrbal(mb, ib, js, jt, (flags & MB_LOOSE),
      t1, tb, s0, s1, &mb->imbal[ib], &mb->haset[ib], MBQ_ET);
}

/* compute the average energy Et at current bin ib */
double mb_calc_et(mb_t *mb, int ib, int flags)
{
  die_if (ib < 0 || ib >= mb->n, "bad ib %d [0, %d).\n", ib, mb->n);
  return mb_etlr(mb, ib, flags);
}

/* recompute all average energy */
void mb_refresh_et(mb_t *mb, int reps)
{
  int i, rep;

  for (rep = 0; rep < reps; rep++) {
    for (i = 0; i < mb->n; i++)
      mb_calc_et(mb, i, MB_LOOSE);
  }
}

static void mb_calc_ehat(mb_t *mb)
{
  int i, js, jt, needcv;
  sm_t *sm0;
  double bet, del, s0[2], s1[2], s2[2], t1[2], t2[2];
  static int once;

  if (once++ == 0) fprintf(stderr, "calc_ehat: fracmin = %g, cvshiftmax = %g\n",
        mb->fracmin, mb->cvshiftmax);

  needcv = (mb->flags & MB_CV);
  for (i = 0; i <= mb->n; i++) {
    mb_mksymwin(mb, i, &js, &jt);
    bet = mb->barr[i];
    if (mb->flags & MB_DAMP) { /* try to use damp sum if possible */
      int ip = (i == mb->n) ? (i-1) : i;
      sm0 = mb->xsums + ip * mb->m - mb->js[ip];
    } else sm0 = mb->sums;
    mb_lrcol(mb, i, js, jt, t1, NULL, s0, s1, sm0, mb->flags & MB_SBCORR);
    if (fabs(s0[0] + s0[1]) < 1e-3) { /* no data in the entire window */
      mb_setqbit(&mb->qua[i], MBQ_EHAT, 0);
      mb_setqbit(&mb->qua[i], MBQ_CV, 0);
      continue;
    }
    mb->ehat[i] = mb_lrbal(mb, i, js, jt, 1,
        t1, 0.0, s0, s1, NULL, NULL, MBQ_EHAT);

    if (needcv) { /* calculate Cv */
      mb_lrcol2(mb, i, js, jt, t2, s0, s2, sm0);
      if (t2[0] * t2[1] > 0) { /* t2[0] and t2[1] share the same sign */
        del = (t2[0] + t2[1]) / (s2[0] + s2[1]);
        if (del < -1) del = -1;
        if (fabs(del) > mb->cvshiftmax)
          del = mb->cvshiftmax * ((del > 0) ? 1.0 : -1.0);
        mb_setqbit(&mb->qua[i], MBQ_CV, 1);
        mb->cvhat[i] = mb->boltz * (bet * bet) *
                       (s2[0]+s2[1]) * (1+del) / (s0[0]+s0[1]);
        /* = (s2[0]+s2[1]+t2[0]+t2[1])/(s0[0]+s0[1])*bet*bet; */
      } else /* normal case */
        mb->cvhat[i] = mb->boltz * (bet * bet) *
          mb_lrbal(mb, i, js, jt, 1, t2, 0.0, s0, s2, NULL, NULL, MBQ_CV);
    }
  }
  /* estimate the partition function */
  for (mb->lnz[0] = 0.0, i = 0; i < mb->n; i++)
    mb->lnz[i+1] = mb->lnz[i] + mb->et[i] * (mb->barr[i] - mb->barr[i+1]);
}

#if (MB_VER == 2)
/* calculate other quantities */
static void mb_calc_hb(mb_t *mb)
{
  int i, js, jt, iv;
  double s0[2], s1[2], t1[2];
  static int once;

  if (mb->vcnt <= 0) return;
  if (once++ == 0) fprintf(stderr, "calc_hb: fracmin = %g\n", mb->fracmin);
  die_if (mb->vsums == NULL || mb->vb == NULL,
    "calc_hb: vsums:%p or vb:%p\n", mb->vsums, mb->vb);

  for (i = 0; i <= mb->n; i++) {
    mb_mksymwin(mb, i, &js, &jt);
    for (iv = 0; iv < mb->vcnt; iv++) {
      mb_lrcol(mb, i, js, jt, t1, NULL, s0, s1,
          mb->vsums + iv * mb->n, mb->flags & MB_SBCORR);
      mb->vb[ iv * (mb->n + 1) + i ] = mb_lrbal(mb,
          i, js, jt, 1, t1, 0, s0, s1, NULL, NULL, 0);
    } /* loop over different quantities */
  } /* loop over bins */
}
#endif

/* estimate `ergt' at the current temperature `bet'
 * return the new bet after integrating Langevin equation,
 * dt is the time step for kT or 1/beta; ndlnw = - d(lnw)/d(beta) */
static double mb_move(mb_t *mb, double erg, double bet, int ib,
    double ndlnw, double *ergt)
{
  double dkt, dktmax, kt1, kt2, dt, bet2, rndmag;

  *ergt = mb_calc_et(mb, ib, 0);
  dt = mb->lgv_dt * mb->boltz;
  dktmax = mb->lgv_dTmax * mb->boltz;
  kt1 = 1.0/bet;
  rndmag = kt1*sqrt(2.0*dt);
  die_if (mb->grand == NULL, "no gaussian RNG\n");
  dkt  = (erg - *ergt + ndlnw)*dt + rndmag * (*mb->grand)();
  if (dkt > dktmax) {
    dkt = dktmax;
    mb->lgv_rej += 1.0;
  } else if (dkt < -dktmax) {
    dkt = -dktmax;
    mb->lgv_rej += 1.0;
  }
  mb->lgv_tot += 1.0;
  kt2 = kt1 + dkt;
  bet2 = 1.0 / kt2;
  return mb->beta = (bet2 < mb->bmax && bet2 > mb->bmin) ? bet2 : bet;
}

/* write various averages to ze_file */
static int mb_wze(mb_t *mb, const char *fname)
{
  int i, ip;
#if MB_VER == 2
  int iv;
#endif
  FILE *fp;

  if (fname == NULL) fname = mb->ze_file;
  die_if (fname == NULL, "file name is NULL");

  mb_winstot(mb);
  mb_calc_ehat(mb);
#if MB_VER == 2
  die_if (mb->vcnt > 0 && mb->vb == NULL, "wze: mb->vb is NULL\n");
  mb_calc_hb(mb);
#endif
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot open %s.\n", fname);
    return 1;
  }
  for (i = 0; i <= mb->n; i++) {
    fprintf(fp, "%16.10f %20.10f %22.10f %22.10f ",
      mb->barr[i], mb->lnz[i], mb->ehat[i], mb->cvhat[i]);
    ip = (i < mb->n) ? i : (i-1); /* for quantities with no [mb->n] */
    fprintf(fp, " %22.10f %s %+10.6f %22.10e %22.10e %22.10e %22.10e",
      mb->et[ip], mb_qua2s(mb->qua[i]), mb->imbal[ip], mb->sums[ip].s,
      mb->vis[ip], mb->shk_gauge[ip], mb->winstot[ip]);
#if MB_VER >=2
    for (iv = 0; iv < mb->vcnt; iv++)
      fprintf(fp, " %22.10f", mb->vb[iv * (mb->n+1) + i]);
#endif
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* prepare and write mb data  */
int mb_write(mb_t *mb)
{
#if MB_VER == 2
  return mb_writebin(mb, mb->av_file, 2);
#else
  return mb_writebin(mb, mb->av_file, 1);
#endif
}

/* use the integral identity to reconstruct an unbiased histogram (robust method) */
static int mb_eh_recon(mb_t *mb, const char *fname)
{
  FILE *fp;
  int ib, j, js, jt, ie, imin, imax, cols, full, keep0;
  double eav, db, x;
  double num, den;
  double del, base, inc;

  if (mb->eh_mode == 0) return 0;
  die_if (mb->eh_mode != 1, "invalid eh_mode %d\n", mb->eh_mode);
  if ((fp = fopen((fname != NULL) ? fname : mb->eh_rfile, "w")) == NULL) {
    fprintf(stderr, "cannot write reconstructed histogram [%s].\n",
        mb->eh_rfile);
    return 1;
  }
  full = mb->flags & MB_EH_KEEPEDGE;
  keep0 = !(mb->flags & MB_EH_NOZEROES);
  del = (mb->flags & MB_EH_ADDAHALF) ? 0.5 : 0; /* for continuous system */
  cols = mb->eh_cnt;
  base = mb->eh_min;
  inc = mb->eh_del;
  db = mb->bdel;
  for (mb->lnz[0] = 0.0, ib = 0; ib < mb->n; ib++) /* build lnZ */
    mb->lnz[ib+1] = mb->lnz[ib] + mb->et[ib]*(mb->barr[ib] - mb->barr[ib+1]);

  /* loop over temperatures, and skip a few intermediate temperatures */
  for (ib = 0; ib <= mb->n; ib += mb->eh_skip) {
    /* reconstruct energy histogram at beta = mb->barr[ib] */
    js = mb->eh_is[ib];
    jt = mb->eh_it[ib];
    die_if(js < 0 || jt > mb->n || js >= jt, "bad window (%d, %d)\n", js, jt);
    /* loop over energy bins */
    for (ie = 0; ie < cols; ie++) {
      eav = base + (ie+del)*inc;
      for (den = 0, j = js; j <= jt; j++) { /* denominator */
        x = mb->ens_w[j] * exp(-eav*db*(j - ib) - mb->lnz[j] + mb->lnz[ib]);
        if (j == js || j == jt) x *= 0.5;
        den += x;
      }
      for (num = 0.0, j = js; j < jt; j++) { /* numerator */
        x = mb->eh_his[j*cols + ie];
        if (fabs(x) < 0.5) continue;
        num += x;
      }
      mb->eh_recon[ie] = num/den;
    }
    /* determine the output range */
    if (full) {
      imin = 0;
      imax = cols;
    } else { /* only print the smallest non-zero region */
      for (ie = cols-1; ie >= 0; ie--) if (mb->eh_recon[ie] > 0) break;
      if ((imax=ie+1) == 0) continue;
      for (imin = 0; imin < imax; imin++) if (mb->eh_recon[imin] > 0) break;
    }
    /* normalize the histogram, and output */
    for (x = 0, ie = imin; ie < imax; ie++) x += mb->eh_recon[ie];
    x = (fabs(x) > 1e-6) ? 1.0/(x*inc) : 1.0;
    for (ie = imin; ie < imax; ie++)
      if (keep0 || mb->eh_recon[ie] > 1e-6)
        fprintf(fp, "%g %.14E %g\n",
          base + (ie+del)*inc, mb->eh_recon[ie]*x, mb->barr[ib]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

#endif /* FILE */


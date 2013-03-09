/* multiple-bin estimators based on integral identities
 * Copyright (c) 2009, 2010 Cheng Zhang */
#ifndef MU_C__
#define MU_C__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_CFG
#define ZCOM_ENDN
#define ZCOM_RNG
#include "zcom2.h"

/* multiple-bin estimator parameters */
typedef struct {
  double  bmin;   /* minimal beta (highest temperature); $usr:cfg; */
  double  bmax;   /* maximal beta (lowest temperature); $usr:cfg; */
  double  bfrac;  /* fraction of beta to be included in mu updating $def: 0.1; $valid: @@ > 0; */
  double  mumin;  /* minimal mu; $key: mu_min; $def: 0; */
  double  mumax;  /* maximal mu; $key: mu_max; $def: 1; */
  double  mudel;  /* bin size of mu; $key: mu_del;  $def: 0.002;  $valid: @@ > 1e-6; */
  int     n;      /* number of temperature bins;
                    $def: (int)((@mumax - @mumin)/@mudel - 1e-5) + 1;
                    $io:bt; $rbverify; */
  double  *muarr; /* temperature array; $cnt: @n+1;
                     $def: @mumin + i * @mudel;
                     $bincnt: @n; $# nasty fix for binary file
                     $binprev: @shk_base; */
  /* $assert: @-check_muarr(@) == 0; check mu array */
  /* $call: @mumax = @mumin + @mudel * @n;  fix mumax to a bin boundary */
  double   mu;    /* current value of mu; $def: 0.5 * (@mumin + @mumax);
                     $io:b; $binprev: @cnt_dbl;
                     $valid: @@ >= @mumin && @@ <= @mumax; */
  int      m;     /* maximal number of bins in a window; $def: 0;
                     $io:bt; $binprev: @n; $rbverify; */
  int      order; /* order, should be 1; $key: mu_order; $def: 1;
                     $io:cbt; $binprev: @m; $valid: @order == 1; */
  unsigned flags; /* combination of flags; $def: 0;  $io:bt; $binprev: @order; */
  int      mwmod; /* $key: mu_mbin_mode; 0: d(mu) 1: dT/T  2: d(kT); $def: 0;
                     $valid: @@ >= 0 && @@ <= 2; */
  double   mwdel; /* delta lnT;  $key: mu_delta_lnmu;  $def: 0.05;
                     $cfgprereq: @mwmod == 1;
                     $valid: @@ > @mudel/pow(@mumin, 1.0); */
                  /* $altvar: @mwdel;
                     delta mu; $key: mu_delta_mu; $def: 0.05;
                     $cfgprereq: @mwmod == 0;
                     $valid: @@ > @mudel/pow(@mumin, 0.0); */
                  /* $altvar: @mwdel;
                     delta kT;   $key: mu_delta_invmu;   $def: 0.1;
                     $cfgprereq: @mwmod == 2;
                     $valid: @@ > @mudel/pow(@mumin, 2.0); */

  /* $flagprefix := MU_; $kprefix := mu_; */
#define MU_DAMP    0x00000001 /* $key: damp;        $def: 1;  use adaptive averaging; */
#define MU_CV      0x00000002 /* $key: needcv;      $def: 1;  compute heat capacity; */
#define MU_SYMWIN  0x00000004 /* $key: sym_mbin;    $def: 1;  use symmetrical window; */
#define MU_ONEBIN  0x00000020 /* $key: single_bin;  $def: 0;  use single bin estimator; */
#define MU_VERBOSE 0x00001000 /* $key: verbose;     $def: 1;  being verbose; */
#define MU_SBCORR  0x00002000 /* $key: sbcorr;      $def: 1;  include energy fluctuation correction due to a small bin width for internal energy, etc. */

  int       *js;    /* lower  boundary of asym. mu windows for vhat (asym.); $cnt: @n+1;  $def: 0; $io:; */
  int       *jt;    /* higher boundary of asym. mu windows for vhat (asym.); $cnt: @n+1;  $def: 0; $io:; */
  int       *jset;  /* lower  boundary of mu windows for vt (usu. sym.), i - jset[i] == jtet[i] - (i+1);
                       $cnt: @n;  $def: 0; $io:; */
  int       *jtet;  /* higher boundary of mu windows for vt (usu. sym.), i - jset[i] == jtet[i] - (i+1);
                       $cnt: @n;  $def: 0; $io:; */
  /* $assert: @-mkwin(@, @mwmod, @mwdel, @js, @jt, @jset, @jtet) == 0;
   * setup the temperature windows (js, jt, jset, jtet)
   * need to setup damp and symwin before calling! */
  /* $call: @m = @-maxwinspan(@); compute m as the maximal window span */

  int      av_nstsave;  /* interval of writing mbav and ze files;
                           $key: mu_nstav;  $def: 10000; */
  int      av_binary;   /* use binary format in mbav file;
                           $key: mu_avbinary;  $def: 1; */
  char     *av_file;    /* name of muav file; $key: mu_avfile;  $def: "mu.av";  */
  char     *zv_file;    /* name of zv file; $key: mu_zvfile;  $def: "ZV";  */
  int      wz_reps;     /* number of iterations before writing zu file; $def: 5; */
  double   *vis;        /* number of visits; $cnt: @n;  $def: 0.0;  $io:none; */
  double   totvis;      /* total number of visits, number of tempering;
                           $def: 0.0;  $io:b; $binprev: @mu; $valid: @@ > 0; */
  double   *winstot;    /* total of sum.s over a multiple-bine temperature window;
                           $cnt: @n;  $def: 0.0; $io:none; */

  /* langevin equation */
  double   (*grand)(void);   /* function pointer to a gaussian random number generator; $def: &grand0;  $io:none; */
  double   lgv_dt;      /* time step for the mu Langevin eq. $key: mu_dt;  $def: 1e-4;  */
  double   lgv_max;  /* maximal amount of temperature change in a step; $key: mu_mvmax;  $def: 0.1;  */

  int      regl;        /* average within a bin first;
                           $key: mu_regularize; $def: 2; */
  double   fracmin;     /* minimal allowable coefficient during left/right combination;
                           $def: 0.0;  */
  double   cvshiftmax;  /* maximal fraction for shift energy fluct.
                           if cv is monotonic, it should be 0.0,
                           for ising model, it can restrain the magnitude
                           $def: 1.0;  */
  /* $kprefix := $0; */

  /* ensemble parameters:
   *  w = 1.0 / mu^exp; */
  double ens_exp; /* ensemble exponent of mu; $key: mu_ensexp;  $def: 0.0;  */
  double ens_amp; /* ensemble focus amplitude;  $key: mu_ensamp;  $def: 0.0;  */
  double *ens_w;  /* array of ensemble weights at bin boundaries; $cnt: @n+1;
                       $def: 1.0/(@-ensinvw(@, @muarr[i])); $io:none; */

  /* shrink parameters */
  double shk_base;    /* current generic shrink amplitude $def: 0.0; $io:b; $binprev: @totvis; */
  int    shk_winadj;  /* adjust shrink according to temperature window width;
                         $key: mu_shkwinadj;  $def: 1;  */
  double shk_max;     /* initial and maximal shrink (adjusted); $key: mu_shk0;
                         $def: 0.01;  $valid: @@ < 0.9 && @@ >= 0.0; */
  double *shk_gauge;  /* array used of modulation shrinking factors; $cnt: @n; $io:none;
                         $def: (@-ensinvw(@, 0.5*(@muarr[i] + @muarr[i+1])) * @m) / (@jtet[i] - @jset[i]);  */
  int    shk_mode;    /* 0: const, 1: amp/t, 2: amp/t^exp; $key: mu_shkmode;  $def: 1;
                         $valid: @shk_mode >= 0 && @shk_mode <= 2; */
  double shk_min;     /* minimal value for enforcing acc. sampling; $key: mu_shkmin;  $def: 0.0;  */
  int    shk_stop;    /* stop shrinking after this number of steps; $key: mu_shkstop;  $def: -1;  */
  double shk_amp;     /* amp t^(-exp); $key: mu_shkamp;  $def: 0.1;  $prereq: @shk_mode >= 1;  */
  double shk_exp;     /* amp t^(-exp); $key: mu_shkexp;  $def: 1.0;  $prereq: @shk_mode >= 2;  */

  /* reconstructed averages $clr := 1; $io := ; */
  double   *lnz;    /* logarithm of the partition function; $cnt: @n+1;  $def: 0.0; */
  double   *vhat;   /* internal energy; $cnt: @n+1;  $def: 0.0;  */
  double   *chat;   /* heat capacity;   $cnt: @n+1;  $def: 0.0;  */
  double   *vt;     /* bin-averaged internal energy; $cnt: @n;  $def: 0.0;
                       $io:b; $binprev: @muarr; */
  double   *imbal;  /* |a+ - a-| / (a+ + a-) for left-right combination;
                       $cnt: @n+1;  $def: 0.0;  */
  int      *hasvt;  /* current vt[i] is reasonably good;
                       $cnt: @n;    $def: 0; $io:none; */
  unsigned *qua;    /* bits represent whether estimated values are unbiased;
                       $cnt: @n+1;  $def: 0; $io:none;  */
  double   *ampf;   /* currently amplification factor for adaptive averaging;
                       $cnt: @n;    $def: 1.0;  */
  /* $io := $0; $clr := $0; */

  sm_t *sums;  /* normal statistics; $io: bt; $binprev: @vt;
                  $obj; $cnt: @n;  $io:bt; $clr; */

  int has_xsums; /* $def: @flags & MU_DAMP; $io:b; $binprev: @flags; */
  int cnt_int; /* number of additional integer variables to be written to binary file
                  $io:b; $def=0; $binprev: @has_xsums; $rbverify; */
  int cnt_dbl; /* number of additional double variables to be written to binary file
                  $io:b; $def=5; $binprev: @cnt_int; $rbverify; */
  int       *idxcnt; /* index count; $cnt: @n;  $def: 0; $io:none; */
  /* compute idxcnt, and adjust m to the maximal value among idxcnt[j].
   * $call: @m = @-maxidxcnt(@);
   * the call is needed before allocating midx, */
  int       *midx;    /* index look-up table;
                         $cnt: @n, @m;  $def: 0;
                         $io:none; */
  /*  $assert: @-mkidx(@) == 0; setup indices  */
  sm_t      *xsums;  /* multiple-bin damping statistics;
                        $obj; $cnt: @n, @m;  $clr;
                        $wbprep: @-normalize(@, -1);
                        $bin_jmin = @js[i]; $bin_jmax = @jt[i+1]; */
} mu_t; /* $ptrname: mb; $mpi: 0; */

#define MUQ_VT    0x00000001  /*  vt quality bit   */
#define MUQ_VHAT  0x00000002  /*  vhat quality bit */
#define MUQ_CV    0x00000004  /*  cv quality bit   */
#define MUQ_LNZ   0x00000008  /*  lnz quality bit  */

#define MU_LOOSE  0x00000010  /*  temporarily allow empty temperature windows */

/* check if mb->muarr is arranged in an ascending order */
static int mu_check_muarr(mu_t *mb)
{
  int i;

  for (i = 0; i <= mb->n; i++)
    if (i > 0 && mb->muarr[i] <= mb->muarr[i-1]) {
      fprintf(stderr, "muarr should ascend: muarr[%d] = %g, muarr[%d] = %g\n",
          i, mb->muarr[i], i-1, mb->muarr[i-1]);
      return 1;
    }
  return 0;
}

/* setup temperature windows, [ajs[i], ajt[i]) for each i in [0..n],
 * and symmetrical ones [ajset[i], ajtet[i]) for each i in [0..n) */
static int mu_mkwin(mu_t *mb, int mwmod, double mwdel,
    int ajs[], int ajt[], int ajset[], int ajtet[])
{
  int i, n, js, jt, idel, di1, di2, di;
  double mv, dmv = 0.0;

  die_if (mb == NULL || ajs == NULL || ajt == NULL,
      "null pointer mb: %p, ajs: %p, ajt: %p", mb, ajs, ajt);
  n = mb->n;
  for (i = 0; i <= n; i++) {
    mv = mb->muarr[i];
    switch (mwmod) {
      case 0: dmv = mwdel; break;
      case 1: dmv = mwdel * mv; break;
      case 2: dmv = mwdel * (mv * mv); break;
      default: die_if(1, "bad mwmod=%d\n", mwmod);
    }
    idel = (int)(dmv/mb->mudel + 0.50000001);
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
    if (mb->flags & MU_SYMWIN) {
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

/* construct a more symmetrical bin for temperature mb->muarr[i] than
 * (mb->js[i], mb->jt[i]) for quantities at a bin boundary, e.g. vhat */
static void mu_mksymwin(mu_t *mb, int i, int *js, int *jt)
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
static int mu_maxwinspan(mu_t *mb)
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
static int mu_maxidxcnt(mu_t *mb)
{
  int i, j, js, jt, max;

  if (!(mb->flags & MU_DAMP)) return mb->m;
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
    if (mb->flags & MU_VERBOSE) printf("mb->m: %d => %d\n", mb->m, max);
    mb->m = max;
  }
  return mb->m;
}

/* set up matrix indices mb->midx, need mb->js, mb->jt */
static int mu_mkidx(mu_t *mb)
{
  int i, j, cnt, js, jt;

  if (!(mb->flags & MU_DAMP)) return 0;
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
static void mu_normalize(mu_t *mb, int i)
{
  int i0, i1, j, js, jt;
  double fac;

  if (!(mb->flags & MU_DAMP)) return;
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

/* return the reciprocal ensemble weight */
double mu_ensinvw(mu_t *mb, double mu)
{
  double invw;
  int ifac;

  mu /= mb->mumax; /* to relative mu */
  ifac = (int)(mb->ens_exp + 0.5); /* round to nearest int */
  if (fabs(mb->ens_exp - ifac) < 1e-5 && ifac >= 0) {
    for (invw = 1.0; ifac > 0; ifac--)
      invw *= mu;
  } else  /* unable to avoid exponential */
    invw = exp(mb->ens_exp * log(mu));

  die_if (invw > 1e6 || invw < 1e-6, "bad invw=%g, mu=%g\n", invw, mu);
  /* printf("calling with mu: %g, factor: %g, return: %g\n",
      mu, mb->ens_exp, invw); */
  return invw;
}

/* compute the temperature-independent shrinking factor */
static double mu_shkbase(mu_t *mb)
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
static double mu_invgam(mu_t *mb, int ib)
{
  double shk;

  shk = mu_shkbase(mb); /* compute the unadjusted */
  if (mb->shk_winadj) { /* multiply the gauge */
    die_if (mb->shk_gauge == NULL, "gauge is null\n");
    die_if (ib < 0 || ib >= mb->n, "index %d out of range\n", ib);
    shk *= mb->shk_gauge[ib];
    if (shk > mb->shk_max) shk = mb->shk_max;
  }
  return 1.0 / (1.0 - shk);
}

/* add energy and mu */
void mu_add(mu_t *mb, double beta, double e, double mu,
    int *pib, double *pinvw, double *ndlnw)
{
  double de = 0.0, dep = 0.0, var = 0.0, ginvw, invw;
  int j, cv = mb->flags & MU_CV;

  *pib = j = (int)((mu - mb->mumin)/mb->mudel);
  die_if (j < 0 || j >= mb->n, "mu = %d, %g, range = (%g, %g, %g)\n",
      j, mu, mb->mumin, mb->mudel, mb->mumax);
  mb->vis[j] += 1.0;
  mb->totvis += 1.0;

  /* add stat. only if beta is close enough to bmin */
  if (beta > mb->bmin + mb->bfrac*(mb->bmax - mb->bmin))
    return;

  *pinvw = invw = mu_ensinvw(mb, mu); /* get weight */
  *ndlnw = mb->ens_exp/mu;
  sm_adde(mb->sums + j, invw, e, &de, &dep, &var, cv);

  if (mb->flags & MU_DAMP) { /* add to damping statistics */
    int i, l, mid, cnt, upd;

    cnt = mb->idxcnt[j];
    for (l = 0; l < cnt; l++) { /* loop over affected estimators */
      mid = mb->midx[j * mb->m + l];
      i = mid / mb->m;
      die_if (i * mb->m + j - mb->js[i] != mid, /* check index */
        "index corruption, i=%d, m=%d, j=%d, js=%d, mid=%d(%d)\n",
            i, mb->m, j, mb->js[i], mid, i * mb->m + j - mb->js[i]);
      /* vt is computed from (jset, jtet), which is a subset of (js, jt),
       * avoid weight updating when j lies outside of the former */
      if (mb->flags & MU_ONEBIN)
        upd = (j == i);
      else if (mb->flags & MU_SYMWIN)
        upd = (j >= mb->jset[i] && j < mb->jtet[i]);
      else upd = 1;
      /* apply adaptive averaging */
      if (upd) mb->ampf[i] *= mu_invgam(mb, i);
      ginvw = mb->ampf[i] * invw; /* multiply accumulated 1/gamma */
      de = 0.0;
      sm_adde(mb->xsums+mid, ginvw, e, &de, &dep, &var, cv);
      /* we call normalization when the weight starts to blow up */
      if (ginvw > 2.0) mu_normalize(mb, i);
    }
  }
}

/* compute total weighted number of visits to T. window of bin i */
static void mu_winstot(mu_t *mb)
{
  int i, j, js, jt, damp;
  double tot;
  sm_t *sm0;

  damp = (mb->flags & MU_DAMP);
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
static void mu_setqbit(unsigned *ptr, unsigned mask, int on)
{
  if (mask) { if (on) *ptr |= mask; else *ptr &= ~mask; }
}

/* translate quality bits into a 0-1 string */
static char *mu_qua2s(unsigned i)
{
  static char buf[64]; /* has to be static to be the return value */
  int  cnt = 0;

  buf[cnt++] = (char)((i & MUQ_VT  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MUQ_VHAT) ? '1' : '0');
  buf[cnt++] = (char)((i & MUQ_CV  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MUQ_LNZ ) ? '1' : '0');
  buf[cnt] = '\0';
  return buf;
}

/* collect moments from i's left and i's right */
static void mu_lrcol(mu_t *mb, int i, int js, int jt,
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
    et0 = mb->hasvt[j] ? mb->vt[j] : el;
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
static void mu_lrcol2(mu_t *mb, int i, int js, int jt,
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
  t2[0] *= -mb->mudel;
  t2[1] *= -mb->mudel;
}

/* return an estimate from combining statistics from left and right
 * *imbal returns the difference between two coefficients
 * *good returns if an estimate is successful */
static double mu_lrbal(mu_t *mb, int ib, int js, int jt, int loose,
    const double t1[], double tb, const double s0[], const double s1[],
    double *imbal, int *good, unsigned qbit)
{
  double a, b;  /* combination coefficients, a for left, b for right */
  double num, den, del, tmp;
  int ip, im;

  a = b = 1.0;
  mu_setqbit(&mb->qua[ib], qbit, 0);

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
    mu_setqbit(&mb->qua[ib], qbit, 1);
  }
FALLBACK:
  /* in case denominator is 0.0, e.g., a*s0[0] + b*s0[1] = 1*0 + 0*3  */
  if (fabs(a * s0[0] + b * s0[1]) < 1e-6) {
    a = b = 1.0;
    mu_setqbit(&mb->qua[ib], qbit, 0);
  }
  if (imbal != NULL) *imbal = (a - b) / (a + b);
  if (good != NULL) *good = 1.0;
  return (a * s1[0] + b * s1[1]) / (a * s0[0] + b * s0[1]);
}

/* update the estimated energy vt of the given ib
 * using a combination of statistics from left and right
 * if MU_LOOSE is in `flags', unknown is allowed */
static double mu_etlr(mu_t *mb, int ib, int flags)
{
  sm_t *sm0;
  int js, jt;
  double s0[2], s1[2], t1[2], tb;  /* first order */

  js = mb->jset[ib];
  jt = mb->jtet[ib];
  die_if (js < 0 || jt <= js || jt > mb->n,
      "invalid indices %d, %d, ib = %d/%d", js, jt, ib, mb->n);

  /* choose statistics */
  sm0 = (mb->flags & MU_DAMP) ? (mb->xsums + ib*mb->m-mb->js[ib]) : mb->sums;

  if (mb->flags & MU_ONEBIN || jt == js + 1) /* window reduced to a bin */
    return mb->vt[ib] = (fabs(sm0[ib].s)>1e-6) ? (sm0[ib].se/sm0[ib].s) : 0.;

  /* using stat. from bins (js, jt) to form two estimates, left & right */
  mu_lrcol(mb, ib, js, jt, t1, &tb, s0, s1, sm0, 1);
  /* compute linear combination coefficients and vt */
  return mb->vt[ib] = mu_lrbal(mb, ib, js, jt, (flags & MU_LOOSE),
      t1, tb, s0, s1, &mb->imbal[ib], &mb->hasvt[ib], MUQ_VT);
}

/* compute the average energy Et at current bin ib */
double mu_calcvt(mu_t *mb, int ib, int flags)
{
  die_if (ib < 0 || ib >= mb->n, "bad ib %d [0, %d).\n", ib, mb->n);
  return mu_etlr(mb, ib, flags);
}

static void mu_calcvhat(mu_t *mb)
{
  int i, js, jt, needcv;
  sm_t *sm0;
  double del, s0[2], s1[2], s2[2], t1[2], t2[2];

  needcv = (mb->flags & MU_CV);
  for (i = 0; i <= mb->n; i++) {
    mu_mksymwin(mb, i, &js, &jt);
    if (mb->flags & MU_DAMP) { /* try to use damp sum if possible */
      int ip = (i == mb->n) ? (i-1) : i;
      sm0 = mb->xsums + ip * mb->m - mb->js[ip];
    } else sm0 = mb->sums;
    mu_lrcol(mb, i, js, jt, t1, NULL, s0, s1, sm0, mb->flags & MU_SBCORR);
    if (fabs(s0[0] + s0[1]) < 1e-3) { /* no data in the entire window */
      mu_setqbit(&mb->qua[i], MUQ_VHAT, 0);
      mu_setqbit(&mb->qua[i], MUQ_CV, 0);
      continue;
    }
    mb->vhat[i] = mu_lrbal(mb, i, js, jt, 1,
        t1, 0.0, s0, s1, NULL, NULL, MUQ_VHAT);

    if (needcv) { /* calculate mu fluctuation */
      mu_lrcol2(mb, i, js, jt, t2, s0, s2, sm0);
      if (t2[0] * t2[1] > 0) { /* t2[0] and t2[1] share the same sign */
        del = (t2[0] + t2[1]) / (s2[0] + s2[1]);
        if (del < -1) del = -1;
        if (fabs(del) > mb->cvshiftmax)
          del = mb->cvshiftmax * ((del > 0) ? 1.0 : -1.0);
        mu_setqbit(&mb->qua[i], MUQ_CV, 1);
        mb->chat[i] = (s2[0]+s2[1]) * (1+del) / (s0[0]+s0[1]);
        /* = (s2[0]+s2[1]+t2[0]+t2[1])/(s0[0]+s0[1])*mu*mu; */
      } else /* normal case */
        mb->chat[i] =
          mu_lrbal(mb, i, js, jt, 1, t2, 0.0, s0, s2, NULL, NULL, MUQ_CV);
    }
  }
  /* estimate the partition function */
  for (mb->lnz[0] = 0.0, i = 0; i < mb->n; i++)
    mb->lnz[i+1] = mb->lnz[i] + mb->vt[i] * (mb->muarr[i] - mb->muarr[i+1]);
}

/* estimate `vav' at the current `mu'
 * return the new mu after integrating Langevin equation */
static double mu_move(mu_t *mb, double v, double mu, int ib,
    double ndlnw, double *vav)
{
  double dmu, max, mu2, dt, rndmag;

  *vav = mu_calcvt(mb, ib, 0);
  dt = mb->lgv_dt;
  max = mb->lgv_max;
  rndmag = sqrt(2.0*dt);
  die_if (mb->grand == NULL, "no gaussian RNG\n");
  dmu  = (*vav - v - ndlnw)*dt + rndmag * (*mb->grand)();
  if (dmu > max) {
    dmu = max;
  } else if (dmu < -max) {
    dmu = -max;
  }
  mu2 = mu + dmu;
  return mb->mu = (mu2 < mb->mumax && mu2 > mb->mumin) ? mu2 : mu;
}

/* write various averages to zv_file */
static int mu_wz(mu_t *mb, const char *fname)
{
  int i, ip;
  FILE *fp;

  if (fname == NULL) fname = mb->zv_file;
  die_if (fname == NULL, "file name is NULL");

  mu_winstot(mb);
  mu_calcvhat(mb);

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot open %s.\n", fname);
    return 1;
  }
  for (i = 0; i <= mb->n; i++) {
    fprintf(fp, "%16.10f %20.10f %22.10f %22.10f ",
      mb->muarr[i], mb->lnz[i], mb->vhat[i], mb->chat[i]);
    ip = (i < mb->n) ? i : (i-1); /* for quantities with no [mb->n] */
    fprintf(fp, " %22.10f %s %+10.6f %22.10e %22.10e %22.10e %22.10e\n",
      mb->vt[ip], mu_qua2s(mb->qua[i]), mb->imbal[ip], mb->sums[ip].s,
      mb->vis[ip], mb->shk_gauge[ip], mb->winstot[ip]);
  }
  fclose(fp);
  return 0;
}

/* prepare and write mb statistics  */
int mu_write(mu_t *mb)
{
  return mu_writebin(mb, mb->av_file, 2);
}

#endif /* FILE */


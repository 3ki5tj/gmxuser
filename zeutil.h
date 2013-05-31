#ifndef ZEUTIL_H__
#define ZEUTIL_H__

/* Utilities for multiple-histogram reweighting
 * used for md1 and md2
 * self-contained (depend only on the standard libraries, but not zcom.h)
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef xnew
#define xnew(x, n) \
  if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); }
#endif

#ifndef xrenew
#define xrenew(x, n) \
  if ((x = realloc(x, (n)*sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); }
#endif

#define LNZERO (-1e8)

#ifndef LNADD_DEFINED
#define LNADD_DEFINED
#define LNBIG 50.0
/* return log(exp(a) + exp(b)) */
static __inline double lnadd(double a, double b)
{
  double c;
  if (a > b) return ( ((c = a-b) > LNBIG) ? a : (a+log(1+exp(-c))) );
  else       return ( ((c = b-a) > LNBIG) ? b : (b+log(1+exp(-c))) );
}
#endif

/* structure data from ZE file */
typedef struct {
  int bcnt, ecnt;
  double *barr, *lnz, *et, *eav;
  double bmin, bmax, bdel;
  double emin, emax, edel;
  double *lneden; /* accumulated weighting factor */
} zedata_t;

/* load the partition function lnz from ZE file `fn' */
static zedata_t *ze_load(const char *fn, double edel)
{
  FILE *fp;
  int n;
  char s[1024];
  double beta, x, eav, e0, e1;
  zedata_t *ze;

  xnew(ze, 1);
  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fn);
    return NULL;
  }
  /* get the number of lines */
  for (n = 0; fgets(s, sizeof s, fp); n++) {
    sscanf(s, "%lf", &beta);
    if (n == 0) {
      ze->bmin = ze->bmax = beta;
      continue;
    }
    if (beta > ze->bmax) ze->bmax = beta;
    if (beta < ze->bmin) ze->bmin = beta;
  }
  rewind(fp); /* go back to the first line */
  ze->bcnt = n - 1;
  ze->bdel = (ze->bmax - ze->bmin)/ze->bcnt;
  printf("beta = (%g, %g) bcnt = %d\n", ze->bmin, ze->bmax, ze->bcnt);
  xnew(ze->barr, ze->bcnt + 1);
  xnew(ze->lnz, ze->bcnt + 1);
  xnew(ze->et, ze->bcnt + 1);
  xnew(ze->eav, ze->bcnt + 1);
  for (n = 0; fgets(s, sizeof s, fp);  n++) {
    if (3 != sscanf(s, "%lf%lf%lf", &beta, &x, &eav)) {
      fprintf(stderr, "cannot read lnz n = %d\n", n);
      fclose(fp);
      return NULL;
    }
    ze->barr[n] = beta;
    ze->lnz[n] = x;
    ze->eav[n] = eav;
  }
  fclose(fp);
  for (n = 0; n < ze->bcnt; n++)
    ze->et[n] = (ze->lnz[n] - ze->lnz[n+1])/ze->bdel;

  /* allocate energy data */
  e0 = ze->eav[ze->bcnt]; /* minimal energy */
  e1 = ze->eav[0]; /* maximal energy */
  ze->edel = edel;
  ze->emin = e0 - (e1-e0)*.1;
  ze->emin = ((int)(ze->emin/ze->edel))*ze->edel;
  ze->emax = e1 + (e1-e0)*.2;
  ze->emax = ((int)(ze->emax/ze->edel))*ze->edel;
  ze->ecnt = (int)((ze->emax - ze->emin)/ze->edel + .5);
  xnew(ze->lneden, ze->ecnt+1);
  printf("energy histogram (%g : %g : %g) ecnt = %d\n",
      ze->emin, ze->edel, ze->emax, ze->ecnt);
  printf("lnz is loaded from %s\n", fn);
  return ze;
}

/* return lnz(bet) */
static double ze_getlnz(zedata_t *ze, double bet)
{
  int ib;
  double db;

  ib = (int)((bet - ze->bmin)/ze->bdel - 1e-6);
  if (bet < ze->bmin || ib > ze->bcnt) {
    fprintf(stderr, "temperature %g too small, ib = %d\n", bet, ib);
  } else if (ib == ze->bcnt) {
    return ze->lnz[ib];
  }
  db = bet - ze->bmin - ze->bdel*ib;
  return ze->lnz[ib] - ze->et[ib]*db;
}

/* return eav(bet) */
static double ze_geteav(zedata_t *ze, double bet)
{
  int ib;
  double db, k;

  ib = (int)((bet - ze->bmin)/ze->bdel - 1e-6);
  if (bet < ze->bmin || ib > ze->bcnt) {
    fprintf(stderr, "temperature %g too small, ib = %d\n", bet, ib);
  } else if (ib == ze->bcnt) {
    return ze->eav[ib];
  }
  db = bet - ze->bmin - ze->bdel*ib;
  k = (ze->eav[ib+1] - ze->eav[ib])/ze->bdel;
  return ze->eav[ib] + k*db;
}

/* free ze */
void ze_free(zedata_t *ze)
{
  if (ze == NULL) return;
  free(ze->barr);
  free(ze->lnz);
  free(ze->et);
  free(ze->eav);
  free(ze->lneden);
  free(ze);
}

#define TBLOCK 64  /* block of frames */

typedef struct {
  double bet, erg, t, lnw;
  int in; /* a frame is included */
  int ie; /* energy index of the frame */
  int dataid;
} frame_t;

typedef struct {
 int n; /* number of frames */
 int nalloc;
 int pos; /* current position */
 int dataid;
 frame_t *fr;
} trj_t;

/* initialize a trajectory structure */
static trj_t *trj_init(void)
{
  trj_t *trj;

  xnew(trj, 1);
  trj->n = 0;
  trj->nalloc = TBLOCK;
  trj->pos = 0;
  xnew(trj->fr, trj->nalloc);
  return trj;
}

/* add a frame to the trajectory */
static int trj_add(trj_t *trj, double bet, double erg, double t)
{
  int i;
  frame_t *fr;

  i = trj->n;
  if (++(trj->n) > trj->nalloc) {
    trj->nalloc += TBLOCK;
    xrenew(trj->fr, trj->nalloc);
  }
  fr = trj->fr + i;
  fr->bet = bet;
  fr->erg = erg;
  fr->t = t;
  fr->lnw = LNZERO;
  fr->in = 1;
  fr->dataid = trj->dataid;
  return 0;
}

/* load trajectory from TRACE file */
static int trj_load(trj_t *trj, const char *fntr)
{
  FILE *fp;
  double t, dT, Ea, Eav, bet;
  int ib, nline;
  static char s[2048];

  if ((fp = fopen(fntr, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fntr);
    return -1;
  }
  for (nline = 0; fgets(s, sizeof s, fp); nline++) {
    if (6 != sscanf(s, "%lf%d%lf%lf%lf%lf",
        &t, &ib, &dT, &Ea, &Eav, &bet)) {
      fprintf(stderr, "file %s, line %d: cannot scan data\ns=%s",
          fntr, nline, s);
      return -1;
    }
    trj_add(trj, bet, Ea, t);
  }
  fclose(fp);
  return 0;
}

/* mark successive frames too close in times as invalid,
 * i.e. trj->fr[i].in = 0 */
static int trj_prune(trj_t *trj, double dt, double tequil, int verbose)
{
  int i;
  double tprev = trj->fr[0].t, hdt = .5*dt;
  frame_t *fr;

  for (i = 1; i < trj->n; i++) {
    fr = trj->fr + i;
    if (fr->t < tprev - 0.01) {
      fprintf(stderr, "frame %d time goes backwards t = %g -> %g\n",
          i, tprev, fr->t);
      exit(1);
    } else if (fr->t < tprev + hdt) {
      if (verbose) {
        fprintf(stderr, "drop frame %d, t = %g, tprev = %g\n",
            i, fr->t, tprev);
      }
      fr->in = 0;
    } else {
      fr->in = 1;
    }
    if (fr->t < tequil)
      fr->in = 0;
    tprev = fr->t;
  }
  return 0;
}

/* compute multiple-histogram weight for a reference temperature
 * save to trj->fr[i].lnw
 * the weight is a only a function of the energy 1.0/exp(ze->lneden[ie])
 * the function should be called after the trajectory is loaded */
static int trj_getweight(trj_t *trj, zedata_t *ze,
    double bref, double delbetmax, int verbose)
{
  int i, j, ie, n = trj->n, en = ze->ecnt, nfr;
  double w, wm, lnzref, bet, erg, lnz, emin, edel, sum;
  const double lnhuge = 1e5;

  for (j = 0; j < en; j++)
    ze->lneden[j] = -(lnhuge+1);
  lnzref = ze_getlnz(ze, bref);
  emin = ze->emin;
  edel = ze->edel;
  nfr = 0;
  for (i = 0; i < n; i++) {
    if (!trj->fr[i].in) continue; /* pre-skip an invalid frame */
    bet = trj->fr[i].bet;
    erg = trj->fr[i].erg;
    trj->fr[i].in = 0;
    if (erg < emin) continue;
    ie = (int)((erg - emin)/edel);
    trj->fr[i].ie = ie;
    if (ie >= en) continue;
    if (fabs(bet - bref) > delbetmax) continue;
    trj->fr[i].in = 1; /* include the frame */
    nfr++;
    lnz = ze_getlnz(ze, bet);
    w = (bref - bet)*erg + lnzref - lnz;
/*
    if (fabs(w) > LNBIG*3.0) {
      if (++nmsg < 1000)
        fprintf(stderr, "w = %g/%g is too large, bref = %g, bet = %g, erg = %g, -dlnz/db = %g\n",
          w, ze->lneden[ie], bref, bet, erg, (lnzref-lnz)/(bet-bref));
    }
*/
    ze->lneden[ie] = lnadd(ze->lneden[ie], w);
  }

  /* normalize energy weighting function, such that wmax = 1.0 */
  for (wm = -1e9, ie = i = 0; i < n; i++) {
    if (!trj->fr[i].in) continue;
    w = -ze->lneden[trj->fr[i].ie];
    if (w > wm) { wm = w; ie = trj->fr[i].ie; }
  }
  printf("lnwmax = %g, wmax = %g erg = %g, ie = %d\n",
      wm, exp(wm), ze->emin + ze->edel*ie, ie);
  for (ie = 0; ie < en; ie++)
    if (ze->lneden[ie] > -lnhuge) /* only for visited energy bins */
      ze->lneden[ie] += wm;
    else ze->lneden[ie] = lnhuge;  /* a not visited bin, just set as a large number for safe */

  /* compute the w at eav(bref) for reference */
  erg = ze_geteav(ze, bref);
  ie = (int)((erg - emin)/edel);
  w = ze->lneden[ie];
  printf("after correction: lnw(eav(bref)) to %g (w = %g), according to bref: %g, erg: %g, ie: %d\n",
      -w, exp(-w), bref, erg, ie);

  if (verbose) {
    /* print out the weight */
    FILE *fp = fopen("wt.dat", "w");
    if (fp != NULL) {
      for (ie = 0; ie < ze->ecnt; ie++)
        fprintf(fp, "%g %g\n", ze->emin + ie*ze->edel, -ze->lneden[ie]);
      fclose(fp);
    }
  }
  /* compute relative weight for each frame in trj */
  for (sum = 0., i = 0; i < n; i++) {
    if (!trj->fr[i].in) continue;
    trj->fr[i].lnw = w = -ze->lneden[trj->fr[i].ie];
    if (w > -10.0) sum += exp(w);
  }
  printf("effective sample size %g\n", sum);
  return nfr;
}

static void trj_free(trj_t *trj)
{
  if (trj == NULL) return;
  free(trj->fr);
  free(trj);
}

/* load a sequence of TRACE files */
static trj_t *trj_loadseq(const char *fntr, int *ndata, double dt,
    double tequil, zedata_t *ze, double bref, double dbet)
{
  trj_t *trj;
  int nfr, i;
  static char ftr[FILENAME_MAX];

  if ((trj = trj_init()) == NULL) {
    fprintf(stderr, "cannot init trajectory\n");
    return NULL;
  }
  /* loop over dataXX/TRACE */
  for (*ndata = -1, i = 1; ; i++) {
    sprintf(ftr, "data%d/%s", i, fntr);
    trj->dataid = i;
    if (trj_load(trj, ftr) != 0) {
      *ndata = i-1;
      break;
    }
  }
  if (*ndata <= 0) {
    fprintf(stderr, "no data folders\n");
    trj_free(trj);
    return NULL;
  }
  trj_prune(trj, dt, tequil, 0); /* remove duplicated frames */
  nfr = trj_getweight(trj, ze, bref, dbet, 0); /* compute weight function */
  printf("%d data folders, %d effective frames for beta = %g\n", *ndata, nfr, bref);
  return trj;
}

/* given a name abc, we search ./abc, ../abc, ../../abc
 * until an existing file is found
 * if not return NULL */
static char *nfexists(const char *fn, int level)
{
  char *path;
  int l, i;
  FILE *fp;

  xnew(path, strlen(fn)+level*3+1);
  for (l = 0; l <= level; l++) { /* l times ../ before file name */
    /* add prefix */
    for (path[0] = '\0', i = 0; i < l; i ++) {
      strcat(path, "../");
    }
    strcat(path, fn);
    if ((fp = fopen(path, "r")) != NULL) {
      fclose(fp);
      return path;
    }
  }
  free(path);
  return NULL;
}

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "xtcio.h"
#include "tpxio.h"     /* topology .tpr */
#include "statutil.h"
#include "physics.h"

/* helper function: read topology and configuration from .tpr file */
static gmx_mtop_t *read_tpx_conf(const char *fn, rvec **x,
    matrix box, int *nat, double *dt)
{
  t_tpxheader hdr;
  int ver, gen, step, natoms;
  real t, lambda;
  gmx_mtop_t *mtop = NULL;
  t_inputrec ir;

  read_tpxheader(fn, &hdr, 1, &ver, &gen);
  snew(*x, hdr.natoms);
  snew(mtop, 1);
#if defined(GMXVERSION) && (GMXVERSION >= 40099)
  read_tpx(fn,
      hdr.bIr ? &ir : NULL, box, &natoms,
      *x, NULL, NULL, mtop);
#else
  read_tpx(fn, &step, &t, &lambda,
      hdr.bIr ? &ir : NULL, box, &natoms,
      *x, NULL, NULL, mtop);
#endif
  *nat = natoms;
  *dt = (hdr.bIr) ? ir.delta_t : 0.002;
  return mtop;
}

#endif


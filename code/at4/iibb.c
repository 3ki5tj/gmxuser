/* integral identity for the joint distribution of two backbone dihedral angles
 * Copyright (C) 2010-2012 Cheng Zhang */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "md2bb.h"

#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"

const char *prog = "iimb";
const char *fnin = "pre_bb.bin"; /* raw input */
const char *fnbout = "iibb.bin";
const char *fntout = "iibb.txt";
const char *fnport = NULL; /* portable 2d distribution */
const char *atcfg = "at.cfg";
int verbose;
int ipver = 0;
int opver = 2;
int use_atcfg = 1;
int multibins = 10;
int dermbin = 0;  /* number of bins to average the force */
int normalize = 1;
int mfcomnz = 1; /* average mean force only to nonempty bins */
int polymer = 0;

/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-i", NULL, &fnin, "input file");
  argopt_add(ao, "-T", NULL, &fntout, "text output");
  argopt_add(ao, "-B", NULL, &fnbout, "binary output");
  argopt_add(ao, "-a", NULL, &fnport, "portable format");
  argopt_add(ao, "-m", "%d", &multibins, "window size for integral identity");
  argopt_add(ao, "-M", "%d", &dermbin, "window size for mean force");
  argopt_add(ao, "-p", "%d", &polymer, "polymer conversion needed");
  argopt_add(ao, "-v", "%b", &verbose, "verbose");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

/* initialize bs_t and load data */
static bbs_t *bbinit(void)
{
  bbs_t *bbs;
  unsigned read_flags = 0;

  xnew(bbs, 1);

  /* load basic information from configure first */
  if (use_atcfg) {
    cfg_t *cfg;

    if ((cfg = cfg_open(atcfg)) == NULL) {
      fprintf(stderr, "cannot read config %s.", atcfg);
      return NULL;
    }
    if (bbs_cfgopen_low(bbs, cfg, NULL) != 0) {
      return NULL;
    }
    cfg_close(cfg);

    read_flags |= BB_VERIFY;
    fprintf(stderr, "configuration file %s is used\n", atcfg);
  }

  /* binary file already has mf, use as is */
  read_flags |= BB_MF;
  if (bbs_readbin(bbs, fnin, &ipver, read_flags) != 0) {
    fprintf(stderr, "error during loading bb binary data\n");
    return NULL;
  }
  fprintf(stderr, "%s version %d\n", fnin, ipver);
  bbs_calcmfpmf(bbs, 0);
  return bbs;
}

INLINE int WRAP(int i, int j, int n)
{ return ((i + n) % n) * n + ((j + n) % n); }

/* recalculate the mean force for the modulating factor */
static int bbmf(bb_t *bb)
{
  int id, i, j, iw, ii, jj, ij;
  int n = bb->bins, n2 = bb->bins2;
  double bl = bb->bl, phi, psi, den;

  iw = dermbin/2;
  printf("average mean force for %dx%d bins\n", iw*2+1, iw*2+1);
  /* compute the mean force */
  for (id = 0; id < n2; id++) {
    i = id / n;
    j = id % n;
    phi = psi = den = 0;

    for (ii = i - iw; ii <= i + iw; ii++)
    for (jj = j - iw; jj <= j + iw; jj++) {
      ij = WRAP(ii, jj, n);
      den += bb->hist[ij];
      phi += bb->sbfphi[ij];
      psi += bb->sbfpsi[ij];
    }

    if (den > .1) {
      bb->mfphi[id] = phi/(bl*den);
      bb->mfpsi[id] = psi/(bl*den);
    } else {
      bb->mfphi[id] = 0.;
      bb->mfpsi[id] = 0.;
    }
  }
  return 0;
}

/* change dihedral to/from polymer convention, +/- pi */
static int pswitch(bb_t *bb)
{
  int i, j, id, i2, j2, id2, n, nh;
  double tmp;

  n = bb->bins;
  nh = n/2;

#define PSWAP(arr) \
  for (id = 0; id < (n*nh); id++) { \
    i = id / n; \
    i2 = (i < nh) ? (i + nh) : (i - nh); \
    j = id % n; \
    j2 = (j < nh) ? (j + nh) : (j - nh); \
    id2 = i2*n+j2; \
    tmp = arr[id], arr[id] = arr[id2], arr[id2] = tmp; }

  PSWAP(bb->hist);
  PSWAP(bb->sbfphi);
  PSWAP(bb->sbfpsi);
  PSWAP(bb->mfphi);
  PSWAP(bb->mfpsi);
  PSWAP(bb->pmf);
  PSWAP(bb->distref);
  PSWAP(bb->distref0);
  return 0;
}

/* integrate over a grid (th0, et0) to (th1, et1) */
static double delphase(double th0, double et0, double th1, double et1)
{
  double ret, s;

  ret = cos(et0 + et1 - (th0 + th1)) * sin(et1 - et0) * sin(th1 - th0)
      + .25 * (th1 - th0) * (sin(4.*et1) - sin(4.*et0));
  s = sin(2.*et1);
  if (fabs(s) > 1e-14)
    ret -= .5 * log(sin(et1 + th1)/sin(et1 + th0)) * (s*s);
  s = sin(2.*et0);
  if (fabs(s) > 1e-14)
    ret += .5 * log(sin(et0 + th1)/sin(et0 + th0)) * (s*s);
  return ret;
}

/* compute phi from x, sin(phi)^2 = x */
static double getphase(double di, double ispan)
{
  double x = di/ispan;

  if (x < 0.) {
    fprintf(stderr, "x must > 0, %g,  %g/%g\n",
        x, di, ispan);
    exit(1);
  }
  return asin(sqrt(x));
}

static double getpmfi(bb_t *bb, int ia, int ja, int jb)
{
  int i, j, ij, sgn, n = bb->bins;
  double v = 0., dx = bb->binw;

  if (ja < jb) {
    sgn = 1;
  } else {
    sgn = -1;
    j = ja; ja = jb; jb = j;
  }
  for (i = ia-1; i < ia+1; i++) {
    for (j = ja; j < jb; j++) {
      ij = WRAP(i, j, n);
      v -= bb->mfphi[ ij ] * dx;
    }
  }
  // 0.5 * (getpmf1(bb, ia, ja, ib+1, jb) + getpmf1(bb, ia, ja, ib-1, jb));
  return .5*v*sgn;
}

static double getpmfj(bb_t *bb, int ia, int ib, int ja)
{
  int i, j, ij, sgn, n = bb->bins;
  double v = 0., dx = bb->binw;

  if (ia < ib) {
    sgn = 1;
  } else {
    sgn = -1;
    i = ia; ia = ib; ib = i;
  }
  for (j = ja-1; j < ja+1; j++) {
    for (i = ia; i < ib; i++) {
      ij = WRAP(i, j, n);
      v -= bb->mfpsi[ ij ] * dx;
    }
  }
  // 0.5*(getpmf1(bb, ia, ja, ib, jb+1) + getpmf1(bb, ia, ja, ib, jb-1));
  return .5*v*sgn;
}

/* dpmf = pmf(ib, jb) - pmf(ia, ja) */
static double getpmf1(bb_t *bb, int ia, int ja, int ib, int jb)
{
  int i, j, is, js, ij, di, dj, n = bb->bins;
  double v = 0.0, dx = bb->binw;
  double xspan, yspan, del1, del2, th0, th1, et0, et1;

  /* a better interpolation is needed */
  if (ia == ib) {
    if (ja == jb) return 0.0;
    else return getpmfi(bb, ia, ja, jb);
  } else if (ja == jb) {
    return getpmfj(bb, ia, ib, ja);
  }

  di = (ib > ia) ? 1 : -1;
  dj = (jb > ja) ? 1 : -1;

  xspan = (ib - ia) * dx;
  yspan = (jb - ja) * dx;
  for (i = ia; (ib - i)*di > 0; i += di) {
    for (j = ja; (jb - j)*dj > 0; j += dj) {
      th0 = getphase(i - ia,      ib - ia);
      th1 = getphase(i - ia + di, ib - ia);
      et0 = getphase(j - ja,      jb - ja);
      et1 = getphase(j - ja + dj, jb - ja);
      del1 = delphase(th0, et0, th1, et1);
      del2 = delphase(et0, th0, et1, th1);
      is = (di > 0) ? i : (i-1);
      js = (dj > 0) ? j : (j-1);
      ij = WRAP(is, js, n);
      v -= del1 * bb->mfphi[ ij ] * xspan;
      v -= del2 * bb->mfpsi[ ij ] * yspan;
    }
  }
  return v;
}

/* handle a single dihedral */
static int dobb(bb_t *bb, int method)
{
  int i, j, id, ii, ia, ib, jj, ja, jb, ij;
  int n = bb->bins, n2 = bb->bins2, iw = multibins/2;
  double bl = bb->bl, sum = 0., fac, bn, num, den0, den, sum1;
  double dpmf;

  if ((iw = multibins/2) < 1) iw = 1;
  bbmf(bb);
  for (i = 0; i < n2; i++) {
    sum += bb->hist[i];
  }
  fac = 1.0/(sum * bb->binw * bb->binw);
  printf("bb %d: %g samples (iw = %d) \n", bb->id, sum, iw);
/*
  // testing case
  for (id = 0; id < n2; id++) {
    double val;
    i = id / n - 75;
    j = id % n - 20;
    bb->mfphi[id] = -0.3*(i+.5);
    bb->mfpsi[id] = -0.5*(j+.5);
  }
*/
  if (method == 1) { /* iteratively compute a self consistent mean force */
    fprintf(stderr, "bad method\n"); exit(1);
  }

  for (id = 0; id < n2; id++) {
    i = id / n;
    j = id % n;
    ia = i - iw;
    ib = i + iw;
    ja = j - iw;
    jb = j + iw;

    num = den = den0 = 0.;
    /* get the number of visits to the region */
    for (ii = ia; ii < ib; ii++)
    for (jj = ja; jj < jb; jj++) {
      ij = WRAP(ii, jj, n);
      num += bb->hist[ij];
    }
    for (ii = ia; ii <= ib; ii++)
    for (jj = ja; jj <= jb; jj++) {
      ij = WRAP(ii, jj, n);
      bn = 1.0;
      if (jj == ja || jj == jb) bn *= .5;
      if (ii == ia || ii == ib) bn *= .5;
      den0 += bn;
      if (method == 0) {
        dpmf = getpmf1(bb, i, j, ii, jj);
      } else {
        dpmf = bb->pmf[ij] - bb->pmf[id];
      }
      den += bn * exp(-dpmf * bl);
    }
    bb->distref[id] = num/den; /* unbiased result */
    bb->distref0[id] = num/den0; /* biased result */
    if (id % 10 == 0)
      printf("progress %8d/%-8d  %g%%  \r", id, n2, 100.0*id/n2);
  }
  for (sum1 = 0.0, i = 0; i < bb->bins2; i++)
    sum1 += bb->distref[i];
  sum1 *= 1.0/sum;
  printf("\nsum: %g\n", sum1);

  if (normalize) {
    double fac2 = fac;

    if (normalize == 2)
      fac2 = fac / sum1;
    for (i = 0; i < n2; i++) {
      bb->hist[i] *= fac;
      bb->sbfphi[i] *= fac;
      bb->sbfpsi[i] *= fac;
      bb->distref[i] *= fac2;
      bb->distref0[i] *= fac;
    }
  }
  if (polymer) pswitch(bb);
  return 0;
}

/* compute the correlation entropy */
static int calcetp(bb_t *bb)
{
  int i, j, n = bb->bins;
  double *rhox, *rhoy, dx = bb->binw, rind, rho;
  double S, Sx, Sy, Sxy;

  /* compute independent distribution along x and y */
  rhox = calloc(n, sizeof(double));
  rhoy = calloc(n, sizeof(double));
  if (rhox == NULL || rhoy == NULL) return -1;
  for (i = 0; i < n; i++) rhox[i] = rhoy[i] = 0.;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      rhox[i] += bb->distref[i*n + j];
      rhoy[j] += bb->distref[i*n + j];
    }
  for (i = 0; i < n; i++) {
    rhox[i] *= dx;
    rhoy[i] *= dx;
  }
  Sx = Sy = 0.0;
  for (i = 0; i < n; i++) {
    if (rhox[i] > 1e-14)
      Sx += -rhox[i]*log(rhox[i]);
    if (rhoy[i] > 1e-14)
      Sy += -rhoy[i]*log(rhoy[i]);
  }
  Sx = Sx*dx + log(2*M_PI);
  Sy = Sy*dx + log(2*M_PI);

  /* compute residue entropy */
  S = Sxy = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      rho = bb->distref[i*n + j];
      if (rho < 1e-14) continue;
      rind = rhox[i]*rhoy[j];
      S += -rho*log(rind/rho);
      Sxy += -rho*log(rho);
    }
  }
  S = S*dx*dx;
  Sxy = Sxy*dx*dx + log(4*M_PI*M_PI);
  printf("information entropy: %g, Sxy %g, Sx %g, Sy %g\n", S, Sxy, Sx, Sy);
  return 0;
}

/* compute the alpha/beta components */
static int calcab(bb_t *bb)
{
  int i, j, n = bb->bins;
  double phi, psi, dx = bb->binw, r;
  double cnt = 0., acnt = 0., bcnt = 0.;

  for (i = 0; i < n; i++)
  for (j = 0; j < n; j++) {
    phi = -M_PI + i*dx;
    psi = -M_PI + j*dx;
    r = bb->hist[i*n + j];
    if (phi < 0.) {
      if (psi < M_PI/2 && psi > -M_PI/2)
        acnt += r;
      else
        bcnt += r;
    }
    cnt += r;
  }
  printf("cnt = %g, alpha = %g%%, beta = %g%%\n",
      cnt*dx*dx, 100.*acnt/cnt, 100.*bcnt/cnt);
  return 0;
}

static int dobs(bbs_t *bbs)
{
  int i;

  for (i = 0; i < bbs->cnt; i++) {
    if (dobb(bbs->arr + i, 0) != 0) return -1;
    calcetp(bbs->arr + i);
    calcab(bbs->arr + i);
  }
  return 0;
}

/* obsolete, to be replace bb_export */
static void export_old(bb_t *bb, const char *fn)
{
  distr2d_t *d = distr2d_open(-M_PI, M_PI, bb->binw, -M_PI, M_PI, bb->binw);
  int i, j, id, n = bb->bins, m = bb->bins;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      distr2dsum_t *ds = d->arr + i * (m + 1) + j;
      id = i * n + j;
      ds->s = bb->hist[id];
      ds->sf = bb->sbfphi[id];
      ds->sg = bb->sbfpsi[id];
      if (ds->s > 0) {
        ds->sf2 = ds->sf * ds->sf / ds->s;
        ds->sg2 = ds->sg * ds->sg / ds->s;
        ds->sfg = ds->sf * ds->sg / ds->s;
      } else {
        ds->sf2 = ds->sg2 = ds->sfg = 0;
      }
    }
  distr2d_save(d, fn);
  distr2d_close(d);
}

int main(int argc, char *argv[])
{
  bbs_t *bbs;
  doargs(argc, argv);
  if ((bbs = bbinit()) == NULL) return 1;
  if (fnport) { /* export to the portable format */
    export_old(bbs->arr, fnport);
  } else {
    if (dobs(bbs) != 0) return 1;
    bbs_writebin(bbs, fnbout, opver);
    bbs_write(bbs, fntout, opver);
  }
  bbs_close(bbs);
  return 0;
}


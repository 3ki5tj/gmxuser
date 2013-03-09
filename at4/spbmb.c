/* multiple-bin estimator for dihedral angle
 * Copyright (C) 2010 Cheng Zhang */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "md2spb.h"

const char *prog = "spbmb";
const char *fnin = "pre.bin";
const char *fnbout = "spbmb.bin";
const char *fntout = "spbmb.txt";
const char *atcfg = "at.cfg";
int verbose;
int ipver = 0;
int opver = 2;
int use_atcfg = 1;
int multibins = 20;
int dermbin = 0;  /* number of bins to average the force */
int normalize = 1;
int mfcomnz = 1; /* average mean force only to nonempty bins */
int polymer = 0;

/* print help message and die */
static void help(void)
{
  static char options[] =
  "OPTIONS:\n"
  " -m: number of bins to average\n"
  " -M: number of bins to average derivatives\n"
  " -r: not to normalize\n"
  " -e: apply average mean force to empty bins as well\n"
  " -p: convert from/to polymer conversion\n"
  " -i: input file\n"
  " -T: text output\n"
  " -B: binary output\n"
  " -v: verbose\n"
  " -h: print this message\n";
  fprintf(stderr, "%s Copyright (C) 2010  Cheng Zhang\n", prog);
  fprintf(stderr, "USAGE:\n"
      "  %s [OPTIONS] spb.bin\n\n%s", prog, options);
  exit(1);
}

/* handling options and get the PDB file name, fname0  */
static void doargs(int argc, char **argv)
{
  int i, j, ch;
  const char *val;

  prog = argv[0];
  verbose = 0;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnin = argv[i];
      continue;
    }
    ch = argv[i][1];
    if (strchr("fbiTBmM", ch)) {
      if (argv[i][2] == '\0') {
        if (i == argc - 1) {
          fprintf(stderr, "miss argument after -%c\n", ch);
          help();
        }
        val = argv[++i];
      } else {
        val = argv[i] + 2;
      }
      switch (ch) {
      case 'i': case 'b': case 'f': fnin = val; break;
      case 'T': fntout = val; break;
      case 'B': fnbout = val; break;
      case 'm':
        multibins = atoi(val);
        printf("average over %d bins\n", multibins);
        break;
      case 'M':
        dermbin = atoi(val);
        printf("derivative average over %d bins\n", dermbin);
        break;
      default:
        fprintf(stderr, "cannot handle option %s\n", argv[i]);
        exit(1);
      }
      continue;
    }

    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'r':
        normalize = 0;
        printf("distributions are not normalized\n");
        break;
      case 'e':
        mfcomnz = 0;
        printf("remove mfcom for empty bins too\n");
        break;
      case 'p':
        polymer = 1;
        printf("polymer conversion needed\n");
        break;
      case 'v': verbose++; break;
      case 'h': default: help();
      }
  }
}

/* initialize spbonds_t and load data */
static spbonds_t *spbinit(void)
{
  spbonds_t *bs;
  unsigned read_flags = SPB_VERBOSE;

  if ((bs = calloc(1, sizeof(*bs))) == NULL) {
    fprintf(stderr, "cannot allocate space\n");
    return NULL;
  }

  /* load basic information from configure first */
  if (use_atcfg) {
    cfgdata_t *cfg;

    if ((cfg = cfgopen(atcfg)) == NULL) {
      fprintf(stderr, "cannot read config %s.", atcfg);
      return NULL;
    }
    if (spbs_cfgopen_low(bs, cfg, NULL) != 0) {
      return NULL;
    }
    cfgclose(cfg);

    read_flags |= SPB_VERIFY;
    fprintf(stderr, "configuration file %s is used\n", atcfg);
  }

  /* binary file already has mf, use as is */
  read_flags |= SPB_MF;
  if (spbs_readbin(bs, fnin, &ipver, read_flags) != 0) {
    fprintf(stderr, "error during loading spb binary data\n");
    return NULL;
  }
  fprintf(stderr, "%s version %d\n", fnin, ipver);
  spbs_calcmfpmf(bs, 0);
  return bs;
}

#define WRAP(i)  ((i + 2*spb->bins) % spb->bins)

/* recalculate the pmf for the modulating factor */
static int spbpmf(spb_t *spb)
{
  int i, cnt;
  double mfcom;

  /* compute the mean force */
  for (i = 0; i < spb->bins; i++) {
    if (spb->hist[i] > .1) {
      spb->mf[i] = (spb->sbf[i] + spb->sdiv[i])/spb->sbl[i];
    } else spb->mf[i] = 0.;
  }

  /* smoothen over to get bins */
  if (dermbin > 0) {
    double *mf1, sm;
    int iw, j, js;

    if ((mf1 = malloc(sizeof(double)*spb->bins)) == NULL) {
      fprintf(stderr, "average of the mean force\n");
      exit(1);
    }
    if ((iw = dermbin/2) < 0) iw = 1;
    for (i = 0; i < spb->bins; i++) {
      for (sm = 0.0, j = i - iw; j <= i + iw; j++) {
        js = (j + spb->bins) % spb->bins;
        sm += spb->mf[js];
      }
      mf1[i] = sm/(iw*2+1);
      //printf("%d, mf = %g, mf1 = %g\n", i, spb->mf[i], mf1[i]);
      //getchar();
    }
    for (i = 0; i < spb->bins; i++)
      spb->mf[i] = mf1[i];
    free(mf1);
  }

  /* Note we should not attempt to try to fill in the missing bins, from the highest peak */

  /* remove the mfcom */
  cnt = 0;
  mfcom = 0.;
  for (i = 0; i < spb->bins; i++) {
    if (fabs(spb->mf[i]) < 1e-14) {
      cnt++;
    } else {
      mfcom += spb->mf[i];
    }
  }
  printf("%d zeroes, mean force residue = %g, ", cnt, mfcom*spb->binw);
  mfcom /= (mfcomnz ? (spb->bins - cnt) : spb->bins);
  for (i = 0; i < spb->bins; i++) {
    if (fabs(spb->mf[i]) > 1e-14 || !mfcomnz)
      spb->mf[i] -= mfcom;
  }

  /* compute pmf */
  for (spb->pmf[0] = 0., i = 0; i < spb->bins; i++)
    spb->pmf[i+1] = spb->pmf[i] - spb->mf[i]*spb->binw;

  return 0.0;
}

static int pswitch(spb_t *spb)
{
  int i, n = spb->bins/2;
  double tmp;

#define PSWAP(arr, add1) \
  for (i = 0; i < n; i++) { \
    tmp = arr[i], arr[i] = arr[i+n], arr[i+n] = tmp; } \
  if (add1) arr[2*n] = arr[0];

  PSWAP(spb->hist, 0);
  PSWAP(spb->sbf, 0);
  PSWAP(spb->sdiv, 0);
  PSWAP(spb->sbl, 0);
  PSWAP(spb->mf, 0);
  PSWAP(spb->pmf, 1);
  PSWAP(spb->distref, 1);
  PSWAP(spb->distref0, 1);
  return 0;
}

/* handle a single dihedral */
static int dospb(spb_t *spb)
{
  int i, j, jr, j0, j1, iw = multibins/2;
  double bl = 1., sum = 0., fac, bn, num, den0, den, sum1;

  spbpmf(spb);
  /* first nonzero bin to compute bl */
  for (i = 0; i < spb->bins; i++)
    if (spb->hist[i] > .5) {
      bl = spb->sbl[i]/spb->hist[i];
      break;
    }
  for (i = 0; i < spb->bins; i++) {
    sum += spb->hist[i];
  }
  fac = 1.0/(sum * spb->binw);
  printf("spb %d: %g samples, ", spb->id, sum);
  for (i = 0; i <= spb->bins; i++) {
    num = den = den0 = 0.;
    j0 = i - iw;
    j1 = i + iw;
    for (j = j0; j < j1; j++) {
      jr = WRAP(j);
      num += spb->hist[jr];
    }
    for (j = j0; j <= j1; j++) {
      jr = WRAP(j);
      bn = (j == j0 || j == j1) ? .5 : 1.0;
      den0 += bn;
      den += bn * exp(-(spb->pmf[jr] - spb->pmf[i]) * bl);
    }
    spb->distref[i] = num/den; /* unbiased result */
    spb->distref0[i] = num/den0; /* biased result  */
  }
  for (sum1 = 0.0, i = 0; i < spb->bins; i++)
    sum1 += spb->distref[i];
  sum1 *= 1.0/sum;
  printf("sum: %g\n", sum1);

  if (normalize) {
    for (i = 0; i <= spb->bins; i++) {
      if (i < spb->bins) {
        spb->hist[i] *= fac;
        spb->sbf[i] *= fac;
        spb->sdiv[i] *= fac;
        spb->sbl[i] *= fac;
      }
      spb->distref[i] *= fac;
      spb->distref0[i] *= fac;
    }
  }
  if (polymer) pswitch(spb);
  return 0;
}

static int dobs(spbonds_t *bs)
{
  int i;

  for (i = 0; i < bs->cnt; i++) {
    if (dospb(bs->arr + i) != 0) return -1;
  }
  return 0;
}

int main(int argc, char *argv[])
{
  spbonds_t *bs;

  doargs(argc, argv);
  if ((bs = spbinit()) == NULL) return 1;
  if (dobs(bs) != 0) return 1;
  spbs_writebin(bs, fnbout, opver);
  spbs_write(bs, fntout, opver, /* as is */1);
  spbs_close(bs);
  return 0;
}


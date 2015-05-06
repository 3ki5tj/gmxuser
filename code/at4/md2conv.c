/*
  converts a spb.bin/spb.cfg file to another

  Copyright (C) 2010  Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define ZCOM_PICK
#define ZCOM_CFG
#define ZCOM_SS
#include "zcom2.h"

#include "md2spb.h"
#include "md2spbx.h" /* additional routines */

const char *prog = "spbconv";
const char *fnin = "pre.bin";
const char *fntxtout = "spb.txt";
const char *fnbinout = "spb.bin";
const char *fnfitlog = "spbfit.log";
const char *atcfg = "at.cfg";
int verbose;
int ipver = 0;
int opver = 2; /* output version */
int use_atcfg = 1; /* search for at.cfg, 1 is more convenient */
int binput = 1; /* binary input */
int recalc_pmf = 0;
int fitting = 1;
int add_corr = 0; /* add a correction from the difference between
                     histogram and target distribution. */
int barrier_margin = 5; /* number of bins in smoothing potential for reducing barriers */
int iter_max = 10000;
/* stop the minimizer if improvement after a round becomes too small */
double impr_min = 1e-6;
int polymer = 0;

/* print error message and die */
static void help(void)
{
  static char options[]=
  "OPTIONS:\n"
  " -a: followed by a configuration file (at.cfg) to supplement additional information.\n"
  " -b: set binary-input mode, followed by a binary spb file (spb.bin)\n"
  " -c: add a correction if histogram deviates the reference distribution (recommended)\n"
  " -I: followed by the maximal rounds of iteration\n"
  " -t: set text-input mode, followed by a text spb file (spb.txt)\n"
  " -B: followed by a binary output file (default spb_out.bin)\n"
  " -T: followed by a text output file, (default spb_out.txt)\n"
  " -m: followed by tolerance of a minimizer\n"
  " -g: followed by file for initial guess\n"
  " -d: followed by # of bins for polishing barrier edges\n"
  " -V: followed by output version\n"
  " -p: switch to/from polymer convention\n"
  " -j: no decomposing, just do read/write\n"
  " -h: print this message\n"
  " -v: verbose, -vv to be more\n"
  "\n";
  fprintf(stderr, "%s  Copyright (C) 2010  Cheng Zhang\n"
  "This program comes with ABSOLUTELY NO WARRANTY. "
  "It is free software, and you are welcome to redistribute "
  "it under certain conditions.\n\n", prog);

  fprintf(stderr, "USAGE:\n");
  fprintf(stderr, "%s [OPTIONS] spb.bin\n\n", prog);
  fprintf(stderr,
          "Example:\n"
          "Using version 0 text file spb.txt with at.cfg:\n"
          "  %s -t spb.txt -a at.cfg\n\n"
          "Using binary file spb.bin\n"
          "  %s -b spb.bin\n\n",
          prog, prog);
  fprintf(stderr, "%s", options);
  fflush(stderr);
  exit(1);
}

/* handling options and get the PDB file name, fname0  */
static int doargs(int argc, char **argv)
{
  int i, j, itmp, ch;
  const char *val;
  double x;

  prog = argv[0];
  verbose = 0;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnin = argv[i];
      continue;
    }
    ch = argv[i][1];
    if (strchr("btaBTgIimVd", ch)) {
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
      case 'b': case 't':
        fnin = val;
        binput = (ch == 'b');
        fprintf(stderr, "Using %s mode input file %s\n",
            binput ? "binary" : "text", fnin);
        break;
      case 'B': case 'T': case 'a': case 'g':
        if (ch == 'a') {
          atcfg = val;
        } else if (ch == 'B') {
          fnbinout = val;
        } else if (ch == 'T') {
          fntxtout = val;
        } else if (ch == 'g') {
          fnfitlog = val;
        }
        fprintf(stderr, "The file for switch -%c is set to %s\n", ch, val);
        break;
      case 'I': case 'i':
        if ((itmp = atoi(val)) > 0) {
          iter_max = itmp;
          printf("maximal # of iteration rounds is set to %d\n", iter_max);
        }
        break;
      case 'm':
        if ((x = atof(val)) > 0 && x < 1.0) {
          impr_min = x;
          printf("minimizer tolerance has been set to %g\n", impr_min);
        }
        break;
      case 'd':
        if ((itmp = atoi(val)) >= 0) {
          barrier_margin = itmp;
          printf("barrier margin is %d\n", barrier_margin);
        }
        break;
      case 'V':
        if ((itmp = atoi(val)) >= 0) {
          opver = itmp;
          printf("output version is set to %d\n", opver);
        }
        break;
      default:
        fprintf(stderr, "cannot handle option %s\n", argv[i]);
        exit(1);
      }
      continue;
    }

    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'j': fitting = 0; break;
      case 'c': add_corr = 1; break;
      case 'p': polymer = 1; break;
      case 'r': recalc_pmf = 1; break;
      case 'v': verbose++; break;
      case 'h': default: help();
      }
  }
  return 0;
}

/* initialized spbonds_t and load data */
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
  } else {
    /* disable verification if no cfg initialization */
    fprintf(stderr, "popscal information needs to be supplemented separately\n");
    exit(1);
  }

  if (binput) {
    /* binary file already has mf, use as is */
    read_flags |= SPB_MF;
    if (spbs_readbin(bs, fnin, &ipver, read_flags) != 0) {
      fprintf(stderr, "error during loading spb binary data\n");
      return NULL;
    }
  } else {
    /* text file already has pmf, but we just recalculate everything */
    /* read_flags |= SPB_NOPMF;*/
    if (spbs_read(bs, fnin, read_flags) != 0) {
      fprintf(stderr, "error during loading spb txt data\n");
      return NULL;
    }
  }

  if (recalc_pmf) {
    fprintf(stderr, "recalculating PMF...\n");
    spbs_calcmfpmf(bs, 0);
  }
  return bs;
}

#define FIT_ORDER 4
#define FART     (1)            /* first parameter for the left side */
#define FALF      (FIT_ORDER)    /* first parameter for the right side */
#define FALAST    (2*FIT_ORDER)  /* the last parameter of serial expansion */
#define ID_XC     (FALAST+1)     /* the center of the peak */
#define ID_A      (FALAST+2)     /* the height of the peak */
#define FTOTCNT   (FALAST+3)     /* # of total parameter */

#define PARTMAX (SPB_PMFSCAL_CNT-1) /* different parts, alpha, beta, turn */
#define FITARR(i)  (f->arr[(i+f->period)%f->period])  /* safer than f->arr[i] */
double fit_penover = 5.0; /* penalty factor for being over than being under */

typedef struct {
  int id, type;
  int i0, i1;
  int period;
  double xpeak0; /* a rough estimate of a[ID_XC] */
  double penover;
  double xmin, xmax, xperiod, dx; /* the whole range */
  double x0, x1; /* fitting range */
  double x;  /* position of the peak */
  double a[FTOTCNT]; /* fitting parameters, the last two xpeak and A */
  double fa[FTOTCNT]; /* the corresponding force */
  double baseline;
  double bl;
  double *arr;
} fit_t;

/* since it's a potential, that's the highest point */
static double fit_calcbaseline(double arr[], int cnt, int doshift)
{
  int i;
  double y, max = -1e9, min = 1e9, baseline;
  for (i = 0; i < cnt; i++) {
    y = arr[i];
    if (y > max) {
      max = y;
    } else if (y < min) {
      min = y;
    }
  }
  baseline = max+(max-min)*1e-8; /* add a little margin */

  /* do a shifting */
  if (doshift)
    for (i = 0; i < cnt; i++)
      arr[i] -= baseline;
  return baseline;
}


/* estimate peak symmetrical parameters */
static int fit_est_peak(fit_t *f)
{
  int i, ic;
  double y, min = 1e9;

  /* locate the minimal value -> center of the peak (negative) */
  for (ic = i = f->i0; i < f->i1; i++) {
    y = FITARR(i);
    if (y < min) {
      min = y;
      ic = i;
    }
  }

  /* TODO: refine the value A and xpeak */
  f->a[ID_A]  = -FITARR(ic);
  f->xpeak0   =  f->xmin+ic*f->dx;
  f->a[ID_XC] =  f->xpeak0;

  /* make sure ic is deeper than its neighbors */
  if (FITARR(ic) > FITARR(ic-1) ||
      FITARR(ic) > FITARR(ic+1)) {
    printf("ERROR!!! peak center=%d is not a peak i0=%d, i1=%d\n", ic, f->i0, f->i1);
    return 1;
  }

  /* moderate the range */
  for (i = ic-1; i >= f->i0; i--) {
    y = FITARR(i);
    if (y < FITARR(i+1) ) {
      fprintf(stderr, "left rebound occurs, truncate at i=%d\n", i%f->period);
      f->i0 = i+1;
      break;
    }
  }
  for (i = ic; i < f->i1; i++) {
    y = FITARR(i);
    if (y > FITARR(i+1) ) {
      fprintf(stderr, "right rebound occurs, truncate at i=%d\n", i%f->period);
      f->i1 = i+1;
      break;
    }
  }

  /* set all series expansion coefficients zero */
  for (i = 0; i <= FALAST; i++)
    f->a[i] = 0.0;

  /* obtain an first guess of a2, from drop to 80% */
  for (i = ic - 1; i >= f->i0; i--)
    if (FITARR(i) > 0.8 * FITARR(ic))
      break;
  y = (i - ic) * f->dx;
  f->a[FALF] = 1.0/(y*y); /* set the first expansion coefficient (right-hand-side) */

  for (i = ic + 1; i < f->i1; i++)
    if (FITARR(i) > 0.8 * FITARR(ic))
      break;
  y = (i - ic) * f->dx;
  f->a[FART] = 1.0/(y*y); /* set the first expansion coefficient (left-hand-side) */


  printf("PEAK at x=%g(%d), height is=%g, baseline at=%g, a2=%g (%d, %d)\n",
      f->a[ID_XC], ic, f->a[ID_A], f->baseline, f->a[FART], f->i0, f->i1);
/*
  { // write a file
    char fname[128];
    FILE *fp;
    sprintf(fname, "arr%d%c", f->type, 'A'+f->id);
    fp=fopen(fname, "w");
    if(fp){
      //for(i=f->i0; i<f->i1; i++){
      for(i=0; i<f->period; i++){
        y=FITARR(i);
        fprintf(fp, "%g %g\n", f->xmin+i*f->dx, y);
      }
      fclose(fp);
    }
    getchar();
  }
*/
  return 0;
}

/* evaluate the value from a single point x = f->min + i*f->dx */
static double fit_eval(fit_t *f, int i, double *pdx, double *pyt0, double *pder, int eonly)
{
  double xspan, x0, xpeak, dxp, dx, dx2, yt, xn2, xn, yt0, der;
  int od, doder = (!eonly), odmin, odmax, is_left;

  /* compute the distance from the center of the peak */
  x0 = i*f->dx + f->xmin;
  xspan = f->xperiod;
  xpeak = f->a[ID_XC];
  /* shift xpeak inside the x-range */
  while (xpeak < f->xmin)
    xpeak += xspan;
  while (xpeak > f->xmax)
    xpeak -= xspan;
  /* now both x0 and xpeak are in the range, dx < f->xperiod */
  dx = x0 - xpeak;

  if (dx > 0) { /* x0 is on the right hand side of the peak */
    is_left = 1;
  } else {
    is_left = 0;
    dx = -dx;
  }
  dxp = xspan - dx; /* the period image of x0 distances from xpeak dx */
  if (dxp < dx) { /* flip side */
    dx = dxp;
    is_left = !is_left;
  }
  /* determin which set of series-expansion coefficients to use */
  odmin = (is_left) ? FART : FALF;
  odmax = odmin + FIT_ORDER - 1;
  if (!is_left)
    dx = -dx;
  if (pdx)
    *pdx = dx;
  dx2 = dx * dx;

  for (yt0 = 0.0, der = 0.0, xn = dx, xn2 = dx2, od = odmin; od < odmax; od++) {
    yt0 += xn2 * f->a[od];
    if (doder) {
      der += xn * f->a[od] * (od*2.0);
      xn *= dx2;
    }
    xn2 *= dx2;
  }
  yt0 = exp(-yt0);
  if (pyt0)
    *pyt0 = yt0;
  if (doder && pder)
    *pder = der;
  yt = -f->a[ID_A] * yt0;
  return yt;
}

#define fit_calcerr(f) fit3_calcerr(f, 1)
#define fit3_calcerr(farr,cnt) fit3_calcerr_x(farr,cnt,0,0)

/* the residue error, or the energy function */
static double fit3_calcerr_x(fit_t farr[], int cnt, int verb, int eonly)
{
  int i, od, id, i0, i1;
  double res, y, yt, yt0[PARTMAX], dy, dy2, dx[PARTMAX], dx2, penover, der[PARTMAX], tmp;
  double fa[PARTMAX][FTOTCNT];
  fit_t *f;

  if (cnt > PARTMAX ) {
    fprintf(stderr, "cannot do so many %d fas\n", cnt);
    return 1e9;
  }
  for (id = 0; id < cnt; id++) {
    /* clear all gradients for all fitting parameters */
    for (od = 0; od < FTOTCNT; od++)
      fa[id][od] = 0.0;
  }

  i0 = farr[0].i0;
  i1 = farr[0].i1;
  for (res = 0.0, i = i0; i < i1; i++) {
    f = farr;
    y = FITARR(i); /* the target value */
    /* adding components from different fitting functions */
    for (yt = 0.0, id = 0; id < cnt; id++) {
      f = farr + id;
      yt += fit_eval(f, i, &dx[id], &yt0[id], &der[id], eonly);
    }
    dy = yt - y;
    dy2 = dy * dy;
    penover = (dy < 0) ? farr->penover : 1.0;
    res += penover * dy2;

    if (verb >= 1) {
      for (id = 0; id < cnt; id++) {
        f = farr + id;
        fprintf(stderr, "CALCERR: type %d, id %c, i %4d: real=%16.8f model=%16.8f, dy=%g\n",
          f->type, f->id+'A', i, y, yt, dy);
      }
    }

    if (eonly) continue;

    /* computing derivatives */
    for (id = 0; id < cnt; id++) {
      int odmin, odmax;

      tmp = 2.0 * penover * dy * yt;
      dx2 = dx[id] * dx[id];
      odmin = (dx[id] > 0) ? FART : FALF;
      odmax = odmin + FIT_ORDER - 1;
      /* we can only set one side of series-expansion parameters */
      for (od = odmin; od < odmax; od++) {
        fa[id][od] += tmp;
        tmp *= dx2;
      }
      fa[id][ID_XC] += -der[id] * 2.0 * penover * dy * yt * dx[id]; /* the xpeak shift force*/
      fa[id][ID_A] += 2.0 * penover * dy * yt0[id]; /* the amplitude force */
    }
  }

  if (!eonly) {
    for (id = 0; id < cnt; id++) { /* components */
      f = farr+id;
      /* multiple all series-expansion coefficents by the interval dx */
      for (od = 0; od <= FALAST; od++)
        f->fa[od] = fa[id][od]*farr[0].dx;
    }
  }
  return res*farr->dx;
}

/* a stupid way of doing linear minimizing,
 * force is not used here */
static double fit3_linmin(fit_t *farr, int cnt, int id, int i,
    double step0, int *pmod, int *peval, int verb)
{
#define HVAL 1e36
  const int t1 = 100;
  const double shrink = 0.08;
  const double x_tol = 1e-14;

  double yp, step, impr = 1.0;
  int t, mod, eval;

  struct tag_searchpt{
    double x;
    double y;
    int valid;
    int hit;
  } now = { 0.0, 0.0, 0, 0 },
    plus = { 0.0, 0.0, 0, 0 }, minus = { 0.0, 0.0, 0, 0 },
    min = {-HVAL, HVAL, 0, -1 }, max = { HVAL, HVAL, 0, 1 };

  /* guess the initial step */
  if (step0 > 0.0) step = step0; else step = 0.2;

  mod = 0; eval = 0;

/* NOTE: we do not use gradient in this routine! So no need to back up force */
/* evaluate and validate a point */
#define lm_eval(pt) { \
  farr[id].a[i]=pt.x; \
  pt.y=fit3_calcerr_x(farr,cnt,verb-3,1); \
  pt.valid=1; \
  eval++; }

/* validated-copy of a point
 * the `src' is first evaluated before copying if it is invalid;
 * `dest' is then copied, and is always valid
 * the hit flag is also assigned
 * */
#define lm_vcopy(dest, src, ht){ \
  if( !(dest.valid=src.valid) ){ lm_eval(src); } \
  dest.x=src.x; \
  dest.y=src.y; \
  dest.valid=1; \
  dest.hit=ht; }

  now.x = farr[id].a[i];
  lm_eval(now);
  yp = now.y;
  if (verb >= 3) {
    fprintf(stderr, "[INIT  ]        id=%3d, i=%3d: a=%14.10f res=%g, step=%g\n",
        id, i, now.x, now.y, step);
  }

  if (i != ID_XC) min.x = 0.0; /* this is true even for A  */

  for (t = 0; t < t1 && step > x_tol && impr > 1e-14; t++) {
    if (max.x-min.x < 2*x_tol) break;

    if (now.hit < 0) {
      lm_vcopy(minus, now, -1);
    } else {
      minus.x = now.x-step;
      minus.hit = 0;
      if (minus.x < min.x+x_tol) {
        lm_vcopy(minus, min, -1);
      } else {
        lm_eval(minus);
      }
    }

    /* if minus is good enough, no need to compute plus */
    if (minus.y < now.y) { plus.y = now.y+1; goto CHOOSE_MINUS; }

    if (now.hit > 0) {
      lm_vcopy(plus, now, 1);
    } else {
      plus.x = now.x+step;
      plus.hit = 0;
      if (plus.x > max.x-x_tol) {
        lm_vcopy(plus, max, 1);
      } else {
        lm_eval(plus);
      }
    }

    if (plus.y >= now.y && minus.y >= now.y) {
      step *= shrink;
      /* update the bracket */
      lm_vcopy(min, minus, -1); /* minus is the new min */
      lm_vcopy(max, plus, 1);  /* plus is the new max */
      if (verb >= 3) {
        printf("[FAILED] t=%3d, id=%3d, i=%3d: (%10.6f,%10.6f) a=%10.6f,step=%12.6e,"
            "hit=%+2d,%+2d,eval=%4d res=%8.3f,incr=%g,%g\n",
            t, farr[id].id, i, min.x, max.x, now.x, step,
            minus.hit, plus.hit, eval,
            now.y, plus.y-now.y, minus.y-now.y);
      }
    } else {
      int ch;
      double prevy;
CHOOSE_MINUS:
      prevy = now.y;
      if (plus.y < minus.y) {
        lm_vcopy(min, now, -1); /* now is the new min */
        lm_vcopy(now, plus, 0); /* plus is the new now */
        ch='+';
      } else {
        lm_vcopy(max, now, 1); /* now is the new max */
        lm_vcopy(now, minus, 0); /* minus is the new now */
        ch='-';
      }
      mod = 1;
      impr = prevy-now.y;
      if (verb >= 2) {
        printf("[MOVE %c] t=%3d, id=%3d, i=%3d: (%10.6f,%10.6f) a=%10.6f,step=%12.6e,"
            "hit=%+2d,%+2d,eval=%4d res=%8.3f,improved=%g\n",
            ch, t, farr[id].id, i, min.x, max.x, now.x, step,
            minus.hit, plus.hit, eval, now.y, impr);
      }
    }
    farr[id].a[i] = now.x;
    if (verb > 5) getchar();
  }
  if (pmod) *pmod = mod;
  if (peval) *peval += eval;
  if (verb >= 2 || (verb >= 1&&mod)) {
    printf("t=%d,step=%g, err=%.8f, improved=%g\n", t, step, now.y, yp-now.y);
  }
  return now.y;
}

/* save fitting data to string */
static int fit3_save(const fit_t *farr, int parts, int id, char *s, int bufsiz)
{
  const int margin = 30;
  int i, j, next, nwr = 0;
  char *p = s;

  sprintf(p, "%d %n", id, &next);
  p += next; if ((nwr += next) >= bufsiz-margin) return -1;
  for (i = 0; i < parts; i++) {
    sprintf(p, "%d %n", i, &next);
    p += next; if ((nwr += next) >= bufsiz-margin) return -1;
    for (j = 0; j < FTOTCNT; j++) {
      sprintf(p, "%20.15f %n", farr[i].a[j], &next);
      p += next; if ((nwr += next) >= bufsiz-margin) return -1;
    }
  }
  return 0;
}

/* load data from string */
static int fit3_load(fit_t *farr, int parts, int id, const char *s)
{
  int i = 0, j, itmp, next;
  const char *p;

  p = s;
  if (1 != sscanf(p, "%d%n", &itmp, &next) || itmp != id) {
    fprintf(stderr, "id = %d vs. %d\n", id, itmp);
    return -1;
  }
  p += next;
  for (i = 0; i < parts; i++) {
    if (1 != sscanf(p, "%d%n", &itmp, &next) || itmp != i) {
      fprintf(stderr, "i = %d vs. %d\n", itmp, i);
      return -1;
    }
    p += next;
    for (j = 0; j < FTOTCNT; j++) {
      if (1 != sscanf(p, "%lf%n", &farr[i].a[j], &next)) {
        fprintf(stderr, "i = %d, data stopped at j = %d\n", i, j);
        return -1;
      }
      p += next;
    }
  }
  return 0;
}

#define LCMAX 4096
static int readlines(char lines[][LCMAX], int cnt, const char *fname)
{
  FILE *fp;
  int i;

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot open %s.\n", fname);
    return -1;
  }
  for (i = 0; i < cnt; i++) {
    if (fgets(lines[i], LCMAX, fp) == NULL) {
      break;
    }
  }
  fprintf(stderr, "load %d lines from %s.\n", i, fname);
  fclose(fp);
  return i;
}

/* print a GNUPLOT friendly output */
static int fit3_print(fit_t *farr, int cnt, double baseline)
{
  char *output;
  char buf[256];
  int id, od, odmin, odmax, chid, do_minus, type;
  fit_t *f;
  char *fname = NULL, *p;
  double xc;
  int inc = (ipver >= 2 ? 8 : 7);

  do_minus = 1;

  if ((output = ssnew(1024)) == NULL)
    return -1;
  sprintf(buf, "dx(x, xc) = x - xc\n");
  sscat(output, buf);
  sprintf(buf, "dx2(x, xc) = dx(x, xc)**2\n");
  sscat(output, buf);

  for (id = 0; id < cnt; id++) {
    f = farr+id;
    chid = id+'A';

    sprintf(buf, "rt%d%c(x) = -%g * exp(-",
        f->type, chid, f->a[ID_A]);
    sscat(output, buf);
    /* print the left hand side of the function */
    odmin = FART;
    odmax = odmin + FIT_ORDER - 1;
    for (od = odmin; od < odmax; od++) {
      if (od < odmax - 1) {
        sprintf(buf, "x*(%.8f+", f->a[od]);
      } else {
        sprintf(buf, "x*%.8f", f->a[od]);
      }
      //printf("BP1: od = %d, buf=%p, %s\n", od, buf, buf); fflush(stdout);
      sscat(output, buf);
    }
    for (od = 1; od < FIT_ORDER; od++)
      sscat(output, ")");
    sscat(output, "\n");

    sprintf(buf, "lft%d%c(x) = -%g * exp(-",
        f->type, chid, f->a[ID_A]);
    sscat(output, buf);
    /* print the right hand side of the function */
    odmin = FALF;
    odmax = odmin + FIT_ORDER - 1;
    for (od = odmin; od < odmax; od++) {
      if (od < odmax - 1) {
        sprintf(buf, "x*(%.8f+", f->a[od]);
      } else {
        sprintf(buf, "x*%.8f", f->a[od]);
      }
      //printf("BP2: od=%d, buf=%p, %s\n", od, buf, buf); fflush(stdout);
      sscat(output, buf);
    }
    for (od = 1; od < FIT_ORDER; od++)
      sscat(output, ")");
    sscat(output, "\n");

    xc = f->a[ID_XC];
    if (polymer) xc = (xc > 0) ? (xc - M_PI) : (xc + M_PI);
    if (xc < 0) xc += f->xperiod;
    sprintf(buf, "f%d%c(x) = (dx(x,%g)>0) ? rt%d%c(dx2(x,%g)) : lft%d%c(dx2(x,%g))\n",
      f->type, chid, xc,
      f->type, chid, xc,
      f->type, chid, xc);
    sscat(output, buf);
    if (do_minus) {
      xc -= f->xperiod;
      sprintf(buf, "f%d%cb(x) = (dx(x,%g)>0) ? rt%d%c(dx2(x,%g)) : lft%d%c(dx2(x,%g))\n",
        f->type, chid, xc,
        f->type, chid, xc,
        f->type, chid, xc);
      sscat(output, buf);
    }
  }

  type = farr->type;
  sprintf(buf,"f%d(x)=f%dA(x)+f%dAb(x)+f%dB(x)+f%dBb(x)+f%dC(x)+f%dCb(x)\n",
      type, type, type, type, type, type, type);
  sscat(output, buf);
  /* guess the best file name */
  sscpy(fname, fnin);
  if (binput) { /* try to transform it to the corresponding text file */
    if ((p = strrchr(fname, '.')) != NULL)
      strcpy(p, ".txt");
    else
      sscat(fname, ".txt");
  }
  sprintf(buf, "plot [-pi:pi]\"%s\" u 1:(-$%d-log($%d)/(%.8f)%+.8f) w l, "
      "f%dA(x)+f%dAb(x), f%dB(x)+f%dBb(x), f%dC(x)+f%dCb(x), f%d(x), 0\n",
      fname, 6+farr->type*inc, 7+farr->type*inc, farr->bl, -baseline,
      farr->type, farr->type,
      farr->type, farr->type,
      farr->type, farr->type, farr->type);
  sscat(output, buf);

  printf("%s", output);
  ssdel(output);
  return 0;
}

#define fit_remove_component(arr, f) fit3_remove_component(arr, arr, f, 1)
static int fit3_remove_component(double dest[], double src[], fit_t *farr, int cnt)
{
  int i, id;
  double yt;

  for (i = 0; i <= farr->period; i++) {
    yt = 0;
    for (id = 0; id < cnt; id++) {
      yt += fit_eval(farr+id, i, NULL, NULL, NULL, 1);
    }
    dest[i] = src[i]-yt;
  }
  return 0;
}

#define fit_minimize(f,verb) fit3_minimize(f,1,1,verb)
static double fit3_minimize(fit_t *farr, int cnt, int restrained, int verb)
{
  double res0, res, resr0, impr, dx;
  int it = 0, t = 0, id, ia, iamax, mod = 0, eval = 0, tmp;

  res = res0 = fit3_calcerr(farr, cnt); eval++;

  for (it = 0; it < iter_max; it++) {
    mod = 0;
    resr0 = res;
    for (id = 0; id < cnt; id++) {
      iamax = (restrained) ? FALAST : (FTOTCNT-1);
      for (ia = 1; ia <= iamax; ia++) { /* fit every parameter */
        res = fit3_linmin(farr, cnt, id, ia, 0.0, &tmp, &eval, verb-1);

        if (ia == ID_XC) {
          dx = farr[id].xmax - farr[id].xmin;
          while (farr[id].a[ia] < farr[id].xmin) farr[id].a[ia] += dx;
          while (farr[id].a[ia] > farr[id].xmax) farr[id].a[ia] -= dx;
        }

        if (tmp)
          mod++;
        t++;
      }
    }
    impr = resr0-res;
    if (verb)
      fprintf(stderr, "iter %6d: %4d modifs, %6d evals, err=%12.8f improved=%20.8e,%12.8f ",
        it, mod, eval, res, resr0-res, res0-res);
    if (verb == 1)
      fprintf(stderr, "\r");
    else if (verb > 1)
      fprintf(stderr, "\n");
    if (mod == 0 || impr < impr_min)
      break;
  }

  res = fit3_calcerr(farr, cnt); eval++;

  printf("MINIMIZE: %d iter., %d linmin, %d evals (%6.3f per linmin), err=%.8f, improved=%14.6e; last mod=%4d\n",
      it, t, eval, eval/(t+1e-8), res, res0-res, mod);
  return res;
}

/* try fitting */
static int fit_peak(fit_t *f, double x0, double x1, int last)
{
  if (f->period < 10) {
    fprintf(stderr, "too few points to perform the fitting\n");
    return 1;
  }

  /* if this is not the last component, we use penover to bias against over-fitting */
  f->penover = last ? 1.0 : fit_penover;

  /* adjust the trial range */
  f->x0 = x0;
  f->x1 = x1;
  f->i0 = (int)( (f->x0-f->xmin)/f->dx );
  f->i1 = (int)( (f->x1-f->xmin)/f->dx )+1;
  if (f->i1 < f->i0) f->i1 += f->period; /* fix periodicity */

  if (f->i0 < 0 || f->i1 < 0) {
    fprintf(stderr, "f array size error (%d, %d) \n", f->i0, f->i1);
    return 1;
  }

  fprintf(stderr, "Start fitting part %d/type %d: (%d,%d)\n", f->id, f->type, f->i0, f->i1);

  /* get a very rough guess */
  if (0 != fit_est_peak(f)) return 1;

  /* refine the guess */
  fit_minimize(f, verbose);

  if (!last) fit_remove_component(f->arr, f);

  return 0;
}

static double fit3_geterr(fit_t farr[], int cnt, double ref[])
{
  int id;

  for (id = 0; id < cnt; id++) {
    /* we use the whole range and eliminate asymmetry */
    farr[id].i0 = 0;
    farr[id].i1 = farr[id].period;
    farr[id].penover = 1.0;
    farr[id].arr = ref;
  }
  return fit3_calcerr(farr, cnt);
}

/* refinement, unconstrained global improvement for fitting */
static int fit3_refine(fit_t farr[], int cnt, double ref[], double remain[])
{
  int id;
  double res;

  for (id = 0; id < cnt; id++) {
    /* we use the whole range and eliminate asymmetry */
    farr[id].i0 = 0;
    farr[id].i1 = farr[id].period;
    farr[id].penover = 1.0;
    farr[id].arr = ref;
  }

  res = fit3_calcerr(farr, cnt);
  fprintf(stderr, "Before optimization: res: %g, %g\n", res, sqrt(res/farr[0].xperiod));

  res = fit3_minimize(farr, cnt, 0/*unrestrained*/, verbose);
  fprintf(stderr, "After minimization: res: %g, %g\n", res, sqrt(res/farr[0].xperiod));

  /* save the remaining component to remain */
  fit3_remove_component(remain, ref, farr, cnt);
  return 0;
}

/* reduce barriers */
static int fit3_barred(int n, double y[], int cut, double max)
{
  const double eps = 1e-12;
  double x, a, yc, xdy, dx, dy;
  int i0, i1, i, il, ir, ic;

  for(;;) {
    /* I. look for a barrier */
    /* I.a. look for a left end to start with */
    for (i0 = 0; i0 < n; i0++)
      if (y[i0] < max - eps) break;
    if (i0 == n) {
      fprintf(stderr, "everything is annihilated\n");
      return -1;
    }
    i1 = i0 + n;

    /* I.b. look for the left end */
    for (; i0 < i1; i0++) {
      if (y[i0 % n] > max + eps) break;
    }
    /* no more barriers */
    if (i0 == i1) break;

    /* I.c. look for the right end */
    for (i = i0; i <= i1; i++) {
      if (y[i % n] < max - eps) break;
    }
    if (i > i1) {
      fprintf(stderr, "cannot find the right end!\n");
      return -1;
    }
    i1 = i;
    fprintf(stderr, "barrier find %d - %d\n", i0, i1);

    /* II. reduce the barrier */
    /* II.a. truncate */
    for (i = i0; i < i1; i++)
      y[i % n] = max;

    /* II.b. smoothen the left corner */
    ic = (i0 + i1) / 2; /* center of the barrier */
    il = i0 - cut;
    if (il <= 0) {
      il += n;
      ic += n;
    }
    if (y[il % n] > max) {
      fprintf(stderr, "bad corner %d -- %d\n", il, ic);
      return -1;
    }
    dx = ic - il;
    yc = max - y[il % n];
    xdy = y[il % n] - y[(il-1) % n];
    if (xdy < 0.) {
      fprintf(stderr, "bad corner derivative %g at il = %d\n", xdy, il);
      return -1;
    }
    a = dx*xdy/yc;
    printf("L: ic = %d (%d, %d), il = %d dy/dx = %g, yc/dx = %g, dx = %g, yc = %g, a = %g\n",
        ic, i0, i1, il, xdy, yc/dx, dx, yc, a);
    xdy *= dx;
    for (i = il; i < ic; i++) {
      x = fabs((i - ic)/dx);
      dy = yc * exp(log(x)*a);
      y[i % n] = max - dy;
    }

    /* II.c. smoothen the right corner */
    ic = (i0 + i1)/2;
    ir = i1 + cut;
    if (y[ir % n] > max) {
      fprintf(stderr, "bad corner %d -- %d\n", ic, ir);
      return -1;
    }
    dx = ir - ic;
    yc = max - y[ir % n];
    xdy = y[ir % n] - y[(ir+1) % n];
    if (xdy < 0.) {
      fprintf(stderr, "bad corner derivative %g at ir = %d\n", xdy, ir);
      return -1;
    }
    a = dx*xdy/yc;
    printf("R: ic = %d (%d, %d), ir = %d dy/dx = %g, yc/dx = %g, dx = %g, yc = %g, a = %g\n",
        ic, i0, i1, ir, xdy, yc/dx, dx, yc, a);
    xdy *= dx;
    for (i = ic + 1; i <= ir; i++) {
      x = fabs((i - ic)/dx);
      dy = yc * exp(log(x)*a);
      y[i % n] = max - dy;
    }
  }
  return 0;
}

/* scale the pmf
 * output/out0 are negative logarithm of the new reference distribution
 * out0 is the desired one,
 * output is the one to foo the program */
static int fit3_scalepmf(fit_t farr[], int cnt, double output[], double out0[],
    double remain[], double popscal[], double corr[], double barred)
{
  int i, id, bins;
  double ytnew, ytold, A, An, bl, max, max0, pmftrunc;
  double scal[SPB_PMFSCAL_CNT];

  if (cnt != SPB_PMFSCAL_CNT-1) {
    fprintf(stderr, "input the wrong # %d of parts\n", cnt);
    return -1;
  }

  bins = farr[0].period;

  /* first translate the population scale to
   * the scaling of the pmf itself */
  bl = farr[0].bl;
  for (id = 0; id < cnt; id++) { /* loop over peaks */
    A = farr[id].a[ID_A];
    if (A < 0.0) {
      fprintf(stderr, "Fatal popscal %g for type %d, part %c peak height %g is negative.\n",
          popscal[id], farr[id].type, farr[id].id+'A', bl*A);
      exit(1);
    }
    if (popscal[id] < 1e-30)
      popscal[id] = 1e-30;
    /* An is the new peak after applying the population scaling popscal[id] */
    An = A + log(popscal[id])/bl;
    if (An < 0.0) {
      fprintf(stderr, "Warning popscal %g for type %d, part %c annihilates the peak! (height is %g, minimal %g)\n",
          popscal[id], farr[id].type, farr[id].id+'A', bl*A, exp(-bl*A));
      An = 0.0;
    }
    /* this is how much we scale the potential */
    scal[id] = An / A;
  }
  /* the last component has special meaning */
  scal[cnt] = popscal[cnt];

  /* we assume it is periodic, output[bins] is corrected later */
  for (max0 = -1e9, i = 0; i < bins; i++) {
    for (ytnew = 0.0, ytold = 0.0, id = 0; id < cnt; id++) { /* loop over peaks */
      double tmp;
      tmp = fit_eval(farr+id, i, NULL, NULL, NULL, 1);
      ytnew += tmp * scal[id];
      ytold += tmp;
    }
    ytnew += scal[cnt]*remain[i];
    ytold += remain[i];
    out0[i] = ytnew;
    if (out0[i] > max0)
      max0 = out0[i];
  }

  /* apply a trimming factor to reduce barriers */
  pmftrunc = log(barred)/bl;
  fprintf(stderr, "truncate pmf by %g, barred = %g\n", pmftrunc, barred);
  fit3_barred(bins, out0, barrier_margin, max0 - pmftrunc);

  /* apply the correction from residue distribution,
   * which is used to foo the program */
  for (max = -1e9, i = 0; i < bins; i++) {
    output[i] = out0[i] + (corr ? corr[i]/bl : 0.0);
    if (output[i] > max)
      max = output[i];
  }

  for (i = 0; i < bins; i++) {
    out0[i] = (out0[i] - max0) * bl;
    output[i] = (output[i] - max) * bl;
  }
  out0[bins] = out0[0];
  output[bins] = output[0];

  return 0;
}

/* correct pmf according to the deviation of histogram from distribution */
static double *fit_applycorr(double npmf[], int n, double hist[], double distref[],
    double bl, int periodic)
{
  double *corr, min, maxcorr;
  int i;

  if ((corr = calloc(n+1, sizeof(double))) == NULL) exit(1);
  for (i = 0; i <= n; i++) {
    /* average histogram to approximate distribution */
    if (i == 0) {
      if (periodic) {
        corr[i] = .5 * (hist[0] + hist[n-1]);
      } else {
        corr[i] = .5 * hist[i];
      }
    } else if (i == n) {
      if (periodic) {
        corr[i] = .5 * (hist[0] + hist[n-1]);
      } else {
        corr[i] = .5 * hist[i-1];
      }
    } else {
      corr[i] = .5 * (hist[i-1] + hist[i]);
    }
    if (corr[i] < 1e-3) {
      fprintf(stderr, "corr[%d] = %g is too small, cannot perform correction.\n",
          i, corr[i]);
      goto ERR;
    }
  }

  /* compute the correction */
  for (min = 1e9, i = 0; i <= n; i++) {
    if (distref[i] < 1e-15) {
      fprintf(stderr, "distref[%d] = %g is too small, no correction\n",
          i, distref[i]);
      goto ERR;
    }
    corr[i] = log(corr[i] / distref[i]);
    if (corr[i] < min) min = corr[i];
  }
  /* shift it properly */
  for (maxcorr = 0, i = 0; i <= n; i++) {
    corr[i] -= min;
    if (corr[i] > maxcorr) maxcorr = corr[i];
  }

  /* correct the pmf */
  for (i = 0; i <= n; i++)
    npmf[i] -= corr[i] / bl;
  fprintf(stderr, "successfully applied correction of max = %g/%g = %g\n",
      maxcorr, bl, maxcorr/bl);
  return corr;
ERR:
  return NULL;
}

static int do_fitting(spbonds_t *bs)
{
  const int parts = 3;
  int id, i, j;
  spb_t *spb;
  fit_t farr[PARTMAX], *f;
  double baseline, bl;
  double *arrbuf, *origin, *output, *out0, *corr;

  FILE *fplog;
  static char initfit[PARTMAX+1][LCMAX]={{'\0'}};
  int hasprev = 0;

  if (parts > PARTMAX) {
    fprintf(stderr, "too many %d parts!\n", parts);
    return -1;
  }
  if (bs->cnt >= PARTMAX+1) {
    fprintf(stderr, "too many spbs %d\n", bs->cnt);
    return -1;
  }

  if (readlines(initfit, PARTMAX+1, fnfitlog) >= bs->cnt+1) {
    hasprev = 1;
  } else {
    fprintf(stderr, "insufficient # of lines in %s\n", fnfitlog);
    hasprev = 0;
  }

  for (id = 0; id < bs->cnt; id++) {
    spb = bs->arr+id;

    /* set up the basic array  */
    bl = spb->sbl[0]/spb->hist[0];
    if ((arrbuf = calloc(spb->bins + 1, sizeof arrbuf[0])) == NULL) { exit(1); }
    if ((origin = calloc(spb->bins + 1, sizeof origin[0])) == NULL) { exit(1); }
    if ((output = calloc(spb->bins + 1, sizeof output[0])) == NULL) { exit(1); }
    if ((out0   = calloc(spb->bins + 1, sizeof out0[0]  )) == NULL) { exit(1); }
    for (i = 0; i <= spb->bins; i++) {
      origin[i] = -spb->pmf[i]-log(spb->distref[i])/bl;
    }
    /* replace distref by the actually observed histogram */
    if (add_corr) {
      corr = fit_applycorr(origin, spb->bins, spb->hist, spb->distref,
          bl, spb->flags & SPB_PERIODIC);
    } else corr = NULL;
    /* calculate the baseline and shift the potential */
    baseline = fit_calcbaseline(origin, spb->bins, 1);
    for (i = 0; i <= spb->bins; i++) arrbuf[i] = origin[i];

    for (j = 0; j < parts; j++) {
      f = farr+j;
      f->id = j;
      f->type = id;
      f->baseline = baseline;
      f->arr = arrbuf; /* set up the array address */
      f->dx = spb->binw;
      f->xmin = spb->min;
      f->xmax = spb->max;
      f->period = spb->bins;
      f->xperiod = f->period*f->dx;
      f->bl = bl;
    }

    /* try to load previous fitting data */
    if (hasprev) {
      if (0 != fit3_load(farr, parts, id, initfit[id+1])) {
        fprintf(stderr, "error in using prev. data %d\n[%s]\n",
          id, initfit[id+1]);
        hasprev = 0;
      } else {
        double err = fit3_geterr(farr, parts, origin);
        if (err > 100.0) {
          fprintf(stderr, "error of previous data %g is too large\n", err);
          hasprev = 0;
          exit(1);
        }
      }
    }

    if (hasprev) {
      fprintf(stderr, "using previous data, id = %d, amp = %g, %g, %g; xc = %g, %g, %g\n",
          id, farr[0].a[ID_A], farr[1].a[ID_A], farr[2].a[ID_A],
          farr[0].a[ID_XC], farr[1].a[ID_XC], farr[2].a[ID_XC]);
    } else {
      /* suggest the initial range */
      if (spb->natoms > 0 && strcmp(spb->atoms[0],"C") == 0) {
        fprintf(stderr, "doing the phi angle ... ");
        for (i = 0; i < spb->natoms; i++)
          fprintf(stderr, "%s%s", (i == 0) ? " " : ":", spb->atoms[i]);
        fprintf(stderr, "\n");

        /* set up the range to find a peak */
        fit_peak(farr+0, 1.6, 3.0, 0);  /* the alpha peak first */
        fit_peak(farr+2, -2.9, -1.0, 0); /* the anti-alpha peak next */
        fit_peak(farr+1, farr->xmin, farr->xmax, 1); /* the rest (turn) is beta */

      } else {
        fprintf(stderr, "doing the psi angle ... ");
        for (i = 0; i < spb->natoms; i++)
          fprintf(stderr, "%s%s", (i == 0) ? " " : ":", spb->atoms[i]);
        fprintf(stderr, "\n");

        fit_peak(farr+0, 1.0, -3.0, 0); /* the alpha peak */
        fit_peak(farr+1, -1.0, 1.0, 0); /* the beta peak is prominent */
        fit_peak(farr+2, 1.0, 0.5, 1); /* the rest (turn), we allow it to wrap around */
      }
    }

    /* use the original buffer as the reference,
     * arrbuf to save the remaining */
    fit3_refine(farr, parts, origin, arrbuf); /* do a refinement */

    /* scale the potential according to the relative scaling parameters
     * `output' is the negative logarithm of the new reference distribution */
    fit3_scalepmf(farr, parts, output, out0, arrbuf, spb->popscal, corr, spb->barred);
    /* we only set the reference distribution and mfref,
     * save it the out put binary file,
     * and the rest can be straighten out by md2 (using md2spb module) */
    for (i = 0; i < spb->bins; i++) {
      spb->mfref[i] = -(output[i+1] - output[i])/spb->binw;
    }
    for (i = 0; i <= spb->bins; i++) {
      spb->distref[i] = exp(-output[i]);
      spb->distref0[i] = exp(-out0[i]);
    }

    fit3_save(farr, parts, id, initfit[id+1], LCMAX);
    fit3_print(farr, parts, farr->baseline);
    fprintf(stderr, "\n\n");

    free(arrbuf);
    free(origin);
    free(output);
    free(out0);
    free(corr);
  }
  fprintf(stderr, "updating mf and pmf...\n");
  /* update mf and pmf, now that we have mfref */
  spbs_calcmfpmf(bs, SPB_PMFALL);

  fprintf(stderr, "writing log %s\n", fnfitlog);
  if ((fplog = fopen(fnfitlog, "w")) == NULL) {
    fprintf(stderr, "cannot write fit log to %s.\n", fnfitlog);
  }
  fprintf(fplog, "# %d\n", bs->cnt);
  for (id = 0; id < bs->cnt; id++) {
    fprintf(fplog, "%s\n", initfit[id+1]);
  }
  fclose(fplog);
  return 0;
}

static void pconvarr(double arr[], int n, int h)
{
  int nh = n/2, i;
  double x;

  for (i = 0; i < nh; i++) {
    x = arr[i], arr[i] = arr[i+nh], arr[i+nh] = x;
  }
  if (!h) arr[n] = arr[0];
}

/* conversion for polymer convention */
static void pconvert(spbonds_t *bs)
{
  spb_t *spb;
  int k, n;

  for (k = 0; k < bs->cnt; k++) {
    spb = bs->arr + k;
    n = spb->bins;
    pconvarr(spb->distref0, n, 1);
    pconvarr(spb->distref, n, 1);
    pconvarr(spb->pmf, n, 1);
    pconvarr(spb->mf, n, 0);
    pconvarr(spb->mfref, n, 0);
    pconvarr(spb->hist, n, 0);
    pconvarr(spb->sbf, n, 0);
    pconvarr(spb->sdiv, n, 0);
    pconvarr(spb->sbl, n, 0);
  }
}

int main(int argc, char **argv)
{
  spbonds_t *bs = NULL;

  if (0 != doargs(argc, argv)) return 1;
  if ((bs = spbinit()) == NULL) return 1;
  if (fitting) do_fitting(bs);
  if (polymer) pconvert(bs);
  spbs_write(bs, fntxtout, opver, /* as is */ 1);
  spbs_writebin(bs, fnbinout, opver);
  spbs_close(bs);
  ssdelall();
  return 0;
}


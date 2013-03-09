/* display dihedral information for an .xtc file
 * Copyright (C) 2010 Cheng Zhang */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "typedefs.h"
#include "smalloc.h"
#include "tpxio.h"     /* topology .tpr */
#include "statutil.h"
#include "physics.h"

#define ZCOM_PICK
#define ZCOM_RV3
#include "zcom.h"

const char *fntrj = "trim.xtc";
const char *fntop = "naked.tpr";
const char *fnhist = "dih.hist";
int polymer = 0;  /* polymer convention */
int genhist = 0;  /* generate histogram */

int fptrj = 0; /* gromacs file handle for trj */
rvec *xref = NULL;
int natoms, epbc;

#define BBMAX 4096
typedef struct {
  int ia[5];
  int resid;
  const char *resname;
} bb_t;
int bbcnt = 0;
bb_t bb[BBMAX];
double phir[BBMAX], psir[BBMAX];  /* reference */
double phin[BBMAX], psin[BBMAX];  /* current */
double phid[BBMAX], psid[BBMAX];  /* sum of deviation from reference */
double phi1[BBMAX], psi1[BBMAX];  /* sum */
double phi2[BBMAX], psi2[BBMAX];  /* sum of square */

int devminid = -1;
double devmin = 1e9, devmin1, devmin2;

#define R2D (180.0/M_PI)

#define BINS 360
#define HMIN (-180.)
#define HMAX 180.
#define HDEL (((HMAX - HMIN)/(BINS))+1e-15)
double Hphi[BINS+1], Hpsi[BINS+1];

/* translate dihedral (in degrees) to bin */
static int x2i(double x)
{
  int i = (int)((x-HMIN)/HDEL);
  if (i < 0 || i >= BINS) {
    fprintf(stderr, "x %g out of range!\n", x);
    exit(1);
  }
  return i;
}

/* write histogram */
int hist(const char *fname)
{
  int i;
  FILE *fp;
  double sphi = 0.0, spsi = 0.0;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  for (sphi = spsi = 0., i = 0; i <= BINS; i++) {
    sphi += Hphi[i];
    spsi += Hpsi[i];
  }
  for (i = 0; i <= BINS; i++) {
    double x = HMIN + HDEL * (i+.5);
    fprintf(fp, "%+8.3f %12.10f %12.10f\n",
        x, Hphi[i]/sphi, Hpsi[i]/spsi);
  }
  fclose(fp);
  return 0;
}

/* add backbone index-set from moltype imt */
static void addbackbone(gmx_mtop_t *mtop, int imt, int ff[],
    int resid, const char *resname, int verbose)
{
  int i, ag, im, imb, nmb = mtop->nmolblock;
  gmx_molblock_t *mb;

  if (verbose)
    printf("%3d: adding %d, %d, %d, %d, %d\n",
      bbcnt+1, ff[0], ff[1], ff[2], ff[3], ff[4]);
  /* search molblock */
  for (ag = 0, imb = 0; imb < nmb; imb++) {
    mb = mtop->molblock + imb;
    if (mb->type == imt) {
      for (im = 0; im < mb->nmol; im++) {
        if (bbcnt >= BBMAX) {
          fprintf(stderr, "too many backbones atoms!\n");
          return;
        }
        for (i = 0; i < 5; i++) {
          bb[bbcnt].ia[i] = ag + im*mb->natoms_mol + ff[i];
          bb[bbcnt].resid = resid;
          bb[bbcnt].resname = resname;
        }
        bbcnt++;
      }
    }
    ag += mb->nmol * mb->natoms_mol;
  }
}

/* obtain backbone dihedral indices */
static void dihindex(gmx_mtop_t *mtop)
{
  static const char *bbatms[5] = {"C", "N", "CA", "C", "N"};
  gmx_moltype_t *mt;
  int start, ff[5], found;
  int imt, i, i2, resid, j;
  char *resname;

  /* loop over moltypes to search a interaction that matches `allatoms' */
  for (imt = 0; imt < mtop->nmoltype; imt++) {
    mt = mtop->moltype + imt;
    for (i = 0; i < mt->atoms.nr; ) {
      start = i;
      /* try to find successive dihedral atom set */
      for (j = 0; j < 5; j++) {
        for (found = 0, i2 = start; i2 < mt->atoms.nr; i2++) {
          if (strcmp(bbatms[j], mt->atoms.atomname[i2][0]) == 0) {
            found = 1;
            break;
          }
        }
        if (!found) break;
        ff[j] = i2;
        start = i2 + 1;
      }
      if (j != 5) break; /* exhausted */
      resid = mt->atoms.atom[ff[2]].resnr;
      resname = mt->atoms.resname[resid][0];
      addbackbone(mtop, imt, ff, resid, resname, 0);
      i = ff[0] + 1;
    }
  } /* loop over moltypes */
}

static double convert(double x)
{
  if (polymer)
    return (x > 0) ? (x - 180.) : (x + 180.);
  else return x;
}

/* compute dihedrals in degree, polymer convension if specified */
static int calcdih(rvec r[], double phi[], double psi[])
{
  int i;

  if (r == NULL) return -1;
  for (i = 0; i < bbcnt; i++) {
    phi[i] = convert( R2D * rv3_calcdihv(NULL, r, bb[i].ia, 0) );
    psi[i] = convert( R2D * rv3_calcdihv(NULL, r, bb[i].ia+1, 0) );
  }
  return 0;
}

/* wrap around the difference */
static double devwrap(double x)
{
  double y;
  if (x < 0) y = 360 + x;
  else y = x - 360;
  return (fabs(x) < fabs(y)) ? x : y;
}

/* turn the difference to string */
static char *dev2s(char *s, double x)
{
  int i, n, max = 7;
  char ch;

  if (x < 0.) {
    ch = '-';
    x = -x;
  } else {
    ch = '+';
  }
  n = (int) (x/30.0 + .5);
  if (n > max) n = max;
  for (i = 0; i < n; i++) s[i] = ch;
  for (i = n; i < max; i++) s[i] = ' ';
  s[max] = '\0';
  return s;
}

/* accumulate data */
static int adddih(void)
{
  static int id;
  int i;
  double x, dx, dev1 = 0., dev2 = 0., dev;

  for (i = 0; i < bbcnt; i++) {
    x = phin[i];
    Hphi[x2i(x)] += 1.0;
    phi1[i] += x;
    phi2[i] += x*x;
    phid[i] += dx = devwrap(x - phir[i]);
    dev1 += fabs(dx);

    x = psin[i];
    Hpsi[x2i(x)] += 1.0;
    psi1[i] += x;
    psi2[i] += x*x;
    psid[i] += dx = devwrap(x - psir[i]);
    dev2 += fabs(dx);
  }
  dev1 /= bbcnt;
  dev2 /= bbcnt;
  dev = dev1 + dev2;
  if (dev < devmin) {
    devminid = id;
    devmin = dev;
    devmin1 = dev1;
    devmin2 = dev2;
  }
  id++;
  return 0;
}

/* print dihedral information */
static void print(int nr)
{
  int i;
  double pc, sc, phifluc, psifluc, dphi, dpsi;
  double dphia = 0.0, dpsia = 0.0, dphib = 0.0, dpsib = 0.0;
  char sphi[16], spsi[16];

  printf("# rid  phiref    phiav   phidev  phifluc  "
             "   psiref    psiav   psidev  psifluc  res  phidev  psidev\n");
  for (i = 0; i < bbcnt; i++) {
    pc = phi1[i]/nr;
    sc = psi1[i]/nr;
    dphi = phid[i]/nr;
    dpsi = psid[i]/nr;
    phifluc = sqrt(phi2[i]/nr - pc*pc);
    psifluc = sqrt(psi2[i]/nr - sc*sc);
    printf("%4d %8.3f %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f %8.3f  %-4s %s %s\n",
        bb[i].resid + 1,
        phir[i], pc, dphi, phifluc,
        psir[i], sc, dpsi, psifluc,
        bb[i].resname,
        dev2s(sphi, dphi), dev2s(spsi, dpsi));
    dphia += dphi;
    dpsia += dpsi;
    dphib += fabs(dphi);
    dpsib += fabs(dpsi);
  }
  fprintf(stderr, "dphi = %g (%g), dpsi = %g (%g)\n",
      dphia/bbcnt, dphib/bbcnt,
      dpsia/bbcnt, dpsib/bbcnt);
  fprintf(stderr, "min = %g, %g, %g (frame %d)\n",
      devmin1, devmin2, devmin, devminid);
}

/* read topology and configuration from .tpr file */
static gmx_mtop_t *read_tpx_conf(const char *fn, rvec **x, matrix box,
    int *n)
{
  t_tpxheader hdr;
  int ver, gen, step;
  real t, lambda;
  gmx_mtop_t *mtop = NULL;

  read_tpxheader(fn, &hdr, 1, &ver, &gen);
  snew(*x, hdr.natoms);
  snew(mtop, 1);
  read_tpx(fn, &step, &t, &lambda, NULL, box, n,
      *x, NULL, NULL, mtop);
  return mtop;
}

/* print usage and then quit */
static void usage(const char *prog)
{
  fprintf(stderr, "%s [OPTIONS] file.xtc\n", prog);
  fprintf(stderr,
      " -s: followed by .tpr\n"
      " -f: followed by trajectory .xtc\n"
      " -p: polymer convention\n"
      " -H: generate histogram\n"
      " -d: followed by the name of histogram file\n"
      "\n");
  exit(1);
}

/* set parameters from command line arguments */
static int doargs(int argc, char **argv)
{
  int i, j, ch;
  const char *param;

  for (i = 1; i < argc;  ) {
    if (argv[i][0] != '-') {
      fntrj = argv[i];
    } else if (argv[i][1] != '-') { /* simple options */
      ch = argv[i][1];
      if (strchr("dfs", ch) != NULL) {
        param = argv[i] + 2;
        if (param[0] == '\0') {
          if (i >= argc - 1)
            usage(argv[0]);
          param = argv[++i];
        }
        if (ch == 'f') {
          fntrj = param;
        } else if (ch == 's') {
          fntop = param;
        } else if (ch == 'd') {
          fnhist = param;
        }
      } else { /* switches */
        for (j = 1; (ch = argv[i][j]) != '\0'; j++)
          switch (ch) {
          case 'H': genhist = 1; break;
          case 'p': polymer = 1; break;
          case 'h':
          default: usage(argv[0]); break;
          }
      }
    }
    i++;
  }
  return 0;
}

int main(int argc, char **argv)
{
  gmx_mtop_t *mtop;
  int nfr, natoms0;
  rvec *xnow = NULL;
  matrix box;
  real t = 0.f;

  (void) argtp;

  if (0 != doargs(argc, argv)) return -1;

  /* read topology */
  mtop = read_tpx_conf(fntop, &xref, box, &natoms0);
  /* get indices for backbone dihedrals */
  dihindex(mtop);
  calcdih(xref, phir, psir); /* compute reference */

  for (nfr = 0; ; nfr++) {
    if (nfr == 0) {
      natoms = read_first_x(&fptrj, fntrj, &t, &xnow, box);
    } else {
      if (!read_next_x(fptrj, &t, natoms, xnow, box))
        break;
    }
    calcdih(xnow, phin, psin);
    adddih();
  }
  close_trj(fptrj);
  /* print stat. results */
  print(nfr);
  hist(fnhist);
  fprintf(stderr, "%d atoms, %d frames\n", natoms, nfr);
  sfree(mtop);
  return 0;
}


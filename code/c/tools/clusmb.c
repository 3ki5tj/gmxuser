/* clusterize information for an .xtc file
 * Copyright (C) 2011 Cheng Zhang */
#include "zeutil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_RNG
#define ZCOM_UTIL
#define ZCOM_RV3
#define ZCOM_SS
#define ZCOM_ARGOPT
#include "zcom2.h"
#include "clus.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnout = "clus.txt";
double tequil = 10000;/* handle arguments */

double Tref = 300; /* reference temperature */
double delbet = 1.0; /* we don't skip frames by default */
double rcutoff = 0.07; /* in nm */
int itermax = 0; /* maximal number of steps */
int every = 100;
int ctmd = 0;
int clusmethod = 1;
double weightmin = 0.05;
int beta_cnt = 100;
double beta_min = 0.3, beta_max = 5.0;
int nstmin = 1000;
int doprev = 0;

static int doargs(int argc, char **argv)
{
  int prev1 = 0, prev2 = 0;
  argopt_t *ao = argopt_open(0);

  ao->desc = "compute root-mean-square deviation from the native";
  argopt_add(ao, "-s", NULL,  &fntop,  "topology file");
  argopt_add(ao, "-f", NULL,  &fnxtc,  "trajectory");
  argopt_add(ao, "-z", NULL,  &fnze,   "ZE file");
  argopt_add(ao, "-t", NULL,  &fntr,   "TRACE file");
  argopt_add(ao, "-o", NULL,  &fnout,  "output file");
  argopt_add(ao, "-T", "%lf", &Tref,   "reference temperature");
  argopt_add(ao, "-e", "%d",  &every,  "skip interval");
  argopt_add(ao, "-c", "%lf", &rcutoff, "cutoff distance mu in nm");
  argopt_add(ao, "-I", "%d",  &itermax, "maximal number of iterations");
  argopt_add(ao, "-m", "%d",  &clusmethod, "cluster sampling method, 0: metropolis, 1: heatbath");
  argopt_add(ao, "-0", "%lf", &beta_min, "minimal beta");
  argopt_add(ao, "-1", "%lf", &beta_max, "maximal beta");
  argopt_add(ao, "-#", "%d",  &beta_cnt, "number of beta steps");
  argopt_add(ao, "-P", "%b",  &prev1, "do previous database, from minimal");
  argopt_add(ao, "-R", "%b",  &prev2, "do previous database, from fresh");
  argopt_add(ao, "-q", "%lf", &tequil, "time of equilibration");
  argopt_add(ao, "-M", "%b", &ctmd, "constant temperature simulation");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
  if (prev1) doprev = 1;
  if (prev2) doprev = 2;

  /* check if it is constant temperature md */
  if (!ctmd && !fexists("data1/ZE")) {
    fprintf(stderr, "no ZE, assume constant temperature md\n");
    ctmd = 1;
  }
  return 0;
}

typedef struct {
  int ia[3];
  int resid;
  const char *resname;
} bb_t;
real *vmass; /* mass vector */

/* a single .xtc frame */
typedef struct {
  real (*x)[3];
  double t; /* time, should be able to identify frame uniquely */
  double w;
  const char *fn; /* file name */
} xframe_t;

/* a dynamics array of frames; */
typedef struct {
  xframe_t *fr;
  int nat; /* # of atoms in each frame */
  int nfr; /* # of frames */
  int cap; /* capacity of frames */
} xmovie_t;
#define FRBLOCSIZ 128

/* allocate initial movie structure */
static xmovie_t *xmov_init(int nat)
{
  xmovie_t *xmov;
  xframe_t *fr;
  int i;

  xnew(xmov, 1);
  xmov->nat = nat;
  xmov->nfr = 0;
  xmov->cap = FRBLOCSIZ;
  xnew(xmov->fr, xmov->cap);
  for (i = 0; i < xmov->cap; i++) {
    fr = xmov->fr + i;
    xnew(fr->x, xmov->nat);
    fr->t = fr->w = 0.;
  }
  return xmov;
}

/* possible increase capacity to contain xmov->nfr + 1 frames */
static void xmov_resize(xmovie_t *xmov)
{
  int i0, i;
  xframe_t *fr;
  if (xmov->nfr + 1 <= xmov->cap) return;
  i0 = xmov->cap;
  xmov->cap += FRBLOCSIZ;
  xrenew(xmov->fr, xmov->cap);
  for (i = i0; i < xmov->cap; i++) {
    fr = xmov->fr + i;
    xnew(fr->x, xmov->nat);
    fr->t = fr->w = 0.;
  }
}

static void xmov_free(xmovie_t *xmov)
{
  int i;

  for (i = 0; i < xmov->cap; i++) {
    if (xmov->fr[i].x != NULL)
      free(xmov->fr[i].x);
  }
  free(xmov->fr);
  memset(xmov, 0, sizeof(*xmov));
  free(xmov);
}

/* add backbone index-set from moltype imt */
static bb_t *addbackbone(bb_t *bb, int *rescnt, gmx_mtop_t *mtop, int imt, int ff[],
    int resid, const char *resname, int verbose)
{
  int i, ag, im, imb, nmb = mtop->nmolblock, nres;
  gmx_molblock_t *mb;

  nres = *rescnt;
  if (verbose)
    printf("%3d: adding %d, %d, %d\n", nres+1, ff[0], ff[1], ff[2]);
  /* search molblock for those whose moltype is imt */
  for (ag = 0, imb = 0; imb < nmb; imb++) {
    mb = mtop->molblock + imb;
    if (mb->type == imt) {
      for (im = 0; im < mb->nmol; im++) {
        xrenew(bb, nres+1);
        for (i = 0; i < 3; i++) {
          bb[nres].ia[i] = ag + im*mb->natoms_mol + ff[i];
          bb[nres].resid = resid;
          bb[nres].resname = resname;
        }
        nres++;
      }
    }
    ag += mb->nmol * mb->natoms_mol;
  }
  *rescnt = nres;
  return bb;
}

/* obtain backbone atom indices */
static bb_t *init_index(int *rescnt, gmx_mtop_t *mtop)
{
  static const char *bbatms[3] = {"N", "CA", "C"};
  gmx_moltype_t *mt;
  int start, ff[3], found;
  int imt, i, i2, resid, j;
  char *resname;
  bb_t *bb = NULL;

  xnew(bb, 1);
  *rescnt = 0;
  /* loop over moltypes to search a interaction that matches `allatoms' */
  for (imt = 0; imt < mtop->nmoltype; imt++) {
    mt = mtop->moltype + imt;
    for (i = 0; i < mt->atoms.nr; ) {
      start = i;
      /* try to find successive dihedral atom set */
      for (j = 0; j < 3; j++) {
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
      if (j != 3) break; /* exhausted */
      resid = mt->atoms.atom[ff[1]].resnr;
      resname = mt->atoms.resname[resid][0];
      bb = addbackbone(bb, rescnt, mtop, imt, ff, resid, resname, 0);
      i = ff[0] + 1;
    }
  } /* loop over moltypes */
  printf("%d residues\n", *rescnt);
  /* initialize the mass matrix */
  xnew(vmass, (*rescnt)*3);
  for (i = 0; i < *rescnt; i++) {
    vmass[3*i] = 14.01;
    vmass[3*i+1] = vmass[3*i+2] = 12.01;
  }
  return bb;
}

/* scan trajectory and accumulate data */
static int dotrj(bb_t *bb, int rescnt, trj_t *trj, const char *fxtc,
    int nat, matrix box, xmovie_t *xmov)
{
  int nfr, i, match, natoms;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL, xc;
  double t1, w;
  real t = 0.f, wtot = 0.f;
  static int ifr = 0;
  char *fxtc_copy;

  xnew(fxtc_copy, strlen(fxtc) + 1);
  strcpy(fxtc_copy, fxtc);
  /* loop over xtc frames */
  for (nfr = 0; ; nfr++) {
    if (nfr == 0) {
      natoms = read_first_x(&fpxtc, fxtc, &t, &x, box);
      if (natoms != nat) {
        fprintf(stderr, "natoms mismatch! top %d, xtc %d\n", nat, natoms);
        return -1;
      }
    } else {
      if (!read_next_x(fpxtc, &t, natoms, x, box))
        break;
    }

    if (ctmd) {
      w = 1.;
    } else {
      /* search trj for a matching frame */
      for (match = 0; trj->pos < trj->n; trj->pos++) {
        t1 = trj->fr[trj->pos].t;
        if (fabs(t1 -t) < 0.5) {
          match = 1;
          break;
        } else if (t1 > t + 0.5) {
          fprintf(stderr, "time difference %g (trace), %g (xtc)\n", t1, t);
          return -1;
        }
      }
      if (!match) {
        fprintf(stderr, "cannot find a frame for t = %g\n", t);
        break;
      }
      if (!trj->fr[trj->pos].in) {
        w = 0.0;
      } else {
        w = exp(trj->fr[trj->pos].lnw);
      }
    }

    /* we drop a frame if its weight is too low
     * we also choose weight in a random fashion */
    if (w > weightmin && (++ifr) % every == 0) {
      xframe_t *fr;
      xmov_resize(xmov);
      fr = xmov->fr + xmov->nfr;
      fr->t = t;
      fr->w = w;
      fr->fn = fxtc_copy;
      /* copy coordinates of 3xrescnt backbone atoms */
      for (i = 0; i < rescnt; i++) {
        rv3_copy(fr->x[3*i],   x[bb[i].ia[0]]);
        rv3_copy(fr->x[3*i+1], x[bb[i].ia[1]]);
        rv3_copy(fr->x[3*i+2], x[bb[i].ia[2]]);
      }
      /* remove the center of mass */
      rv3_zero(xc);
      for (wtot = 0.f, i = 0; i < rescnt * 3; i++) {
        rv3_sinc(xc, fr->x[i], vmass[i]);
        wtot += vmass[i];
      }
      rv3_smul(xc, 1.0f/wtot);
      for (i = 0; i < rescnt * 3; i++)
        rv3_dec(fr->x[i], xc);
      xmov->nfr++;
    }
  }
  close_xtc(fpxtc);
  return 0;
}

/* compute the root-mean-square deviation matrix (backbone)
 * mat[i][j] with i < j */
static float **calcrmsmat(xmovie_t *xmov)
{
  float **mat;
  int i, j, k, nat, nfr, npr, np;
  real (*x1)[3], (*x2)[3];

  /* allocate space for the distance matrix
   * only upper-triangle of the matrix */
  mat = dismat_alloc(xmov->nfr);

  /* compute the distance matrix */
  nat = xmov->nat;
  nfr = xmov->nfr;
  npr = nfr*(nfr-1)/2;
  for (np = 0, i = 0; i < nfr-1; i++) {
    x1 = xmov->fr[i].x;
    for (j = i+1; j < nfr; j++) {
      x2 = xmov->fr[j].x;
/*
      // the code only works if the center of mass of the two structures coincide
      for (k = 0; k < nat; k++)
        rv3_copy(x1b[k], x1[k]);
      do_fit(nat, vmass, x2, x1b);
      mat[i][j] = rmsdev(nat, vmass, x2, x1b);
*/
      mat[i][j] = rv3_rmsd(x1, NULL, x2, vmass, nat, 0, NULL, NULL);
      if (++np % 1000 == 0)
        printf("computing rmsd for %d, %d %d/%d = %g%%; \r", i, j, np, npr, 100.0*np/npr);
    }
  }
  printf("finished computing rmsd %d.          \n", nfr);
  return mat;
}

/* write frame information in cluster output */
static void writeclsheader(FILE *fp, const clsys_t *cls, const void *data)
{
  int i;
  const xmovie_t *xmov = (const xmovie_t *)data;
  const char *fn;
  xframe_t *fr;

  die_if (xmov->nfr != cls->np, "# of atoms mismatch, %d vs. %d", xmov->nfr, cls->np);
  for (i = 0; i < xmov->nfr; i++) {
    fr = xmov->fr + i;
    fn = (i == 0 || fr->fn != fr[-1].fn) ?  fr->fn : "";
    fprintf(fp, "# %d %g %.0f %s\n", i, fr->w, fr->t, fn);
  }
}

/* read frame information from file */
static void readclsheader(FILE *fp, clsys_t *cls, void *pdata)
{
  xmovie_t *xmov = (xmovie_t *)pdata;
  char buf[1024], *p, *fnprev = NULL;
  int i, itmp, next;
  xframe_t *fr;

  xmov->nfr = cls->np;
  xmov->cap = xmov->nfr;
  xnew(xmov->fr, xmov->cap);
  for (i = 0; i < cls->np; i++) {
    fr = xmov->fr + i;
    fr->x = NULL;
    die_if (NULL == fgets(buf, sizeof buf, fp),
      "cannot read %d frame\n", i);
    strip(buf);
    die_if (3 != sscanf(buf, "#%d%lf%lf%n", &itmp, &fr->w, &fr->t, &next),
        "error read %d frame [%s]\n", i, buf);
    die_if (i != itmp, "frame mismatch %d vs. %d", i, itmp);
    p = buf+next;
    strip(p);
    if (*p != '\0') fnprev = ssdup(p);
    fr->fn = fnprev;
  }
}

/* optimize cluster */
static void opticlus(clsys_t *cls, xmovie_t *xmov)
{
  int verbose = 2;
  cls_zalgo(cls, itermax, clusmethod, beta_min, beta_max, beta_cnt, nstmin, verbose);
  printf("after optimization %d clusters, energy = %g\n", cls->nc, cls->ene);
  /* also find centroid */
  cls_write(cls, fnout, &writeclsheader, xmov, 0);
  cls_free(cls, 1);
}

/* do clustering from saved frames */
static int doclus(xmovie_t *xmov)
{
  float **rmsmat;
  double *warr;
  clsys_t *cls;
  int i;

  /* constructing RMSD */
  printf("RMSD matrix size = %gM (n = %d)\n",
      1e-6*sizeof(rmsmat[0][0])*xmov->nfr*(xmov->nfr+1)/2, xmov->nfr);
  rmsmat = calcrmsmat(xmov);
  /* make a weight array */
  xnew(warr, xmov->nfr);
  for (i = 0; i < xmov->nfr; i++) warr[i] = xmov->fr[i].w;

  /* do the clustering */
  cls = cls_init(rmsmat, warr, xmov->nfr, rcutoff);
  opticlus(cls, xmov);

  return 0;
}

static int domoltrj(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat, rescnt;
  char fxtc[FILENAME_MAX], *ftop;
  matrix box;
  rvec *xref = NULL;
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  bb_t *bb;
  double dt, bref = 1.0/(BOLTZ*Tref), wtot;
  xmovie_t *xmov;

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  bb = init_index(&rescnt, mtop); /* get indices for backbone dihedrals */

  if (!ctmd) { /* initialize weights for tempering by loading TRACE */
    if ((ze = ze_load(fnze, 10)) == NULL) {
      fprintf(stderr, "cannot load %s\n", fnze);
      return -1;
    }
    /* we first load TRACE, compute weight of each frame */
    if ((trj = trj_loadseq(fntr, &ndata, dt, tequil,
            ze, bref, delbet)) == NULL) {
      fprintf(stderr, "failed to load trace\n");
      return -1;
    }
  } else {
    ndata = 10000000;
  }

  /* we now successively load xtc, try to match the TRACE frames */
  if (!ctmd) trj->pos = 0;
  xmov = xmov_init(rescnt*3);
  for (i = 1; i <= ndata; i++) {
    sprintf(fxtc, "data%d/%s", i, fnxtc);
    if (!fexists(fxtc)) { /* test if the file exists */
      die_if(!ctmd, "cannot read %s\n", fnxtc);
      break;
    }
    if (0 != dotrj(bb, rescnt, trj, fxtc, nat, box, xmov)) {
      fprintf(stderr, "error doing %s\n", fnxtc);
      return -1;
    }
  }
  for (wtot = 0., i = 0; i < xmov->nfr; i++) wtot += xmov->fr[i].w;
  printf("%d frames (wtot = %g) for clustering...\n", xmov->nfr, wtot);

  /* compute rmsd and do clustering */
  doclus(xmov);

  xmov_free(xmov);
  sfree(mtop);
  if (!ctmd) {
    if (trj) trj_free(trj);
    if (ze) ze_free(ze);
  }
  return 0;
}

/* load a previous database, then do optimization */
static void optiprev(void)
{
  clsys_t *cls;
  xmovie_t *xmov;

  xnew(xmov, 1);
  cls = cls_read(fnout, readclsheader, xmov);
  die_if (cls == NULL, "cannot load frames\n");
  printf("loaded prev. cluster from %s, %d clusters, ene = %g, mu = %g, np %d, wtot %g\n",
      fnout, cls->nc, cls->ene, cls->mu0, cls->np, cls->wtot);
  if (rcutoff > 0. && fabs(cls->mu0-rcutoff) > 1e-6) { /* reset mu */
    printf("warning: changing clustering mu from %g to %g\n", cls->mu0, rcutoff);
    cls_changemu(cls, rcutoff);
  }
  if (doprev == 2) {
    cls_evensplit(cls, cls->nc);
  }
  opticlus(cls, xmov);
  xmov_free(xmov);
}

int main(int argc, char **argv)
{
  (void) argtp;
  doargs(argc, argv);
  if (doprev) {
    optiprev();
  } else {
    domoltrj();
  }
  return 0;
}


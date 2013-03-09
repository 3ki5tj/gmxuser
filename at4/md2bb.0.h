/* a pair of phi-psi backbone dihedrals */
#ifndef BB_C__
#define BB_C__

#include <stdio.h>
#include <string.h>
#include <math.h>

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_CFG
#define ZCOM_RV3 /* define real */
#define ZCOM_ENDN /* binary i/o */
#define ZCOM_SS
#define ZCOM_DISTR
#include "zcom2.h"

#define BB_ATOMS 5
/* five atoms: i, j, k, l, h */

#define BB_MF        0x0010
#define BB_PMF       0x0040
#define BB_PMFALL    (BB_MF|BB_PMF)
#define BB_VERIFY    0x0080

/* buffer information before the force is available */
typedef struct {
  double phi, psi, mfphi, mfpsi, g2phi, g2psi;
  int idx[BB_ATOMS];
  real gi[3], gl[3]; /* gradient for phi */
  real gj[3], gh[3]; /* gradient for psi */
} bbbuf_t; /* $skip; */

#define BB_PMFSCALS 9

#if defined(GMX_MPI)
  #ifndef BB_MPI
    #define BB_MPI 1
  #endif
#endif

/* if GMX_THREADS is defined, tmpi.h should be used, GROMACS 4.1 */
#if ( defined(BB_MPI) && !defined(GMX_THREADS) )
#include "mpi.h"
#endif

/* let object generator know */
#define USE_MPI BB_MPI

typedef struct {
  /* $kprefix := bb%d_; $kargs := @id;
   * $# cfg key looks like bb0_xxx, bb1_xxx */
  bbs_t *bbs; /* $usr: parent; $# pointer to parent, passed during init */
  int id; /* $usr:cfg; index of this bb */
  unsigned flags; /* flags for global settings $io:; */

  double distmin;    /* minimal relative value for distribution
                        $def: @~; */
  double distmax;    /* maximal relative value for distribution
                        $def: @~; */
  double mfmax;      /* maximal allowable mean force $def: @~;
                        $io:bt; $binprev: @bins; $binprereq: ver > 0; */
  int    pmfside;    /* 1: pmf >= 0;  0: pmf[0] = 0.0;  -1: pmf <= 0; $def: @~;
                        $io:bt; $binprev: @mfmax; $binprereq: ver > 0; */
  double popscal[BB_PMFSCALS];
                     /* population scaling factor for different components,
                        alpha, beta, turn, rest
                        $def: @~[i]; $io:c;
                        $complete; $valid: @@[i] >= 0.0; */
  double bl;  /* beta*lambda, only useful in prep run,
                 which should be constant; $def: @~; */
  /* $io := $0; $# reset */

  /* $io:=$0;  $# reset */
  unsigned rflags;  /* $usr:rb; */
  int    bins;      /* how many bins for phi or psi
                       $def: @~; $io:cbt; $binprev: @id;
                       $rbverify: @rflags & BB_VERIFY; */
  int    bins2;     /* bins^2; $def: @bins * @bins; $io:none; */
  double binw;      /* bin width $def:2*M_PI/@bins; $io:none; */

  /* $prereq := $ismaster; $cnt := @bins2; $clr := 1; */
  double *hist; /* histogram $binprev: @pmfside; */
  double *sbfphi; /* sum of beta * phi force;  $binprev: @hist; */
  double *sbfpsi; /* sum of div.v;  $binprev: @sbfphi; */
  double *sbfphi2;  /* phi^2; $binprev: @sbfpsi; */
  double *sbfpsi2;  /* psi^2; $binprev: @sbfphi2; */
  double *sbfphipsi; /* phi*psi; $binprev: @sbfpsi2; */

#ifdef BB_MPI
  /* local data accumulator, for MPI communication
   * $io := ; $cnt := @bins2; $clr := 1;
   * $mpi := alloc; $sync := 1; $redtmp := @dbltmp; */
  double *hist_l;       /* $reduce: @hist; */
  double *sbfphi_l;     /* $reduce: @sbfphi; */
  double *sbfpsi_l;     /* $reduce: @sbfpsi; */
  double *sbfphi2_l;    /* $reduce: @sbfphi2; */
  double *sbfpsi2_l;    /* $reduce: @sbfpsi2; */
  double *sbfphipsi_l;  /* $reduce: @sbfphipsi; */
  /* $cnt := $0; $io := $0; $clr := $0;
   * $mpi := $0; $sync := $0; $redtmp := $0; $# reset  */
#endif

  double *distref;  /* target distribution,
                       $io:b;
                       $cnt: @bins2; $def: 0.0;
                       $binprev: @sbfphipsi;
                       $rbvalid: @-calcmfref(@) == 0; $# compute mean-force component from distref;
                       $mpi; */
  double *distref0; /* the actual target distribution without correction,
                       $io:b;
                       $cnt: @bins2; $def: 0.0;
                       $binprev: @distref; */
  double *mfphiref; /* target mean force of phi: negative of the logarithm of
                       the target distribution,
                       $cnt: @bins2; $def: 0.0; $prev: @distref;  $io:; $mpi; */
  double *mfpsiref; /* target mean force of psi;
                       $cnt: @bins2; $def: 0.0; $prev: @mfphiref; $io:; $mpi; */

  /* $clr := 1; */
  double *dbltmp;  /* for temporary use, i.e., synchronizing, io
                      $cnt: @bins2; $io:; $prereq = $ismaster; */
  double *mfphi; /* phi force $cnt: @bins2; $io: bt;
                    $binprev: @mfpsiref; $mpi; $bcast; */
  double *mfpsi; /* psi force $cnt: @bins2; $io: bt;
                    $binprev: @mfphi; $mpi; $bcast; */
  double *pmf;   /* p.m.f. $io:none; $cnt = @bins2; $mpi; */

  /* $clr := $0; */

  /* buffer list (temporary atom indices, force, etc)
   * $io := ; $# turn off io */
  int bufcnt; /* number of items in buf */
  int bufcap; /* allocated memory */
  bbbuf_t *buf; /* $cnt=0; $# don't allocation space at the beginning */

  /* GROMACS extensions, they are here to avoid another struct */
  int type;     /* the interaction *type* in gromacs interaction list */
  int exclcnt;  /* # of items */
  int exclcap;  /* capacity of excl[] */
  int *excl;    /* interactions (atom indices) that are considered as bad
                   $cnt=0; $mpicnt = @exclcap * BB_ATOMS; $mpi; */
} bb_t;

typedef struct {
  /* $kprefix := bb_; */
  const char  *iopfx;    /* $usr:cfgref; */
  char  *bin_file;  /* name of the binary file for IO (master)
                       $def = $> @@ = ssdup("")\;\nif (@iopfx != NULL) sscpy(@@, @iopfx)\;\nsscat(@@, "bb.bin");
                       $prereq: $ismaster; $cfgprereq:1; */
  char  *txt_file;  /* name of the text file for IO (master)
                       $def = $> @@ = ssdup("")\;\nif (@iopfx != NULL) sscpy(@@, @iopfx)\;\nsscat(@@, "bb.txt");
                       $prereq: $ismaster; $cfgprereq:1; */

  /* default global settings */
  double distmin;     /* minimal relative value for bb distribution
                         $def=1e-4; */
  double distmax;     /* maximal relative value for bb distribution
                         $def=1.0; $valid: @distmin < @distmax; */
  double mfmax;       /* maximal magnitude */
  double popscal[BB_PMFSCALS]; /* $def = 1.0; $io:c; */
  int pmfside; /* method of determine the baseline of bb pmf
                  $def = -1;  */
  int bins;   /* number of bins $def: 60; */
  double bl;  /* beta * lambda $def: 1.0; */

  int    cnt; /* the number of special bonds (input)
                 $io:cbt; $rbverify; */
  bb_t *arr;  /* array of individual bbs
                 $obj; $cnt: @cnt;
                 $cfgargs: @, i;
                 $rbargs: @rflags;
                 $clr; $mpi; $reduce; $bcast; */
  unsigned rflags; /* $usr:rb; */
} bbs_t; /* $ptrname: bbs; $fprefix: bbs_; */

/* return index from (phi, psi) */
static int bb_idx(bb_t *bb, double phi, int *i, double psi, int *j)
{
  *i = (int)((phi + M_PI) / bb->binw);
  if (*i < 0) *i = 0; else if (*i >= bb->bins) *i = bb->bins - 1;
  *j = (int)((psi + M_PI) / bb->binw);
  if (*j < 0) *j = 0; else if (*j >= bb->bins) *j = bb->bins - 1;
  return (*i) * bb->bins + (*j);
}

/* add an entry, fphi: -dH/dphi, fpsi: -dH/dpsi
 * if the Hamiltonain H = H0 + H1, the function should be
 * called *after* the contribution from H1 is included */
static int bb_add(bb_t *bb, double phi, double fphi,
    double psi, double fpsi, double beta)
{
  int i, j, id;
  double *hist, *sbfphi, *sbfpsi, *sbfphi2, *sbfpsi2, *sbfphipsi;

  id = bb_idx(bb, phi, &i, psi, &j);
#ifdef BB_MPI
  /* in MPI, we add to the local buffer, sync later */
  hist   = bb->hist_l;
  sbfphi = bb->sbfphi_l;
  sbfpsi = bb->sbfpsi_l;
  sbfphi2 = bb->sbfphi2_l;
  sbfpsi2 = bb->sbfpsi2_l;
  sbfphipsi = bb->sbfphipsi_l;
#else
  hist   = bb->hist;
  sbfphi = bb->sbfphi;
  sbfpsi = bb->sbfpsi;
  sbfphi2 = bb->sbfphi2;
  sbfpsi2 = bb->sbfpsi2;
  sbfphipsi = bb->sbfphipsi;
#endif
  hist[id]      += 1.0;
  sbfphi[id]    += (fphi *= beta);
  sbfpsi[id]    += (fpsi *= beta);
  sbfphi2[id]   += fphi * fphi;
  sbfpsi2[id]   += fpsi * fpsi;
  sbfphipsi[id] += fphi * fpsi;
  return 0;
}

/* compute component of mean force from reference distribution */
static int bb_calcmfref(bb_t *bb)
{
  int i, j, ip, jp, id, n = bb->bins;
  double v, vi, vj, vij, fp, fm, dx;

  die_if(bb->bins % 2 != 0, "bins %d should be even\n", bb->bins);

  dx = bb->binw;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      id = i*n + j;
      ip = (i+1) % n;  /* periodic boundary */
      jp = (j+1) % n;
      v   = -log(bb->distref[i*n + j]);
      vj  = -log(bb->distref[i*n + jp]);
      vi  = -log(bb->distref[ip*n + j]);
      vij = -log(bb->distref[ip*n + jp]);
      fp = v - vij;
      fm = vi - vj;
      bb->mfphiref[id] = .5*(fp - fm)/dx;
      bb->mfpsiref[id] = .5*(fp + fm)/dx;
    }
  return 0;
}

/* compute mean force of bin i */
static int bb_calcmf(bb_t *bb)
{
  double num, den, mf;
  int i;

  for (i = 0; i < bb->bins2; i++) {
    den = bb->bl*bb->hist[i];

    /* first offset the intrinsic pmf component
     * then add a bias from the reference distribution */
    num = bb->mfphiref[i]*bb->hist[i] - bb->sbfphi[i];
    mf = (fabs(den) > 0.001) ? (num/den) : 0.0;
    if (mf > bb->mfmax) mf = bb->mfmax;
    else if (mf < -bb->mfmax)  mf = -bb->mfmax;
    bb->mfphi[i] = mf;

    num = bb->mfpsiref[i]*bb->hist[i] - bb->sbfpsi[i];
    mf = (fabs(den) > 0.001) ? (num/den) : 0.0;
    if (mf > bb->mfmax) mf = bb->mfmax;
    else if (mf < -bb->mfmax)  mf = -bb->mfmax;
    bb->mfpsi[i] = mf;
  }
  return 0;
}

/* compute potential of mean force (approximate) */
static int bb_calcpmf(bb_t *bb)
{
  int i, j, im, n = bb->bins;
  double dx = bb->binw;

  bb->pmf[0] = 0.0;
  for (i = 0; i < n; i++) {
    if (i < n - 1) {
      bb->pmf[(i+1)*n] = bb->pmf[i*n]
        - .5*(bb->mfphi[i*n] + bb->mfphi[i*n + n - 1])*dx;
    }
    im = (i-1+n) % n;
    for (j = 0; j < n - 1; j++) {
      bb->pmf[i*n + j + 1] = bb->pmf[i*n + j]
        - .5*(bb->mfpsi[i*n+j] + bb->mfpsi[im*n+j])*dx;
    }
  }
  return 0;
}

static double bb_getpmf(bb_t *bb, double phi, double *mfphi,
    double psi, double *mfpsi)
{
  int i = 0, j = 0, id;
  double dx, dy;

  id = bb_idx(bb, phi, &i, psi, &j);
  dx = phi - (-M_PI + i*bb->binw);
  dy = psi - (-M_PI + j*bb->binw);
  *mfphi = bb->mfphi[id];
  *mfpsi = bb->mfpsi[id];
  return bb->pmf[id] - dx*(*mfphi) - dy*(*mfpsi);
}

/* calculate the mean force and its potential of all bins
 * set flags to a combination of BB_PMF, BB_MF */
static void bbs_calcmfpmf(bbs_t *bbs, unsigned flags)
{
  int id;
  bb_t *bb;

  if (bbs->cnt <= 0) return;

  for (id = 0; id < bbs->cnt; id++) {
    bb = bbs->arr+id;
    if (flags & BB_MF) bb_calcmf(bb); /* mean force */
    if (flags & BB_PMF) bb_calcpmf(bb); /* pmf */
  }
}

/* GROMACS extension, add an atom-set to the exclusion list */
int bb_excl(bb_t *bb, int *ids, int offset)
{
  const int bufsiz = 16;
  int j;

  if (bb == NULL) return -1;
  if (bb->exclcnt >= bb->exclcap) {
    bb->exclcap += bufsiz;
    bb->excl = realloc(bb->excl, bb->exclcap * sizeof(*ids) * BB_ATOMS);
    die_if (bb->excl == NULL, "no memory for excl list %d, size: %d\n",
        bb->id, bb->exclcap);
  }
  for (j = 0; j < BB_ATOMS; j++)
    bb->excl[bb->exclcnt * BB_ATOMS + j] = ids[j] + offset;
  bb->exclcnt++;
  return 0;
}

/* add atom indices, grad, ... to buffer
 * cannot use bb_add now, because force is not yet available */
int bb_buf(bb_t *bb, const int ia[],
    double phi, double mfphi, double g2phi, real gi[3], real gl[3],
    double psi, double mfpsi, double g2psi, real gj[3], real gh[3])
{
  const int blocksiz = 16;
  int j;
  bbbuf_t *cc;

  if (bb == NULL) return -1;
  if (bb->bufcnt >= bb->bufcap) {
    bb->bufcap += blocksiz;
    /* expand the current buffer, if it is out of capacity */
    bb->buf = realloc(bb->buf, bb->bufcap * sizeof(bb->buf[0]));
    die_if (bb->buf == NULL,
      "no memory for buf, %d, size=%d\n", bb->id, bb->bufcap);
  }
  cc = bb->buf + bb->bufcnt;
  cc->phi = phi;
  cc->psi = psi;
  cc->mfphi = mfphi;
  cc->mfpsi = mfpsi;
  cc->g2phi = g2phi;
  cc->g2psi = g2psi;
  for (j = 0; j < BB_ATOMS; j++) {
    cc->idx[j] = ia[j];
  }
  rv3_copy(cc->gi, gi);
  rv3_copy(cc->gl, gl);
  rv3_copy(cc->gj, gj);
  rv3_copy(cc->gh, gh);
  bb->bufcnt++;
  return 0;
}

/* add bufferred items, with force to local data */
int bbs_debuf(bbs_t *bbs, real f[][3], double beta, double lambda)
{
  int id, ic;
  bb_t *bb;
  bbbuf_t *cc;
  double g2, fvphi, fvpsi;

  if (bbs->cnt <= 0) return 0;

  for (id = 0; id < bbs->cnt; id++) {
    bb = bbs->arr + id;
    if (bb->bufcnt <= 0) continue;
    for (ic = 0; ic < bb->bufcnt; ic++) {
      /* compute div + beta * grad(H).grad(x) and beta * lambda */
      cc = bb->buf + ic;

      /* phi use the first atom, psi use the last item
       * to completely eliminate correlation */

      g2 = rv3_dot(cc->gi, cc->gi);
      if (g2 <= 1e-20) continue; /* skip ? */
      fvphi = rv3_dot(f[ cc->idx[0] ], cc->gi)/g2;
      fvphi -= cc->mfphi * lambda;

      g2 = rv3_dot(cc->gh, cc->gh);
      if (g2 <= 1e-20) continue;
      fvpsi = rv3_dot(f[ cc->idx[4] ], cc->gh)/g2;
      fvpsi -= cc->mfpsi * lambda;
      bb_add(bb, cc->phi, fvphi, cc->psi, fvpsi, beta);
    }
    bb->bufcnt = 0; /* clear buffer */
  }
  return 0;
}

/* export to compatiable format */
int bb_export(bb_t *bb, const char *fn)
{
  distr2d_t *d = distr2d_open(-M_PI, M_PI, bb->binw, -M_PI, M_PI, bb->binw);
  int i, j, id, n = bb->bins, m = bb->bins;

  die_if (n != bb->bins || m != bb->bins, "n %d, m %d, bins %d\n", n, m, bb->bins);
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      distr2dsum_t *ds = d->arr + i * (m + 1) + j;
      id = i * n + j;
      ds->s = bb->hist[id];
      ds->sf = bb->sbfphi[id];
      ds->sg = bb->sbfpsi[id];
      ds->sf2 = bb->sbfphi2[id];
      ds->sg2 = bb->sbfpsi2[id];
      ds->sfg = bb->sbfphipsi[id];
    }
  distr2d_save(d, fn);
  distr2d_close(d);
  return 0;
}

int bbs_write(bbs_t *bbs, const char *fn, int ver)
{
  int i, j, id, ib, n;
  double phi, psi;
  bb_t *bb;
  FILE *fp;

  /* compatible format */
  bb_export(bbs->arr, "bbcompat.txt");

  if (fn == NULL) fn = bbs->txt_file;
  xfopen(fp, fn, "w", return -1);

  fprintf(fp, "# %d %d 0 | ", bbs->cnt, ver);
  for (ib = 0; ib < bbs->cnt; ib++) {
    bb = bbs->arr + ib;
    fprintf(fp, "%d %20.10e %d | ",
        bb->bins, bb->mfmax, bb->pmfside);
  }
  fprintf(fp, "\n");

  for (ib = 0; ib < bbs->cnt; ib++) {
    bb = bbs->arr + ib;
    n = bb->bins;

    for (i = 0; i < n; i++) {
      phi = -M_PI + i * bb->binw;
      for (j = 0; j < n; j++) {
        double x, fphi = 0, fpsi = 0, fphi2 = 0, fpsi2 = 0, fphipsi = 0;

        psi = -M_PI + j * bb->binw;
        id = i*n + j;
        if ((x = bb->hist[id]) > .0) {
          fphi = bb->sbfphi[id]/x;
          fpsi = bb->sbfpsi[id]/x;
          fphi2 = bb->sbfphi2[id]/x - fphi*fphi;
          fpsi2 = bb->sbfpsi2[id]/x - fpsi*fpsi;
          fphipsi = bb->sbfphipsi[id]/x - fphi*psi;
        }
        fprintf(fp, "%9.6f %9.6f ", phi, psi);
        fprintf(fp, "%24.17e %24.17e %24.17e %24.17e %24.17e %24.17e "
            "%10.6f %24.17e %24.17e %10.6f %10.6f\n",
            x, fphi, fpsi, fphi2, fpsi2, fphipsi,
            bb->pmf[id], bb->mfphi[id], bb->mfpsi[id],
            bb->distref[id], bb->distref0[id]);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* collects dihedral data, master calculates the mf,
 * bcast mf, every node computes the pmf, mfcom */
int bbs_syncdata(bbs_t *bbs)
{
#ifdef BB_MPI
  int id;
#endif

  if (bbs->cnt <= 0) return 0;

#ifndef BB_MPI /* single processor case */
  bbs_calcmfpmf(bbs, BB_PMFALL);
#else /* MPI version */

  if (0 != bbs_reduce(bbs)) { /* collect local data */
    fprintf(stderr, "error during bbs_reduce: %d\n", bbs->mpi_rank);
    return -1;
  }

  if (bbs->mpi_rank == 0) { /* master computes the mean force */
    for (id = 0; id < bbs->cnt; id++)
      bb_calcmf(bbs->arr + id);
  }

  if (0 != bbs_bcast(bbs)) {  /* bcast sync'ed result */
    fprintf(stderr, "error during bbs_bcast: %d\n", bbs->mpi_rank);
    return -1;
  }

  for (id = 0; id < bbs->cnt; id++) /* compute mf from pmf */
    bb_calcpmf(bbs->arr + id);
#endif /* BB_MPI */
  return 0;
}

/* bb_cleardata: clear bb_t data */
static void bb_cleardata(bb_t *bb)
{
  int i;

  if (bb->mpi_rank == 0) {
    for (i = 0; i < bb->bins2; i++) {
      bb->hist[i] = 0.0;
      bb->sbfphi[i] = 0.0;
      bb->sbfpsi[i] = 0.0;
      bb->sbfphi2[i] = 0.0;
      bb->sbfpsi2[i] = 0.0;
      bb->sbfphipsi[i] = 0.0;
    }
  }
#ifdef BB_MPI
  for (i = 0; i < bb->bins2; i++) {
    bb->hist_l[i] = 0.0;
    bb->sbfphi_l[i] = 0.0;
    bb->sbfpsi_l[i] = 0.0;
    bb->sbfphi2_l[i] = 0.0;
    bb->sbfpsi2_l[i] = 0.0;
    bb->sbfphipsi_l[i] = 0.0;
  }
#endif
}

static void bbs_cleardata(bbs_t *bbs)
{
  int i;
  for (i = 0; i < bbs->cnt; i++)
    bb_cleardata(bbs->arr + i);
}

#endif /* BB_C__ */




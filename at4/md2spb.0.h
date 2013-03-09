/*  handle special bonds */
#ifndef MD_SPB_C__
#define MD_SPB_C__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_CFG
#define ZCOM_RV3 /* define real */
#define ZCOM_ENDN /* binary i/o */
#define ZCOM_SS
#include "zcom2.h"

#ifndef SPB_XFUNCS
#define SPBSTRCLS static
#else
#define SPBSTRCLS
#endif

#ifdef MAXATOMLIST /* GROMACS include/types/idef.h */
#define SPB_MAXNATOMS MAXATOMLIST
#else
#define SPB_MAXNATOMS 6  /* maximal # of atoms involved in an interaction */
#endif

#define SPB_VERBOSE   0x0001
#define SPB_MF        0x0010
#define SPB_MFCOM     0x0020
#define SPB_PMF       0x0040
#define SPB_PMFALL    (SPB_MF|SPB_MFCOM|SPB_PMF)
#define SPB_VERIFY    0x0080
#define SPB_VERREF    0x0100

#define SPB_CSORDER   6

#if defined(GMX_MPI)
  #ifndef SPB_MPI
    #define SPB_MPI 1
  #endif
#endif

/* if GMX_THREADS is defined, tmpi.h should be used, GROMACS 4.1 */
#if ( defined(SPB_MPI) && !defined(GMX_THREADS) )
#include "mpi.h"
#endif

/* let object generator know */
#define USE_MPI SPB_MPI

typedef struct { /* needs to have real defined */
  double x, div, grad2, mfph;
  int    idx[SPB_MAXNATOMS];
  real   grad[SPB_MAXNATOMS][3];
} spbbuf_t; /* $skip; */

#define SPB_PMFSCAL_CNT 4
typedef struct {
  /* $kprefix := spb%d_; $kargs := @id;
   * $# cfg key will look like spb0_xxx, sbp1_xxx, ... */
  spbonds_t *bs;  /* $usr: parent;
                     $# pointer to parent, which is to be passed during init. */
  int   id;       /* $usr:cfg; index of this spb: 0, 1, ...  $io:; */
  int   natoms;   /* # of atoms involved (calculate from input)
                     $io:bt; $binprev: @pmfside; $binprereq: ver > 0; */
  char* atoms[SPB_MAXNATOMS]; /* string of atoms involved (converted from input)
                                 master only, parsed from atmbuf $io:; */
  int   atmsiz; /* size of atmbuf, computed from atmbuf
                   used for binary IO $io:; */
  char* atmbuf; /* the original string of atoms, managed by ss
                   $key: @<atoms; $cfgprereq:1;
                   $valid: @-parse_atoms_txt(@, 1) == 0;
                   $cnt: @atmsiz; $io: cbt;
                   $prereq: $ismaster;
                   $rbcnt: @atmsizb;
                   $rbvalid: @-parse_atoms_bin(@, 1) == 0;
                   $binprev: @atmsizb;
                   $binprereq := ver > 0 && atmsizb > 0;
                   $# initialize natoms, atoms, atmsiz */
  /* $# prepare buffer for writing binary */
  int   atmsizb;    /* multiple of 16 that is not less than atmsiz
                       $usr:bintmp; $io:b; $binprev: @natoms;
                       $def: ((@atmsiz + 15)/16)*16; $rbdef: 0;
                       $binprereq: ver > 0;
                       $rbvalid: (@atmsiz = @@) >= 0 && ((@atmbuf = ssnew(@@)) != NULL); */
  char  atmbufr[16];/* '\0' $def='\0';
                       $usr:wbtmp; $bincnt: @atmsizb - @atmsiz;
                       $io:b; $binprev: @atmbuf; */
  /* $binprereq := $0; */

  unsigned flags;    /* flags for global settings $io:; */
#define SPB_PERIODIC 0x0004 /* $key: seamless; $def: 1; spb is periodic, like an angle or a dihedral */
  /* target distribution tuning, used to override global values
   * $io := c; */
  double barred;     /* factor to reduce barrier, used in md2conv;
                        $def: @~; $valid: @@ >= 1.0; */
  double distmin;    /* minimal relative value for distribution
                        $def: @~; */
  double distmax;    /* maximal relative value for distribution
                        $def: @~; */
  double distpwrscl; /* scaling factor to control the distribution
                        $def: @~; */
  double mfmax;      /* maximal allowable mean force $def: @~;
                        $io:bt; $binprev: @max; $binprereq: ver > 0; */
  int    pmfside;    /* 1: pmf >= 0;  0: pmf[0] = 0.0;  -1: pmf <= 0; $def: @~;
                        $io:bt; $binprev: @mfmax; $binprereq: ver > 0; */
  double popscal[SPB_PMFSCAL_CNT];
                     /* population scaling factor for different components,
                        alpha, beta, turn, rest
                        $def: @~[i]; $io:c;
                        $complete; $valid: @@[i] >= 0.0; */

  int     formula; /* use cosine-sine series; $def: @~; */
  double paramc[SPB_CSORDER]; /* parameters for cosine series $def = 0.0; $io:c; $mpi; */
  double params[SPB_CSORDER]; /* parameters for sine series $def = 0.0; $io:c; $mpi; */
  /* $io:=$0;  $# reset */

  unsigned rflags;  /* $usr:rb; $# parameter of readbin() */
  int    bins;      /* how many bins (input)
                       $io:bt; $binprev: @id;
                       $rbverify: @rflags & SPB_VERIFY; */
  double binw;      /* bin width = 2*pi/bins $io:none; */
  double min;       /* lower boundary of the angle, $def: -M_PI;
                       $io:bt; $binprev: @bins;
                       $rbverify: @rflags & SPB_VERIFY; */
  double max;       /* higher boundary of the angle, $def: M_PI;
                       $io:bt; $binprev: @min;
                       $rbverify: @rflags & SPB_VERIFY; */
  int binrefcnt;    /* $usr:bintmp; $def: @bins + (ver>0); $binprev: @mf; $io:; */
  double *distref;  /* target distribution,
                       $io:b; $rbverify: @rflags & SPB_VERREF;
                       $cfgcnt: 0; $# avoid memory allocate during cfgopen, don't have @bins before parsing sdistref;
                       $cnt: @bins + 1;
                       $binprev: @binrefcnt; $bincnt: @binrefcnt;
                       $rbvalid: @-calcmfref(@) == 0; $# compute mean-force component from distref;
                       $mpi; */
  double *distref0; /* the actual target distribution without correction,
                       $io:b;  $cfgcnt: 0; $# avoid allocation;
                       $cnt: @bins + 1;
                       $binprereq: ver >= 2; $binprev: @distref; */
  /* $rbcall: if (ver < 2) {\nint k\;\n\nfor (k = 0\; k <= @bins\; k++)\n\t@distref0[k] = @distref[k]\;\n}
   * $binprev: @distref0; */
  char   *sdistref; /* $key: @<distref; $io:c; $must;
                       $prereq: $ismaster; $cfgprereq:1;
                       $valid: 0 == @-initdistref(@, @@);
                       $# initialize bins, binw, min, max, distref; */
  double *mfref;   /* target mean force: negative of the logarithm of
                      the target distribution,
                      $cnt:@bins;
                      $def: 0.0;
                      $cfgvalid: 0 == @-calcmfref(@); $cfgprev: @distref;
                      $io:; $mpi; */

  /* $prereq := $ismaster; $cnt := @bins; $clr := 1; */
  double *hist; /* histogram
                   $binprev: @atmbufr; */
  double *sbf;  /* sum of beta* (projection of force to phi)
                   $binprev: @hist; */
  double *sdiv; /* sum of div.v
                   $binprev: @sbf; */
  double *sbl;  /* denominator, sum of Beta*Lambda
                   $binprev: @sdiv; */
  /* $prereq := $0; $clr := $0;  $# reset */

#ifdef SPB_MPI
  /* local data accumulator, for MPI communication
   * $io := ; $clr := 1;
   * $mpi := alloc; $sync := 1; $redtmp := @dbltmp; */
  double *hist_l; /* $reduce: @hist; */
  double *sbf_l;  /* $reduce: @sbf; */
  double *sdiv_l; /* $reduce: @sdiv; */
  double *sbl_l;  /* $reduce: @sbl; */
  /* $cnt := $0; $io := $0; $clr := $0;
   * $mpi := $0; $sync := $0; $redtmp := $0; $# reset  */
#endif


  /* $clr := 1; */
  double *dbltmp;  /* for temporary use, i.e., synchronizing, io
                      $cnt = @bins + 1; $io:;
                      $prereq = $ismaster; */
  double *mf;      /* array of mean force $cnt = @bins;
                      $io:bt; $binprev: @sbl; $mpi; $bcast;
                      $rbvalid: @-calcmfcom(@) == 0 && @-calcpmf(@) == 0;
                      */
  double *pmf;     /* array of p.m.f $io:none;
                      $cnt = @bins + 1; $mpi; */
  /* $clr := $0; */

  double mfcom;   /* average force $io:none; */
  double pmfmin;  /* minimal of pmf, $io:none; */
  double pmfmax;  /* maximal of pmf, $io:none; */

  /* buffer list (temporary atom indices, force, etc)
   * $io := ; $# turn off io */
  int bufcnt; /* number of items in buf */
  int bufcap; /* allocated memory */
  spbbuf_t *buf; /* $cnt=0; $# don't allocation space at the beginning*/

  /* GROMACS extensions, they are here to avoid another struct */
  int type;     /* the interaction *type* in gromacs interaction list */
  int exclcnt;  /* # of items */
  int exclcap;  /* capacity of excl[] */
  int *excl;    /* interactions (atom indices) that are considered as bad
                   $cnt=0; $mpicnt = @exclcap * @natoms; $mpi; */
  /* $reducevalid: !($ismaster) || @-calcmf(@) == 0;
   *               master computes mean force */
  /* $bcastvalid: @-calcmfcom(@) == 0 && @-calcpmf(@) == 0;
   *              every node computes mfcom and pmf from mf*/
} spb_t; /* $private */

typedef struct {
  /* $kprefix := spb_; */
  const char  *iopfx;    /* $usr:cfgref; */
  char  *bin_file;  /* name of the binary file for IO (master)
                       $def = $> @@ = ssdup("")\;\nif (@iopfx != NULL) sscpy(@@, @iopfx)\;\nsscat(@@, "spb.bin");
                       $prereq: $ismaster; $cfgprereq:1; */
  char  *txt_file;  /* name of the text file for IO (master)
                       $def = $> @@ = ssdup("")\;\nif (@iopfx != NULL) sscpy(@@, @iopfx)\;\nsscat(@@, "spb.txt");
                       $prereq: $ismaster; $cfgprereq:1; */

  /* default global settings that applies to every spb */
  double barred;     /* factor to reduce barrier, used in md2conv;
                        $def: 1.0; $valid: @@ >= 1.0; */
  double distmin;     /* minimal relative value for spb distribution
                         $def=1e-4; */
  double distmax;     /* maximal relative value for spb distribution
                         $def=1.0; $valid: @distmin < @distmax; */
  double distpwrscl;  /* spb scaling factor to adjust the input distribution */
  double mfmax;       /* spb maximal magnitude */
  double popscal[SPB_PMFSCAL_CNT]; /* $def = 1.0; $io:c; */
  int pmfside; /* method of determine the baseline of spb pmf
                  $def = -1;  */

  int formula; /* use cosine-sine series $def = 0; */

  int    cnt;       /* the number of special bonds (input)
                       $io:cbt; $rbverify; */
  spb_t *arr;       /* array of individual spbs
                       $obj; $cnt: @cnt;
                       $cfgargs: @, i;
                       $rbargs: @rflags;
                       $clr; $mpi; $reduce; $bcast; */
  unsigned rflags; /* $usr:rb; */
} spbonds_t; /* $ptrname: bs; $fprefix: spbs_; */

#define SPB_CHECKNULL     0x10000  /* not null */
#define SPB_CHECKMASTER   0x20000  /* must be master */
#define spbs_check(bs, flags) spbs_check_(bs, flags, __FILE__, __LINE__)
static int spbs_check_(const spbonds_t *bs, unsigned int flags, const char *file, int line)
{
  if (flags & SPB_CHECKNULL && bs == NULL) {
    fprintf(stderr, "null pointer to spbonds_t at file %s, line %d\n",
            file, line);
    exit(1);
  }
  if (flags & SPB_CHECKMASTER && bs->mpi_rank != 0) {
    fprintf(stderr, "#%d: not master, at file %s, line %d\n",
            bs->mpi_rank, file, line);
    exit(1);
  }
  return 0;
}

/* return bin id for a given ang */
static int spb_x2i(spb_t *spb, double x)
{
  int i = (int)((x - spb->min) / spb->binw);
  /* fix a wrong bin index due to round-off error */
  if (i == -1 && (spb->min - x) < 0.01*spb->binw) {
    i = 0;
  } else if (i == spb->bins && (x - spb->max) < 0.01*spb->binw) {
    i = spb->bins - 1;
  }
  die_if(i < 0 || i >= spb->bins, "bad index at %g (%d), min=%g, max=%g\n",
    x, i, spb->min, spb->max);
  return i;
}

/* add an entry, f: -dH/da, div: divergence
 * if the Hamiltonain H = H0 + H1, the function should be
 * called *after* the contribution from H1 is included */
static int spb_add(spb_t *spb, double x, double f, double divr,
    double beta, double lambda)
{
  int i;
  double *hist, *sbf, *sbl, *sdiv;

  i = spb_x2i(spb, x);
#ifdef SPB_MPI
  /* in MPI, we add to the local buffer, sync later */
  hist = spb->hist_l;
  sbf  = spb->sbf_l;
  sdiv = spb->sdiv_l;
  sbl  = spb->sbl_l;
#else
  hist = spb->hist;
  sbf  = spb->sbf;
  sdiv = spb->sdiv;
  sbl  = spb->sbl;
#endif
  hist[i] += 1.0;
  sdiv[i] += divr;
  sbf[i]  += beta * f;
  sbl[i]  += beta * lambda;
  return 0;
}

/* compute mean force of bin i */
static int spb_calcmf(spb_t *spb)
{
  double num, den, mf;
  int i;

  for (i = 0; i < spb->bins; i++) {
    num = spb->mfref[i]*spb->hist[i] - (spb->sdiv[i] + spb->sbf[i]);
    den = spb->sbl[i];
    mf = (fabs(den) > 0.001) ? (num/den) : 0.0;
    if (mf > spb->mfmax)
      mf = spb->mfmax;
    else if (mf < -spb->mfmax)
      mf = -spb->mfmax;
    spb->mf[i] = mf;
  }
  return 0;
}

/* for a periodic spb, compute the average force component,
 * to be subtracted to avoid discontinuity at the boundaries
 * mf of an unvisited bin should be zero */
static int spb_calcmfcom(spb_t *spb)
{
  double mfcom;
  int i;

  if (spb->flags & SPB_PERIODIC) {
    for (mfcom = 0.0, i = 0; i < spb->bins; i++)
      mfcom += spb->mf[i];
    mfcom /= spb->bins;
  } else mfcom = 0.0;
  spb->mfcom = mfcom;
  return 0;
}

/* compute potential of mean force at bin boundaries spb->pmf[]
 * also compute spb->pmfmin and spb->pmfmax */
static int spb_calcpmf(spb_t *spb)
{
  double emax, emin, ene;
  int binid;

  /* loop over every bin to get the maximum/minimum */
  emax = emin = spb->pmf[0] = ene = 0.0;
  for (binid = 0; binid < spb->bins; binid++) {
    ene -= (spb->mf[binid] - spb->mfcom) * spb->binw;
    spb->pmf[binid + 1] = ene;
    if (ene < emin)
      emin = ene;
    else if (ene > emax)
      emax = ene;
  }
  /* avoid numerical error due to round off */
  if (spb->flags & SPB_PERIODIC)
    spb->pmf[spb->bins] = spb->pmf[0];
  spb->pmfmin = emin;
  spb->pmfmax = emax;
  return 0;
}

/* compute the reference mean force from distref */
int spb_calcmfref(spb_t *spb)
{
  int i;

  for (i = 0; i < spb->bins; i++)
    spb->mfref[i] = log(spb->distref[i+1]/spb->distref[i]) / spb->binw;
  return 0;
}

/* calculate the mean force and its potential of all bins
 * set flags to a combination of SPB_MFCOM,SPB_PMF, SPB_MF
 * for periodic spb (SPB_PERIODIC), mf = spb->mf[] - spb->mfcom */
static void spbs_calcmfpmf(spbonds_t *bs, unsigned flags)
{
  int id;
  spb_t *spb;

  spbs_check(bs, SPB_CHECKNULL);
  if (bs->cnt <= 0) return;

  for (id = 0; id < bs->cnt; id++) {
    spb = bs->arr+id;
    if (flags & SPB_MF) spb_calcmf(spb); /* mean force */
    if (flags & SPB_MFCOM) spb_calcmfcom(spb); /* average mean force */
    if (flags & SPB_PMF) spb_calcpmf(spb); /* pmf */
  }
}

/* compute potential of mean force, interpolate within bin, subtract baseline */
static double spb_pmflin(const spb_t *spb, double x, int i)
{
  double v, dx;

  v  = spb->pmf[i];
  dx = x - (spb->min + i*spb->binw);
  v += -(spb->mf[i] - spb->mfcom)*dx; /* last bin correction */
  if (spb->pmfside != 0) /* subtract the base line */
    v -= (spb->pmfside > 0) ? spb->pmfmin : spb->pmfmax;
  return v;
}

/* compute pmf from a cosine-sine series */
static double spb_pmfcs(spb_t *spb, double x, double *mf)
{
  double v, dv, ck, cx, cx2, sx, a, b;
  int k;

  die_if (!spb->formula, "must use formula mode");

  cx = cos(x);
  cx2 = cx * cx;
  sx = sin(x);
  /* v = a0 + b0*sin(x), dv = b0*cos(x); */
  a = spb->paramc[0];
  b = spb->params[0];
  v = a + b * sx;
  dv = b * cx;
  ck = 1.f;
  for (k = 1; k < SPB_CSORDER; k++) {
    /* cos part: v = ak cos(x)^k, dv = - ak k cos(x)^(k-1) sin(x)
     * sin part: v = bk cos(x)^k sin(x), dv = bk cos(x)^(k-1) [-k + (k+1)*cos(x)^2] */
    a = spb->paramc[k];
    b = spb->params[k];
    dv -= ck*(a*k*sx + b*(k - (k+1)*cx2));
    ck *= cx;
    v += ck*(a + b*sx);
  }
  *mf = -dv;
  return v;
}

/* get the calculated mf & pmf for a single spb,
 * spb->mf[], spb->mfcom, spb->pmf[] must exist! */
SPBSTRCLS double spb_getpmf(spb_t *spb, double x, double *mf)
{
  if (spb->formula) { /* use cosine-sine series */
    return spb_pmfcs(spb, x, mf);
  } else {
    int i = spb_x2i(spb, x);
    if (mf != NULL) *mf = spb->mf[i] - spb->mfcom;
    return spb_pmflin(spb, x, i);
  }
}

/* round spb->min and spb->max to multiples of M_PI, if possible */
static void spb_round2pi(spb_t *spb)
{
  const double tol = 1e-2;

  if (fabs(spb->min + M_PI) < tol) {
    spb->min = -M_PI;
    fprintf(stderr, "spb->min is fixed to %+.10f\n", spb->min);
  }
  if (fabs(spb->max - M_PI) < tol) {
    spb->max = +M_PI;
    fprintf(stderr, "spb->max is fixed to %+.10f\n", spb->max);
  }
  spb->binw = (spb->max-spb->min)/spb->bins;
}


/* parse cfg. string of atom names, separated by space, colon, comma, or dash
 * one character per delimiter, every delimiter is replaced by '\0'
 * assign spb->atoms, spb->atmsiz, spb->natoms */
static int spb_parse_atoms_txt(spb_t *spb, int verbose)
{
  int j;
  char *p;
  const char *delim = " :-,"; /* space, coma, minus and : can be used as delimiters */

  p = spb->atmbuf;
  spb->atmsiz = (p != NULL) ? (strlen(p) + 1) : 0;
  for (j = 0; j < SPB_MAXNATOMS && p != NULL; j++) {
    spb->atoms[j] = p;
    p = strpbrk(p, delim);
    if (p != NULL)
      *p++ = '\0';
    if (verbose)
      fprintf(stderr, "%d:%s ", j, spb->atoms[j]);
  }
  if (verbose)
    fprintf(stderr, "\n");
  spb->natoms = j;
  return 0;
}

/* parse a string of atom names separated by null characters
 * the string length is specified by spb->atmsiz,
 * the number of items is limited by spb->natoms
 * assign spb->atoms */
static int spb_parse_atoms_bin(spb_t *spb, int verbose)
{
  int j;
  char *p, *q;

  p = spb->atmbuf;
  q = p + spb->atmsiz;
  for (j = 0; j < spb->natoms && p < q; j++) {
    spb->atoms[j] = p;
    p += strlen(p)+1;
    if (verbose) fprintf(stderr, "%d:%s ", j, spb->atoms[j]);
  }
  if (verbose) fprintf(stderr, "\n");
  return 0;
}

/* parse the string input followed "spbN_distref = "
 * in configuration file, it should look like
 *
 *   bins min max | val_1 val_2 ... val_bins
 *
 * where val_1 ... val_bins specified the full distribution.
 * for a uniform distribution, this can be omitted, thus, just
 *
 *   bins min max
 *
 * If min = -pi and max = pi, they can be omitted as well, so just
 *
 *   bins
 *
 * spb->bins, sbp->min, spb->max and spb->distref will be set.
 * other properties are initialized by spb_initdefault()
 * */
static int spb_parse_sdistref(spb_t *spb, char *input)
{
  int i, npmf, periodic = spb->flags & SPB_PERIODIC;
  char *p, *q;
  double x, y;

  i = (int) strtol(input, &p, 10);
  die_if (i <= 0 || i > 100000, "bad # of bins %d for spb %d [%s]\n",
          i, spb->id, input);
  spb->bins = i;
  /* try to get the range (if it exists) */
  spb->min = -M_PI;
  spb->max =  M_PI;
  if (2 == sscanf(p, "%lf%lf", &x, &y)) {
    spb->min = x;
    spb->max = y;
  }
  spb_round2pi(spb); /* fix spb->min/max to multiples of PI, if possible */
  spb->binw = (spb->max-spb->min)/spb->bins;
  die_if ((spb->distref = calloc(spb->bins + 1, sizeof(double))) == NULL,
      "no memory for distref");
  die_if ((spb->distref0 = calloc(spb->bins + 1, sizeof(double))) == NULL,
      "no memory for distref0");
  /* now parse the distribution line, it looks like
   * spb1_dist = bins min max | val_1 val_2 ... val_bins */
  if ((q = strchr(input, '|')) == NULL) goto NO_DISTREF;
  q++; /* q now points to the start of distribution data */
  /* skip spaces and unnecessary delimiters
   * q is the start of the distribution (just after `|') */
  for (; *q && strchr(";,| \t\n\r", *q) != NULL; ) q++;
  if (*q == '\0') goto NO_DISTREF;

  npmf = spb->bins +  (periodic ? 0 : 1);
  /* starting from the character after `|', read the distribution */
  for (p = q, i = 0; i < npmf; ) {
    x = strtod(p, &q);
    spb->distref[i++] = x;
    if (*q == '\0') break; /* no further entry */
    p = q;
  }
  if (i < npmf) {
    fprintf(stderr, "Warning: insufficient items for distribution of spb %d,"
        " after %d/%d, zeros assumed.\n", spb->id, i, npmf);
    for (; i < npmf; i++) spb->distref[i] = 0.0;
  }
  if (periodic) spb->distref[spb->bins] = spb->distref[0];
  for (i = 0; i <= spb->bins; i++)
    spb->distref0[i] = spb->distref[i];
  printf("\n");
  for (i = 0; i <= spb->bins; i++) printf("%d: %g, ", i, spb->distref[i]);
  printf("\n");
  return 0;
NO_DISTREF:
  for (i = 0; i <= spb->bins; i++) {
    spb->distref0[i] = spb->distref[i] = 1.0;
  }
  fprintf(stderr, "No ref. dist. for spb %d, buf=[%s]\n"
        "  assume a uniform one\n", spb->id, input);
  return 0;
}

/* tune the reference distribution (if it's not uniform) */
static int spb_tunedistref(spb_t *spb)
{
  const double tol = 1e-5;
  double x;
  int i, verbose = 1, isuniform;

  /* 0. check to see if it is a uniform distribution */
  isuniform = (spb->distref[0] > 1e-9);
  if (isuniform) {
    for (i = 1; i <= spb->bins; i++)
      if (fabs(spb->distref[i] - spb->distref[0]) > tol*spb->distref[0]) {
        isuniform = 0;
        break;
      }
    if (isuniform) { /* indeed uniform */
      if (verbose)
        fprintf(stderr, "The target distribution is uniform at %g, "
            "no tuning necessary.\n", spb->distref[0]);
      return 0;
    }
  }
  /* 1. apply the potential scaling factor */
  if (fabs(spb->distpwrscl - 1.0) > 1e-6){
    for (i = 0; i <= spb->bins; i++) {
      x = spb->distref[i];
      spb->distref[i] = (x > 0.0) ? exp(log(x) * spb->distpwrscl) : 0.0;
    }
  }
  /* 2. get the maximal value */
  for (x = 0, i = 0; i <= spb->bins; i++)
    if (spb->distref[i] > x)
      x = spb->distref[i];
  if (x < 1e-9) { /* we do not allow an empty distribution */
    fprintf(stderr, "invalid distribution for spb %d\n", spb->id);
    return -1;
  }
  /* 3. normalize and truncate */
  for (i = 0; i <= spb->bins; i++) {
    spb->distref[i] /= x;
    if (spb->distref[i] < spb->distmin) {
      spb->distref[i] = spb->distmin;
    } else if (spb->distref[i] > spb->distmax) {
      spb->distref[i] = spb->distmax;
    }
    if (verbose) {
      fprintf(stderr, "%3d: %8.6f  ", i, spb->distref[i]);
      if ((i+1) % 5 == 0) fprintf(stderr, "\n");
    }
  }
  if (verbose) fprintf(stderr, "\n");
  /* wrap around for periodic spb */
  if (spb->flags & SPB_PERIODIC) {
    fprintf(stderr, "Enforcing pmf for spb %d at %.14f\n",
        spb->id, spb->distref[0]);
    spb->distref[ spb->bins ] = spb->distref[0];
  }
  return 0;
}

static int spb_initdistref(spb_t *spb, char *input)
{
  spb_parse_sdistref(spb, input);
  spb_tunedistref(spb);
  return 0;
}

/* write as text */
SPBSTRCLS int spbs_write(spbonds_t *bs, const char *fname, int ver, int pmf_asis)
{
  int j, id, econ = 1, nrows = 0;
  spb_t *spb0, *spb;
  FILE *fp;

  spbs_check(bs, SPB_CHECKNULL | SPB_CHECKMASTER);

  if (bs->cnt <= 0) return 0;

  if (!pmf_asis) spbs_calcmfpmf(bs, SPB_PMFALL);

  spb0   = bs->arr; /* use the first spb as a reference */
  nrows = spb0->bins; /* # of rows to print = maximal # of bins */
  for (id = 1; id < bs->cnt; id++) {
    spb = bs->arr + id;
    if ( fabs(spb->min - spb0->min) > 1e-8
      || fabs(spb->max - spb0->max) > 1e-8
      || spb->bins != spb0->bins) {
      econ = 0; /* cannot use economic (multiple-column) mode */
    }
    if (spb->bins > nrows) nrows = spb->bins;
  }
  nrows++;

  if (fname == NULL) fname = bs->txt_file;
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }

  /* specify version and if economic mode is used */
  fprintf(fp, "# %d %d %d | ", bs->cnt, ver, econ);

  for (id = 0; id < bs->cnt; id++) {
    spb = bs->arr + id;
    fprintf(fp, "%d %20.14e %20.14e %20.14e %d ",
        spb->bins, spb->min, spb->max, spb->mfmax, spb->pmfside);
    fprintf(fp, "%d ", spb->natoms);
    for (j = 0; j < spb->natoms; j++) {
      if (j > 0) fprintf(fp, "-");
      fprintf(fp, "%s", spb->atoms[j]);
    }
    fprintf(fp, " | ");
  }
  fprintf(fp, "\n");

  for (j = 0; j < nrows; j++) {
    for (id = 0; id < bs->cnt; id++) {
      /* a quick initialization for default values */
      double arr[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

      spb = bs->arr + id;
      arr[0] = spb->min + j*spb->binw;
      if (j < spb->bins) {
        arr[1] = spb->hist[j];
        if (fabs(spb->sbl[j]) > 0.0) {
          arr[2] = spb->sbf[j] / spb->sbl[j];
          arr[3] = spb->sdiv[j] / spb->sbl[j];
        }
        if (spb->hist[j] > 0.0)
          arr[4] = spb->sbl[j] / spb->hist[j];
        arr[7] = spb->mfref[j];
      }
      if (j <= spb->bins) {
        arr[5] = spb->pmf[j];
        arr[6] = spb->distref[j];
        arr[8] = spb->distref0[j];
      }
      if (id == 0 || !econ)
        fprintf(fp, "%8.5f ", arr[0]);
      /* if j >= spb->bins, zeroes will be written */
      fprintf(fp, "%13.3f %24.17e %24.17e %20.16f %10.6f %10.6f %10.6f ",
              arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7]);
      if (ver >= 2)
        fprintf(fp, "%10.6f ", arr[8]);
    } /* loop over spbs */
    fprintf(fp, "\n");
  } /* loop over rows */
  fclose(fp);
  return 0;
}

/* GROMACS extension, add an atom-set to the exclusion list */
SPBSTRCLS int spb_excl(spb_t *spb, const int *ids, int offset)
{
  const int bufsiz = 16;
  int j;

  if (spb == NULL) return -1;
  if (spb->exclcnt >= spb->exclcap) {
    spb->exclcap += bufsiz;
    spb->excl = realloc(spb->excl, spb->exclcap*sizeof(*ids)*spb->natoms);
    die_if (spb->excl == NULL, "no memory for excl list %d, size: %d\n",
        spb->id, spb->exclcap);
  }
  for (j = 0; j < spb->natoms; j++)
    spb->excl[spb->exclcnt * spb->natoms + j] = ids[j] + offset;
  spb->exclcnt++;
  return 0;
}

/* add spb atom indices, grad, ... to buffer
 * cannot use spb_add now, because force is not yet available */
SPBSTRCLS int spb_buf(spb_t *spb, const int ia[], real grad[][3],
    double x, double divr, double grad2, double mfph)
{
  const int blocksiz = 16;
  int j;
  spbbuf_t *cc;

  if (spb == NULL) return -1;
  if (spb->bufcnt >= spb->bufcap) {
    spb->bufcap += blocksiz;
    /* expand the current buffer, if it is out of capacity */
    spb->buf = realloc(spb->buf, spb->bufcap * sizeof(spb->buf[0]));
    die_if (spb->buf == NULL,
      "no memory for buf, %d, size=%d\n", spb->id, spb->bufcap);
  }
  die_if (spb->natoms > SPB_MAXNATOMS,
    "too many atoms spb %d, n: %d\n", spb->id, spb->natoms);
  cc = spb->buf + spb->bufcnt;
  cc->x = x;
  cc->div = divr;
  cc->grad2 = grad2;
  cc->mfph = mfph;
  for (j = 0; j < spb->natoms; j++) {
    cc->idx[j] = ia[j];
    rv3_copy(cc->grad[j], grad[j]);
  }
  spb->bufcnt++;
  return 0;
}

/* add bufferred items */
SPBSTRCLS int spbs_debuf(spbonds_t *bs, real f[][3], double beta, double lambda)
{
  int id, ic, j;
  spb_t *spb;
  spbbuf_t *cc;
  double fv, g2, fvself;

  spbs_check(bs, SPB_CHECKNULL);
  if (bs->cnt <= 0) return 0;

  for (id = 0; id < bs->cnt; id++) {
    spb = bs->arr + id;
    if (spb->bufcnt <= 0) continue;
    for (ic = 0; ic < spb->bufcnt; ic++) {
      /* compute div + beta * grad(H).grad(x) and beta * lambda */
      cc = spb->buf + ic;
      g2 = cc->grad2;
      if (g2 <= 1e-20) continue;
      fv = 0.0;
      for (j = 0; j < spb->natoms; j++) {
        int ia = cc->idx[j];
        fv += rv3_dot(f[ia], cc->grad[j]);
      }
      fv /= g2;
      fvself = cc->mfph * lambda;
      fv -= fvself;
      spb_add(spb, cc->x, fv, cc->div, beta, lambda);
    }
    spb->bufcnt = 0; /* clear buffer */
  }
  return 0;
}

/* collects histogram-like data hist/sdiv/sbf/sbl, master calculates the mf,
 * bcast mf, everyone computes the pmf, mfcom */
SPBSTRCLS int spbs_syncdata(spbonds_t *bs)
{
#ifndef SPB_MPI
  /* single processor case */
  spbs_calcmfpmf(bs, SPB_PMFALL);
#else
  if (bs->cnt <= 0) return 0;
  /* the MPI version */
  if (0 != spbs_reduce(bs)) { /* collect local data */
    fprintf(stderr, "error during spbs_reduce: %d\n", bs->mpi_rank);
    return -1;
  }
  if (0 != spbs_bcast(bs)) {  /* bcast sync'ed result */
    fprintf(stderr, "error during spbs_bcast: %d\n", bs->mpi_rank);
    return -1;
  }
#endif /* SPB_MPI */
  return 0;
}
#endif /* SPB_C__ */


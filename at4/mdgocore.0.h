/* Go model
 * Copyright (c) 2011 Cheng Zhang */
#define ZCHAVEREAL 1
#define ZCOM_PICK
#define ZCOM_ENDN
#define ZCOM_RNG
#define ZCOM_LOG
#define ZCOM_CFG
#define ZCOM_PDB
#define ZCOM_ROTFIT
#define ZCOM_AV
#define ZCOM_TMH
#include "zcom.h"

#ifndef BOLTZ
#define BOLTZ  8.314511212e-3 /* Boltzmann constant */
#endif

/* define a long long int type */
#ifndef I32ONLY
typedef long long  llong_t;
#define llfmt "%10lld"
#else
typedef long llong_t;
#define llfmt "%10ld"
#endif

/* for object-generator */
#define USE_MPI GMX_MPI

/* simple neighbor list for an atom */
#define ATNBMAX 255
typedef struct {
  int id;    /* index of this atom */
  int nbcnt; /* number of neighbors */
  int nb[ATNBMAX];   /* index of neighboring atoms; */
} atnbls_t; /* $skip; $mpi:0; */

typedef struct {
  double tp0;         /* copy from GROMACS thermostat temperature;
                         $def: 300.0; $usr:cfg; $io:; */
  double beta0;       /* $def: 0.405; $io:none; */
  double tmstep;      /* copy from GROMACS MD integration step, for convenience;
                         $def: 0.001; $usr:cfg; $io:; */
  int    mode;        /* mode; $usr:cfg; */
  double boltz;       /* Boltzmann constant $usr:cfg; */

  /* Go model parameter */
  char *fnref;        /* file name for the reference pdb; $key: refpdb; $def: "ref.pdb"; */
  float nbe;          /* energy unit for a contact; $def: 1.0f; $io:c; */
  float rc;           /* cutoff of defining contacts; $def: 5.0f; $io:c; */
  int nat;            /* # of atoms */
  int nres;           /* # of residues */
  rv3_t *xref;        /* reference structure;
                         $io:none; $cnt:0; $# to avoid allocation during init.;
                         $mpi; $mpicnt:@nat; */
  rv3_t *fgo;         /* Go force $io:none; $cnt:0; $mpi; $mpicnt:@nat; */
  int   *isct;        /* nat X nat array of contact;
                         $io:none; $cnt:0; $mpi; $mpicnt:@nat*@nat; */
  real *drref;        /* nat X nat of reference distances;
                         $io:none; $cnt:0; $mpi; $mpicnt:@nat*@nat; */
  real *wt;           /* weight of heavy atom $io:none; $cnt:0; $mpi; $mpicnt: @nat; */

  int atcnt;          /* number of home atoms; $io:none; */
  atnbls_t *atnb;     /* neighbor list of home atoms $cnt:0; */

  float edih;         /* energy for each dihedral pair; $def: 1.0f;  $io:c; */
  float dihbias;      /* bias of signed dihedral; $io:c; $def: 0.0f; */
  int dihcnt;         /* number of dihedral restraints $io:0; */
  int *dihid;         /* dihedral indices; $cnt:0; $mpi; $mpicnt:(4 * @dihcnt); */
  real *dihref;       /* reference dihedral values; $cnt:0; $mpi; $mpicnt:@dihcnt; */
  real *dihdis;       /* 1-4 distance; $cnt:0; $mpi; $mpicnt:@dihcnt; */
  int dihlist;        /* how to make dihedral list $def:0;
                         0: heavy atoms; 2: backbone */
  float dihlistrc;    /* cutoff distance in making dihlist, $def: 10.0f; */
  char *helixrange;   /* string of helices index ranges,
                         e.g., "4-21, 28-45, 53-70";
                         $def: NULL; $io:c; */

  pdbmodel_t *pm;     /* pdbmodel */

  int nstgoadj;       /* frequency of adjusting lambda $io:c; $def:10; */
  int nstgorep;       /* frequency of reporting lambda $io:c; $def:100; */
  double H0;          /* original Hamiltonian, pot. energy; $def: 0.0;  $io:none; */
  double goep;        /* energy from Go model ; $def: 0.0;  $io:none; */
  double goepvdw, goepdih; /* vdw and dih energy; $io:none; */
  double goepdihbias; /* bias energy that is always there; $io:none; */
  double rmsd;        /* RMSD from the reference structure; $io:none; */
  double lambda;      /* current $io:none; */

  /* parameters */
  double lambda_init; /* initial lambda $io:c; */
  double lambda_min; /* smallest (least constrained) lambda   $io:c; $def:-100.0; */
  double lambda_max; /* largest (most constrained) lambda   $io:c; $def:100.0; */
  double dlambda; /* bin size: $io:c; $def: -0.01; */
  av_t avlamb[1]; /* average lambda for a H1-feedback simulation $io:none; */
  av_t avrmsd[1]; /* average rmsd for a constant-lambda simulation $io:none; */
  av_t avgoep[1]; /* average epot for a constant-lambda simulation $io:none; */
  double goeptarget; /* mode 1, target go energy $def:0; $io:c; */
  double lambdadt; /* mode 1 magnitude of adjusting lambda, $def:1e-3; $io:c; */

  int nstcheckf; /* frequency of checking forces $io:c; $def: 1000; */

  /* TMH stuff: lambda */
  tmh_t *tmh;
  double tmh_tp0;     /* $def: 0.10; */
  double tmh_tp1;     /* $def: 0.20; */
  double tmh_dtp;     /* $def: 0.0001; */
  double tmh_tpc;     /* $def: 0.10; */
  double tmh_epot0;   /* $def: -160; */
  double tmh_epot1;   /* $def:  -60; */
  double tmh_depot;   /* $def:    5; */
  double tmh_epotmin; /* $def: -200; */
  double tmh_epotmax; /* $def:    0; */
  double tmh_epotdel; /* $def:   0.1; */

  double tmh_elimit;  /* $def: 1e9; */
  double tmh_springk; /* $def: 0.0; */

  char  *tmh_fndhde;  /* $def: "mdgo.e"; */
  char  *tmh_fnehis;  /* $def: "mdgo.ehis"; */
  char  *tmh_fntp;    /* $def: "mdgo.t"; */

  double tmh_ampmax; /* $def: 1e-4; */
  double tmh_ampmin; /* $def: 0.0; */
  double tmh_ampc;   /* $def: 1.0; */

  double tmh_ampfac;    /* $def: 0.316227766016838; */
  double tmh_percutoff; /* $def: 0.999999; */
  int    tmh_updampc; /* update ampc according to Wang-Landau; $def: 0; */

  int    tmh_entropic;  /* entropic sampling; $def: 0 */
  double tmh_entampmax; /* $def: 1e-4; */
  double tmh_entampmin; /* $def: 0.0; */
  double tmh_entampc;   /* $def: 100.0 */

  double tmh_lgvdt;  /* $def: 1e-5; */
  double tmh_dhde;   /* $io:none; $def: 1; */
  double tmh_dhdemin; /* $def: 0.1; */
  double tmh_dhdemax; /* $def: 10.0; */
  int    tmh_dhdeorder; /* $def: 1; */

  int tmh_nstsave;   /* $def: 1000000; */

  logfile_t *log;
  char *fnlog; /* $def: "mdgo.tr"; */
  char *fnrng; /* $def: "MTSEED"; */
} ago_t; /* $cfgopen; $rb:0; $wb:0; $reduce:0; */

/** returns if an atom is a backbone atom */
static int isbackbone(const char *atnm)
{
  if (strcmp(atnm, "CA") == 0) return 1;
  if (strcmp(atnm, "N")  == 0) return 1;
  if (strcmp(atnm, "C")  == 0) return 1;
  return 0;
}

/* make a list of constraint dihedral angles
 * based on backbone atoms (i, k1, k2, j)
 * For every pair r(i, j) < dihlistrc
 * choose k1 and k2, such that
 * 1. r(k1, k2) > dihlistrc
 * 2. r(i, k1) + r(j, k2) is minimum.
 * swap k1 and k2 to ensure k1 < k2  */
static void ago_mkdih14(ago_t *ago, const int *se, int bbonly)
{
  int i, j, k1, k2, k1s, k2s, ir, jr, kr1, kr2, did, nat = ago->nat;
  real rik1, rk2j, r2, rmin;
  double rc = ago->dihlistrc * 0.1; /* translate into nm */
  pdbmodel_t *pm = ago->pm;

  for (i = 0; i < nat; i++) {
    if (ago->wt[i] <= 0.0) continue;
    if (bbonly && !isbackbone(pm->atm[i].atnm)) continue;
    ir = pm->atm[i].rid;
    if (ir < se[0] || ir >= se[1]) continue;

    for (j = i+1; j < nat; j++) {
      if (ago->wt[j] <= 0) continue;
      /* make sure r(i, j) < rc */
      if (rv3_dist(ago->xref[i], ago->xref[j]) > rc) continue;
      if (bbonly && !isbackbone(pm->atm[j].atnm)) continue;
      jr = pm->atm[j].rid;
      if (jr < se[4] || jr >= se[5]) continue;

      /* searching for a non-contacting k1s, k2s such that
       * r(i, k1s) + r(k2s, j) is minimum */
      rmin = 1e9;
      k1s = k2s = -1;
      for (k1 = i+1; k1 < j; k1++) {
        if (ago->wt[k1] <= 0) continue;
        if (bbonly && !isbackbone(pm->atm[k1].atnm)) continue;
        kr1 = pm->atm[k1].rid;
        if (kr1 < se[2] || kr1 >= se[3]) continue;
        rik1 = rv3_dist(ago->xref[i], ago->xref[k1]);

        for (k2 = i+1; k2 < j; k2++) {
          if (k2 == k1 || ago->wt[k2] <= 0) continue;
          if (bbonly && !isbackbone(pm->atm[k2].atnm)) continue;
          if (rv3_dist(ago->xref[k1], ago->xref[k2]) <= rc) continue; /* k1 and k2 should be apart */
          kr2 = pm->atm[k2].rid;
          if (kr2 < se[2] || kr2 >= se[3]) continue;
          rk2j = rv3_dist(ago->xref[j], ago->xref[k2]);

          r2 = rik1 + rk2j;
          if (r2 > rmin) continue;
          rmin = r2;
          k1s = k1;
          k2s = k2;
        }
      }

      die_if (k1s < 0, "cannot find pivot for atoms %d and %d\n", i+1, j+1);
      did = ago->dihcnt;
      ago->dihcnt = did + 1;
      xrenew(ago->dihid, ago->dihcnt*4);
      xrenew(ago->dihref, ago->dihcnt);
      xrenew(ago->dihdis, ago->dihcnt);
      if (k1s > k2s) k1 = k1s, k1s = k2s, k2s = k1;
      ago->dihid[did*4 + 0] = i;
      ago->dihid[did*4 + 1] = k1s;
      ago->dihid[did*4 + 2] = k2s;
      ago->dihid[did*4 + 3] = j;
      ago->dihref[did] = rv3_calcdihv(NULL, ago->xref, ago->dihid + 4*did, DIH_POLYMER);
      ago->dihdis[did] = rv3_dist(ago->xref[i], ago->xref[j]);
    }
  }
}

/* make a list of constraint dihedral angles */
static void ago_mkdihhvy(ago_t *ago, const int *se)
{
  int i, j, k1, k2, ir, jr, kr1, kr2;
  pdbmodel_t *pm = ago->pm;

  for (k1 = 0; k1 < ago->nat; k1++) {
    if (ago->wt[k1] <= 0.0) continue;
    kr1 = pm->atm[k1].rid;
    if (kr1 < se[2] || kr1 >= se[3]) continue;

    for (k2 = k1 + 1; k2 < ago->nat; k2++) {
      if (ago->wt[k2] <= 0.0) continue;
      if (!ago->isct[k1*ago->nat + k2]) continue;
      kr2 = pm->atm[k2].rid;
      if (kr2 < se[2] || kr2 >= se[3]) continue;

      for (i = 0; i < k1; i++) {
        if (ago->wt[i] <= 0.0) continue;
        if (!ago->isct[i*ago->nat + k1]) continue;
        ir = pm->atm[i].rid;
        if (ir < se[0] || ir >= se[1]) continue;

        for (j = k2+1; j < ago->nat; j++) {
          int did;
          if (ago->wt[j] <= 0.0) continue;
          if (!ago->isct[k2*ago->nat + j]) continue;
          jr = pm->atm[j].rid;
          if (jr < se[4] || jr >= se[5]) continue;

          did = ago->dihcnt;
          ago->dihcnt = did + 1;
          xrenew(ago->dihid, ago->dihcnt*4);
          xrenew(ago->dihref, ago->dihcnt);
          xrenew(ago->dihdis, ago->dihcnt);
          ago->dihid[did*4 + 0] = i;
          ago->dihid[did*4 + 1] = k1;
          ago->dihid[did*4 + 2] = k2;
          ago->dihid[did*4 + 3] = j;
          ago->dihref[did] = rv3_calcdihv(NULL, ago->xref, ago->dihid + 4*did, DIH_POLYMER);
          ago->dihdis[did] = rv3_dist(ago->xref[i], ago->xref[j]);
          //printf("%d: %d - %d - %d - %d, %g\n", did+1, i+1, k1+1, k2+1, j+1, ago->dihref[did]);
        }
      }
    }
  }
}

/* parse the helix range string to */
static int ago_parsehelixrange_cfg(ago_t *ago, int **pse)
{
  char *s = ago->helixrange, *p, *q, *r;
  const int SIZE = 16;
  int ng = 0, ngcap = SIZE, *se, i, j;

  if (s == NULL) return 0;
  strip(s);
  if (s[0] == '\"') s++;
  p = s + strlen(s) - 1;
  if (*p == '\"') *p = '\0';
  *pse = NULL;
  xnew(se, 2*ngcap);
  for (p = s; ;) {
    q = strchr(p, ',');
    if (q) *q = '\0';
    strip(p); /* p is something like "4-20" */
    if ((r = strchr(p, '-')) == NULL) {
      fprintf(stderr, "failed parse item %s\n", p);
      free(se);
      return 0;
    }
    *r = '\0';
    i = atoi(p);
    j = atoi(r+1);
    se[2*ng] = i - 1;
    se[2*ng+1] = j - 1 + 1;
    if (++ng >= ngcap) {
      ngcap += SIZE;
      xrenew(se, 2*ngcap);
    }
    if (q == NULL) break;
    else p = q + 1;
  }
  *pse = se;
  return ng;
}

/* parse helix */
static int ago_parsehelixrange(ago_t *ago, pdbaac_t *pc, int **pse)
{
  int i, ng = ago_parsehelixrange_cfg(ago, pse);
  if (ng == 0) {
    fprintf(stderr, "use the default algorithm to parse helix\n");
    ng = pdbaac_parsehelices(pc, pse);
  }
  for (i = 0; i < ng; i++) {
    printf("Group %2d: %5d - %5d\n", i, (*pse)[2*i], (*pse)[2*i + 1]);
  }
  return ng;
}

/* load reference structure
 * make contact map */
static int ago_loadxref(ago_t *ago, const char *fn)
{
  int i, j, ir, jr, prid, prid2;
  int ng, *se, k1, k2, kr1, kr2;
  real rc, dr;
  pdbmodel_t *pm;
  pdbaac_t *pc;

  die_if ((pm = pdbm_read(fn, 0)) == NULL, "cannot readpdb %s\n", fn);
  ago->pm = pm;
  ago->nat = pm->natm;
  ago->nres = pm->nres;
  xnew(ago->xref, ago->nat);
  xnew(ago->fgo, ago->nat);
  xnew(ago->wt, ago->nat);

  for (i = 0; i < ago->nat; i++) rv3_zero(ago->fgo[i]);

  /* get xref */
  for (i = 0; i < ago->nat; i++) {
    /* copy from pm->atm (A) to ago->xref (nm), x0.1 */
    rv3_smul(rv3_copy(ago->xref[i], pm->atm[i].x), 0.1f);
    j = pm->atm[i].atnm[0];
    if (j == 'H' || j == 'D') ago->wt[i] = 0.0f;
    else if (j == 'N') ago->wt[i] = 14.0f;
    else if (j == 'O') ago->wt[i] = 16.0f;
    else if (j == 'S') ago->wt[i] = 32.0f;
    else ago->wt[i] = 12.0f;
  }

  /* computing reference distances */
  xnew(ago->drref, ago->nat * ago->nat);
  xnew(ago->isct, ago->nat * ago->nat);
  rc = ago->rc * .1f; /* to nm */
  for (i = 0; i < ago->nat; i++) {
    for (j = i + 1; j < ago->nat; j++) {
      prid  = i * ago->nat + j;
      prid2 = j * ago->nat + i;
      dr = rv3_dist(ago->xref[i], ago->xref[j]);
      ago->drref[prid2] = ago->drref[prid] = dr;
      ago->isct[prid2] = ago->isct[prid] = (dr < rc)
          && ago->wt[i] > 0. && ago->wt[j] > 0.;
    }
    prid = i * ago->nat + i;
    ago->drref[prid] = 0;
    ago->isct[prid] = 0;
  }

  /* compute dihedral sets */
  die_if ((pc = pdbaac_parse(pm, 0)) == NULL, "cannot parse protein chain, %s\n", fn);
  ng = ago_parsehelixrange(ago, pc, &se);
  ago->dihcnt = 0;
  xnew(ago->dihid, 4);
  xnew(ago->dihref, 1);
  xnew(ago->dihdis, 1);
  if (ng == 3) { /* build dihedrals */
    if (ago->dihlist == 1) {
      ago_mkdih14(ago, se, 1); /* 1-4 method backbone */
    } else if (ago->dihlist == 2) {
      ago_mkdih14(ago, se, 0); /* 1-4 heavy */
    } else {
      ago_mkdihhvy(ago, se);
    }
  }
  free(se);
  pdbaac_free(pc);
  //pdbm_free(pm);
  return 0;
}

/* calculate CA rmsd */
static real ago_carmsd(ago_t *ago, rv3_t x[])
{
  return rotfit3(x, NULL, ago->xref, ago->wt, ago->nat, NULL, NULL);
}

/* load previous mean force data */
static int ago_loaddata(ago_t *ago, int isctn)
{
  if (!isctn) return 0;
  /* load dhde data */
  if (ago->mode == 2) {
    double amp, t0;
    tmh_load(ago->tmh, ago->tmh_fnehis, ago->tmh_fndhde, &amp, &t0);
    if (ago->tmh->wl) {
      if (ago->tmh_ampc <= 0.0) { /* completely disable updating if ampc = 0.0 */
        ago->tmh->wl->lnfwl = 0;
        ago->tmh->wl->lnfc = 0;
        amp = 0;
      }
      ago->tmh->wl->lnf = amp;
    }
  }
  fprintf(stderr,"successfully load previous data\n");
  return 0;
}

/* basic initialization steps (master only)
 * 1. load configuration file
 * 2. load data
 * */
static ago_t *ago_open(const char *fname, int isctn, double tp0, double tmstep, int mode)
{
  ago_t *ago;

  /* this will also initialize settings for member objects such as ago->bbs */
  ago = ago_cfgopen((fname != NULL) ? fname : "mdgo.cfg", tp0, tmstep, mode, BOLTZ);
  die_if (ago == NULL, "failed to load configuration file %s\n", fname);
  ago->lambda = ago->lambda_init; /* initially uses no lambda */
  ago_loadxref(ago, ago->fnref);
  ago->beta0 = 1.0/(BOLTZ * ago->tp0); /* simulation temperature */
  ago->log = log_open(ago->fnlog);
  if (mode == 2) {
    double ensexp = 2;
    int dhdeorder = ago->tmh_dhdeorder;

    if (isctn) { /* load range from dhde file */
      tmh_loaderange(ago->tmh_fndhde, &ago->tmh_tp0, &ago->tmh_tp1, &ago->tmh_dtp,
          &ago->tmh_epot0, &ago->tmh_epot1, &ago->tmh_depot,
          &ago->tmh_epotmin, &ago->tmh_epotmax, &ago->tmh_epotdel,
          &ensexp, &dhdeorder);
    }
    die_if ((ago->tmh = tmh_open(
        ago->tmh_tp0, ago->tmh_tp1, ago->tmh_dtp,
        ago->tmh_epot0, ago->tmh_epot1, ago->tmh_depot,
        ago->tmh_epotmin, ago->tmh_epotmax, ago->tmh_epotdel,
        ensexp, dhdeorder)) == NULL, "cannot open tmh, dtp %g\n", ago->tmh_dtp);
    ago->tmh->elimit = ago->tmh_elimit;
    ago->tmh->springk = ago->tmh_springk;
    if (ago->tmh_entropic) {
      tmh_initwlcvg(ago->tmh, ago->tmh_entampc, ago->tmh_entampmax,
          ago->tmh_ampfac, ago->tmh_percutoff, ago->tmh_entampmin,
          ago->tmh_updampc ? WLCVG_UPDLNFC : 0);
    } else {
      tmh_initwlcvg(ago->tmh, ago->tmh_ampc, ago->tmh_ampmax,
          ago->tmh_ampfac, ago->tmh_percutoff, ago->tmh_ampmin,
          ago->tmh_updampc ? WLCVG_UPDLNFC : 0);
    }
    ago->tmh->scl = ago->beta0;
    ago->tmh->dhdemin = ago->tmh_dhdemin;
    ago->tmh->dhdemax = ago->tmh_dhdemax;
  }
  ago_loaddata(ago, isctn);
  return ago;
}

/* given dr2, return energy, *amp is the force amplitude */
static __inline real ago_pairene(real dr2, real drref, real *amp)
{
  real dr2ref = drref * drref, invr2, dr6, del, dr;

  if (dr2 > dr2ref) { /* long range: standard Lennard-Jones */
    invr2 = 1.0f/dr2;
    dr2 = dr2ref * invr2;
    dr6 = dr2 * dr2 * dr2;
    *amp = 12.f * (dr6 - 1.f) * dr6 * invr2;
    return (dr6 - 2.f) * dr6;
  } else { /* short range: harmonic */
    dr = (real) sqrt(dr2);
    del = 1.0f - dr/drref;
    *amp = 2.f * del/(drref * dr);
    return del * del - 1.f;
  }
}

/* given dih, dr2 (1-4), return energy
 * ene = -cos(dih - dihref) * f(r14)
 * f(r14) is given by ago_pairene(r14) */
static __inline real ago_dihene(real dih, real dihref, real dr2, real drref,
    real *ampdih, real *ampr)
{
  real rfac, dihdel, dfac;

  dihdel = dih - dihref;
  dfac = (real) cos(dihdel);
  rfac = ago_pairene(dr2, drref, ampr);
  *ampdih = (real) sin(dihdel) * rfac;
  *ampr *= dfac;
  return dfac*rfac;
}

/* compute the native contact energy
 * use for checking and a single processor */
static real ago_goepot(ago_t *ago, rv3_t x[])
{
  int i, j, id, nat = ago->nat, prid, nct = 0;
  real dx[3], dr2, nbe = ago->nbe, ene = 0;
  real dih, edih = ago->edih, amp, ampdih;

  /* van der Waales energy */
  for (i = 0; i < nat; i++) {
    for (j = i + 1; j < nat; j++) {
      prid = i * nat + j;
      if (!ago->isct[prid]) continue;
      dr2 = rv3_sqr( rv3_diff(dx, x[i], x[j]) );
      ene += nbe * ago_pairene(dr2, ago->drref[prid], &amp);
      nct++;
      //printf("%3d-%-3d: r = %g, %8.4f %6d\n", ir, jr, sqrt(dr2/dr2ref), ene, nct);
    }
  }

  /* contact dihedral energy assuming no pbc  */
  for (id = 0; id < ago->dihcnt; id++) {
    dih = rv3_calcdihv(NULL, x, ago->dihid + 4 * id, DIH_POLYMER);
    i = ago->dihid[4*id];
    j = ago->dihid[4*id+3];
    dr2 = rv3_sqr( rv3_diff(dx, x[i], x[j]) );
    ene += edih * ago_dihene(dih, ago->dihref[id], dr2, ago->dihdis[id], &ampdih, &amp);
  }
  return ene;
}

/* dump ago_t to file fn
 * for an array, only the first and last `arrmax' elements are dumped */
static int ago_dump(ago_t *ago, rv3_t *x, const char *fn, int arrmax)
{
  FILE *fp;
  int i, j, nat, pr, nct;
  real rij;

  xfopen(fp, fn, "w", return -1);
  ago_manifest(ago, fp, arrmax); /* call the standard manifest function */

  nat = ago->nat;

  /* print complete CA index */
  die_if (ago->pm == NULL, "%p PDB model\n", ago->pm);

  /* print contact map */
  for (nct = 0, i = 0; i < nat; i++) {
    for (j = i + 1; j < nat; j++) {
      pr = i * nat + j;
      if ( !ago->isct[pr] ) continue;
      rij = rv3_dist(x[i], x[j]);
      fprintf(fp, "pair %6d: %4d - %4d: %8.4f A, init %8.4f A\n",
          nct+1, i+1, j+1, ago->drref[pr]*10.0, rij*10.0);
      nct++;
    }
  }

  /* print dihedrals */
  fprintf(fp, "Contact dihedrals:\n");
  for (i = 0; i < ago->dihcnt; i++) {
    real dih = rv3_calcdihv(NULL, x, ago->dihid + 4*i, DIH_POLYMER);
    fprintf(fp, "dih %4d: %4d %4d %4d %4d %8.4f, init %8.4f\n", i+1,
        ago->dihid[4*i]+1, ago->dihid[4*i+1]+1, ago->dihid[4*i+2]+1, ago->dihid[4*i+3]+1,
        ago->dihref[i], dih);
  }

  fprintf(fp, "%d residues, %d atoms, rmsd %g A, init.goep %g, # of contacts %d, dihedrals %d\n",
      ago->nres, ago->nat, ago->rmsd*10.0, ago->goep, nct, ago->dihcnt);
  fclose(fp);
  return 0;
}



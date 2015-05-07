/* compute root-mean-square deviation distributions
 * Copyright (C) 2011 Cheng Zhang */
#include "zeutil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_RV3
#define ZCOM_HIST
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#include "zcom2.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fndist = "rms.dist";
const char *fndist2 = NULL;
const char *fnlog = NULL;
double tequil = 10000;
double Tref = 300; /* reference temperature */
double delbet = 1e3;
double xmax = 10.0;
double xdel = 0.001;
int ctmd = 0;

/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "compute root-mean-square deviation from the native";
  argopt_add(ao, "-s", NULL,  &fntop, "topology file");
  argopt_add(ao, "-f", NULL,  &fnxtc, "trajectory");
  argopt_add(ao, "-z", NULL,  &fnze,  "ZE file");
  argopt_add(ao, "-o", NULL,  &fndist, "output distribution");
  argopt_add(ao, "-a", NULL,  &fndist2, "second distribution");
  argopt_add(ao, "-g", NULL,  &fnlog, "log file");
  argopt_add(ao, "-q", "%lf", &tequil, "time of equilibration");
  argopt_add(ao, "-m", "%lf", &xmax,   "maximal rmsd");
  argopt_add(ao, "-d", "%lf", &xdel,   "bin for rmsd");
  argopt_add(ao, "-M", "%b", &ctmd, "constant temperature simulation");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

/* backbone information */
typedef struct {
  int ia[3];
  int resid;
  const char *resname;
} bb_t;
real *vmass; /* mass vector */

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
  xnew(vmass, (*rescnt + 3)*3);
  for (i = 0; i < *rescnt; i++) {
    vmass[3*i] = 14.01;
    vmass[3*i+1] = vmass[3*i+2] = 12.01;
  }
  return bb;
}

/* from full actom index to backbone index */
static void xgetbb(rv3_t *xbb, rv3_t *x, bb_t *bb, int rescnt)
{
  int i;

  for (i = 0; i < rescnt; i++) {
    rv3_copy(xbb[3*i],   x[ bb[i].ia[0] ]);
    rv3_copy(xbb[3*i+1], x[ bb[i].ia[1] ]);
    rv3_copy(xbb[3*i+2], x[ bb[i].ia[2] ]);
  }
}

/* scan trajectory and accumulate data */
static int dotrj(hist_t *hsrms, bb_t *bb, int rescnt, FILE *fplog,
    trj_t *trj, const char *fxtc, rvec *xrefbb,
    int nat, matrix box, double *wtot)
{
  int nfr, match, natoms;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL, *xbb;
  double t1, w, rmsd;
  real t = 0.f;

  xnew(xbb, rescnt * 3);
  /* loop over xtc frames */
  for (nfr = 0; ; nfr++) {
    if (nfr == 0) {
      /* first frame, get the file handle fpxtc, allocate x
       * and check the number of atoms */
      natoms = read_first_x(&fpxtc, fxtc, &t, &x, box);
      if (natoms != nat) {
        fprintf(stderr, "natoms mismatch! top %d, xtc %d\n", nat, natoms);
        return -1;
      }
    } else {
      if (!read_next_x(fpxtc, &t, natoms, x, box))
        break;
    }

    if (ctmd) { /* each frame has a weight */
      w = 1.;
    } else {
      /* for each xtc frame, we search a corresponding
       * frame in trj at the same time */
      for (match = 0; trj->pos < trj->n; trj->pos++) {
        t1 = trj->fr[trj->pos].t;
        if (fabs(t1 - t) < 0.5) {
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

    /* skip trivial frames */
    if (fplog == NULL && w < 1e-6) continue;

    /* extract backbone indices */
    xgetbb(xbb, x, bb, rescnt);

    /* fit x to xref */
    rmsd = rv3_rmsd(xbb, NULL, xrefbb, vmass, rescnt*3, 0, NULL, NULL);

    if (fplog) {
      fprintf(fplog, "%10.0f %8.4f %g %d\n", t, rmsd, w, trj->fr[trj->pos].dataid);
    }
    hs_add(hsrms, &rmsd, w, HIST_VERBOSE);
    *wtot += w;
  }
  free(xbb);
  return 0;
}

/* run simulation */
int run(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat, rescnt = 0;
  char fxtc[FILENAME_MAX] = "";
  char *ftop;
  matrix box;
  rvec *xref = NULL, *xrefbb = NULL;
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  double dt, bref = 1.0/(BOLTZ*Tref), wtot = 0.0;
  hist_t *hsrms = NULL;
  bb_t *bb = NULL;
  FILE *fplog = NULL;

  memset(box, 0, sizeof(box));

  /* read topology and the reference structure */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  printf("using the structure in %s as the reference\n", ftop);

  bb = init_index(&rescnt, mtop); /* get indices for backbone dihedrals */
  xnew(xrefbb, rescnt * 3);
  xgetbb(xrefbb, xref, bb, rescnt);

  hsrms = hs_open(1, 0.0, xmax, xdel);

  if (fnlog) {
    xfopen(fplog, fnlog, "w", exit(1));
  }

  if (fndist2 != NULL) {
    /* if there's a second database, combine it with the first database */
    die_if (0 != hs_load(hsrms, fndist, HIST_VERBOSE),
        "cannot load master %s\n", fndist);
    die_if (0 != hs_load(hsrms, fndist2, HIST_VERBOSE|HIST_ADDITION),
        "cannot add the second %s\n", fndist2);

  } else  {
    /* if there's no second database, we accumulate trajectories */

    /* 1. initialize multiple histogram weights for tempering */
    if (!ctmd) {
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
      /* for regular MD, set ndata to a large number */
      ndata = 10000000;
    }

    /* 2. we now successively load xtc from each data dir,
     * try to match the TRACE frames */
    if (!ctmd) trj->pos = 0;

    for (i = 1; i <= ndata; i++) {
      sprintf(fxtc, "data%d/%s", i, fnxtc);
      if (!fexists(fxtc)) { /* test if the file exists */
        die_if(!ctmd, "cannot read %s\n", fnxtc);
        break;
      }
      if (0 != dotrj(hsrms, bb, rescnt, fplog,
            trj, fxtc, xrefbb, nat, box, &wtot)) {
        fprintf(stderr, "error doing %s\n", fnxtc);
        return -1;
      }
    }
  }
  hs_save(hsrms, fndist, HIST_ADDAHALF|HIST_VERBOSE);
  hs_close(hsrms);
  sfree(mtop);
  if (!ctmd) {
    if (trj) trj_free(trj);
    if (ze) ze_free(ze);
  }
  free(xrefbb);
  return 0;
}

int main(int argc, char **argv)
{
  (void) argtp;
  doargs(argc, argv);
  run();
  return 0;
}


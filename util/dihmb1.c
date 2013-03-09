/* display dihedral information for an .xtc file
 * Copyright (C) 2011 Cheng Zhang */
#include "zeutil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_RV3
#define ZCOM_HIST
#define ZCOM_UTIL
#include "zcom.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnout = "dihp.txt";
const char *fndist = "dih.dist";
const char *fndist2 = NULL;
const char *fnreg = "dih.reg";
const char *fnreg2 = NULL;
double tequil = 10000;
double Tref = 300; /* reference temperature */
double delbet = 0.2;
int ctmd = 0;
const char *fnlog = NULL;
FILE *fplog = NULL;
const char *fnbracket = NULL; /* "fold.bra"; */

/* print usage and then quit */
static void usage(const char *prog)
{
  fprintf(stderr, "%s [OPTIONS] file.xtc\n"
      "It output phi/psi dihedral distributions, \n"
      "and faction of being different dihedral conformations,\n"
      "and a log file\n",
      prog);
  fprintf(stderr,
      " -s: followed by .tpr\n"
      " -f: followed by trajectory .xtc\n"
      " -z: followed by ZE file\n"
      " -t: followed by TRACE file\n"
      " -o: followed by output file (correct residue fraction)\n"
      " -D: followed by output distribution\n"
      " -a: followed by second distribution\n"
      " -g: followed by output log file\n"
      " -q: followed by time of equilibration\n"
      " -T: followed by reference temperature\n"
      " -d: followed by delta beta\n"
      " -b: followed by the bracket file\n"
      " -M: constant temperature MD\n"
      "\n");
  exit(1);
}

/* set parameters from command line arguments */
static int doargs(int argc, char **argv)
{
  int i, j, ch;
  const char *param;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnxtc = argv[i];
      continue;
    }
    ch = argv[i][1];
    if (strchr("fstzoDqdThgbarR", ch) != NULL) {
      param = argv[i] + 2;
      if (param[0] == '\0') {
        if (i >= argc - 1) {
          fprintf(stderr, "need arg. after %s\n", argv[i]);
          usage(argv[0]);
        }
        param = argv[++i];
      }
      if (ch == 'f') fnxtc = param;
      else if (ch == 's') fntop = param;
      else if (ch == 't') fntr = param;
      else if (ch == 'z') fnze = param;
      else if (ch == 'o') fnout = param;
      else if (ch == 'D') fndist = param;
      else if (ch == 'a') fndist2 = param;
      else if (ch == 'q') tequil = atof(param);
      else if (ch == 'd') delbet = atof(param);
      else if (ch == 'T') Tref = atof(param);
      else if (ch == 'g') fnlog = param;
      else if (ch == 'b') fnbracket = param;
      else if (ch == 'r') fnreg = param;
      else if (ch == 'R') fnreg2 = param;
      continue;
    }
    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'M': ctmd = 1; break;
      case 'h': usage(argv[0]); break;
      default:
        fprintf(stderr, "cannot handle %s\n", argv[i]);
        usage(argv[0]);
      }
  }
  if (fnlog != NULL) delbet = 1e3;
  /* check if it is constant temperature md */
  if (!ctmd && !fexists("data1/ZE")) {
    fprintf(stderr, "no ZE, assume constant temperature md\n");
    ctmd = 1;
  }
  return 0;
}

#define R2D (180/M_PI)
enum {ALPHA, BETA, LHELIX, LOOP, REGIONS};

typedef struct {
  int ia[5];
  int resid;
  const char *resnm;
  double phir, psir;
  int regr; /* region r */
  double hist[REGIONS];
  double smphi, smpsi, smphi2, smpsi2;
} bbdih_t;

/* note we use degree instead of rad */
#define XMAX 180
#define XMIN (-XMAX)
#define XDEL 1.0

/* traslate phi and psi to dihedral region integer */
static int dih2reg(double phi, double psi)
{
  if (phi < 0) {
    if (psi > -100 && psi < 80) return ALPHA;
    else return BETA;
  } else {
    if (psi > -80 && psi < 100) return LHELIX;
    else return LOOP;
  }
}

/* compute dihedrals in degree */
static int calcdih(rvec r[], int *idx, double *phi, double *psi)
{
  const double xm = 180-1e-12;
  *phi = R2D * (rv3_calcdihv(NULL, r, idx, 0));
  if (*phi < -xm) *phi = -xm; else if (*phi > xm) *phi = xm;
  *psi = R2D * (rv3_calcdihv(NULL, r, idx+1, 0));
  if (*psi < -xm) *psi = -xm; else if (*psi > xm) *psi = xm;
  return dih2reg(*phi, *psi);
}

/* add backbone index-set from moltype imt */
static bbdih_t *addbackbone(bbdih_t *bb, int *bbcnt, gmx_mtop_t *mtop, int imt, int ff[],
    int resid, const char *resnm, int verbose)
{
  int i, ag, im, imb, nmb = mtop->nmolblock, nbb;
  gmx_molblock_t *mb;

  nbb = *bbcnt;
  if (verbose)
    printf("%3d: adding %d, %d, %d, %d, %d\n",
      nbb+1, ff[0], ff[1], ff[2], ff[3], ff[4]);
  /* search molblock for those whose moltype is imt */
  for (ag = 0, imb = 0; imb < nmb; imb++) {
    mb = mtop->molblock + imb;
    if (mb->type == imt) {
      for (im = 0; im < mb->nmol; im++) {
        xrenew(bb, nbb+1);
        for (i = 0; i < 5; i++) {
          bb[nbb].ia[i] = ag + im*mb->natoms_mol + ff[i];
          bb[nbb].resid = resid;
          bb[nbb].resnm = resnm;
        }
        nbb++;
      }
    }
    ag += mb->nmol * mb->natoms_mol;
  }
  *bbcnt = nbb;
  return bb;
}

/* obtain backbone dihedral indices */
static bbdih_t *dihindex(int *bbcnt, gmx_mtop_t *mtop)
{
  static const char *bbatms[5] = {"C", "N", "CA", "C", "N"};
  gmx_moltype_t *mt;
  int start, ff[5], found;
  int imt, i, i2, resid, j;
  char *resnm;
  bbdih_t *bb = NULL;

  xnew(bb, 1);
  *bbcnt = 0;
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
      resnm = mt->atoms.resname[resid][0];
      bb = addbackbone(bb, bbcnt, mtop, imt, ff, resid, resnm, 0);
      i = ff[0] + 1;
    }
  } /* loop over moltypes */

/*
  for (i = 0; i < *bbcnt; i++) {
    for (j = 0; j < REGIONS; j++)
      bb[i].hist[j] = 0.0;
  }
*/
  return bb;
}

static __inline double DIHWRAP(double dx) {
  double pi = 180;
  return ((dx) > pi) ? (dx - 2*pi) : ((dx) < -pi) ? (dx + 2*pi) : (dx);
}

double (*bras)[2];
int nbras;

static int inbras(double t)
{
  int i;
  if (bras == NULL || nbras == 0) return 1; /* no brackets defined */
  for (i = 0; i < nbras; i++)
    if (t >= bras[i][0] && t < bras[i][1])
      return 1;
  return 0;
}

/* load brackets */
static int loadbras(const char *fn)
{
  int i;
  FILE *fp;
  static char buf[1024];

  nbras = 0;
  bras = NULL;
  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fn);
    return -1;
  }
  /* count the number of lines */
  for (nbras = 0; fgets(buf, sizeof buf, fp); ) {
    strip(buf);
    if (buf[0] == '\0') break;
    nbras++;
  }
  xnew(bras, nbras+1);
  rewind(fp);
  for (i = 0; i < nbras; i++) {
    fgets(buf, sizeof buf, fp);
    strip(buf);
    die_if (2 != sscanf(buf, "%lf%lf", &(bras[i][0]), &(bras[i][1])),
      "cannot scan bracket from line [%s]\n", buf);
  }
  fclose(fp);
  printf("brackets are loaded successfully\n");
  return 0;
}

/* scan trajectory and accumulate data */
static int dotrj(bbdih_t *bb, int bbcnt,
    hist_t *hs_fs, hist_t *hs_reg,
    trj_t *trj, const char *fxtc,
    int nat, matrix box, int *satmax, int *helmax, double *wtot)
{
  int nfr, i, match, reg, natoms, sat, hel;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL;
  double t1, w, phi, psi;
  real t = 0.f;

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
        /* when logging, we need a zero-w frames */
        if (fplog == NULL) continue;
      } else {
        w = exp(trj->fr[trj->pos].lnw);
      }
    }
    if (!inbras(t)) {
      continue;
    }

    /* compute dihedrals of backbone atoms, and the region where it belongs */
    for (hel = sat = 0, i = 0; i < bbcnt; i++) { /* add to the region histogram */
      reg = calcdih(x, bb[i].ia, &phi, &psi);
      //bb[i].hist[reg] += w;

      bb[i].smphi  += phi*w;
      bb[i].smphi2 += phi*phi*w;
      bb[i].smpsi  += psi*w;
      bb[i].smpsi2 += psi*psi*w;

      if (reg == bb[i].regr) sat++;
      if (reg == 0) hel++;

      hs_add1(hs_reg, 0, 1.0*reg + .5, w,  HIST_VERBOSE);
      hs_add1(hs_reg, i+1, 1.0*reg + .5, w,  HIST_VERBOSE);

      hs_add1(hs_fs, 0, phi, w, HIST_VERBOSE);
      hs_add1(hs_fs, 1, psi, w, HIST_VERBOSE);
      hs_add1(hs_fs, 2*i+2, phi, w, HIST_VERBOSE);
      hs_add1(hs_fs, 2*i+3, psi, w, HIST_VERBOSE);
    }
    *wtot += w;
    if (sat > *satmax) {
      *satmax = sat;
      //fprintf(stderr, "max sat %d: %s frame %d, t = %g\n", sat, fxtc, nfr+1, t);
    }
    if (hel > *helmax) {
      *helmax = hel;
      //fprintf(stderr, "max hel. %d: %s frame %d, t = %g\n", hel, fxtc, nfr+1, t);
    }
    if (fplog) {
      fprintf(fplog, "%.3f %g %g\n", t, 1.0*sat/bbcnt, 1.0*hel/bbcnt);
    }
  }
  return 0;
}

/* to one a letter */
char aa2letter(const char *resnm)
{
  struct { const char *str; char ch; } map[] = {
    {"GLY", 'G'}, {"ALA", 'A'}, {"VAL", 'V'}, {"LEU", 'L'}, {"ILE", 'I'},
    {"PHE", 'F'}, {"TYR", 'Y'}, {"TRP", 'W'}, {"PRO", 'P'},
    {"SER", 'S'}, {"THR", 'T'}, {"MET", 'M'}, {"CYS", 'C'}, {"CYX", 'C'}, {"CYM", 'C'},
    {"LYS", 'K'}, {"LYP", 'K'}, {"ARG", 'R'},
    {"HIS", 'H'}, {"HID", 'H'}, {"HIP", 'H'}, {"HIE", 'H'},
    {"GLU", 'E'}, {"GLN", 'Q'}, {"ASP", 'D'}, {"ASN", 'N'},
    {NULL, 'X'}};
  int i;

  if (strlen(resnm) == 4) resnm++;
  for (i = 0; map[i].str != NULL; i++) {
    if (strcmp(resnm, map[i].str) == 0)
      return map[i].ch;
  }
  fprintf(stderr, "Error: cannot abbreviate %s\n", resnm);
  return 'X';
}

/* recalc bb statistics */
static void bbrecalc(bbdih_t *bb, int bbcnt, hist_t *hs_fs)
{
  int i, j, i1, i2;
  double x, smphi, smphi2, smpsi, smpsi2, wf, ws;

  die_if (2*(bbcnt+1) != hs_fs->rows, "size mismatch %d != %d\n", 2*bbcnt+2, hs_fs->rows);
  for (i = 0; i < bbcnt; i++) {
    i1 = i*2 + 2;
    i2 = i*2 + 3;
    smphi = smphi2 = 0.;
    smpsi = smpsi2 = 0.;
    for (j = 0; j < hs_fs->n; j++) {
      x = hs_fs->xmin + (j+.5)*hs_fs->dx;
      wf = hs_fs->arr[i1*hs_fs->n + j];
      smphi += wf*x;
      smphi2 += wf*x*x;
      ws = hs_fs->arr[i2*hs_fs->n + j];
      smpsi += ws*x;
      smpsi2 += ws*x*x;
    }
    bb[i].smphi = smphi;
    bb[i].smphi2 = smphi2;
    bb[i].smpsi = smpsi;
    bb[i].smpsi2 = smpsi2;
  }
}

/* write dihp.txt */
static int dihwrite(bbdih_t *bb, int bbcnt,
    hist_t *hs_reg, const char *fn)
{
  FILE *fp;
  int i, j, reg, iregmax;
  double sum, phiav, phidev, psiav, psidev;
  double cfrac = 0., cfrac1 = 0., cnt = 0, *arr;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open file %s.\n", fn);
    return -1;
  }
  for (i = 0; i < bbcnt; i++) {
    fprintf(fp, "%3d %c ", bb[i].resid+1, aa2letter(bb[i].resnm));
    arr = hs_reg->arr + i*REGIONS;
    for (sum = 0., j = 0; j < REGIONS; j++)
      sum += arr[j];
    /* probability of being in each dihedral region */
    for (j = 0; j < REGIONS; j++)
      fprintf(fp, "%g ", arr[j]/sum);
    reg = bb[i].regr;
    /* the probability of being in the right dihedral region */
    fprintf(fp, "%g %g %d %g ",
        bb[i].phir, bb[i].psir, reg, arr[reg]/sum);
    for (iregmax = 0, j = 1; j < REGIONS; j++)
      if (arr[j] > arr[iregmax])
        iregmax = j;
    if (iregmax == reg) cfrac1 += 1.0;
    cfrac += arr[reg]/sum;
    cnt += 1;
    phiav  = bb[i].smphi/sum;
    phidev = sqrt(bb[i].smphi2/sum - phiav*phiav);
    psiav  = bb[i].smpsi/sum;
    psidev = sqrt(bb[i].smpsi2/sum - psiav*psiav);
    fprintf(fp, "%g %g %g %g\n", phiav, phidev, psiav, psidev);
  }
  printf("correct fraction %g, %g, cnt %g\n", cfrac/cnt, cfrac1/cnt, cnt);
  fprintf(fp, "# %g %g %g\n", cfrac/cnt, cfrac1/cnt, cnt);
  fclose(fp);
  return 0;
}

int run(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat, bbcnt;
  int satmax = 0, helmax = 0; /* maximal number of correct residues and maximal number helices */
  char fxtc[FILENAME_MAX];
  char *ftop;
  matrix box;
  rvec *xref = NULL;
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  bbdih_t *bb;
  double dt, bref = 1.0/(BOLTZ*Tref), wtot = 0.0;
  hist_t *hs_fs = NULL, *hs_reg = NULL;

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  bb = dihindex(&bbcnt, mtop); /* get indices for backbone dihedrals */
  for (i = 0; i < bbcnt; i++) /* compute bb[0..bbcnt-1].regr */
    bb[i].regr = calcdih(xref, bb[i].ia, &(bb[i].phir), &(bb[i].psir));
  hs_fs = hs_open(2+bbcnt*2, XMIN, XMAX, XDEL);
  hs_reg = hs_open(1+bbcnt, 0, REGIONS, 1.0);

  if (fndist2 != NULL) {
    die_if (0 != hs_load(hs_fs, fndist, HIST_VERBOSE),
        "cannot load master %s\n", fndist);
    die_if (0 != hs_load(hs_fs, fndist2, HIST_VERBOSE|HIST_ADDITION),
        "cannot add the second %s\n", fndist2);
    bbrecalc(bb, bbcnt, hs_fs);
  }
  if (fnreg2 != NULL) {
    die_if (0 != hs_load(hs_reg, fnreg, HIST_VERBOSE),
        "cannot load master %s\n", fnreg);
    die_if (0 != hs_load(hs_reg, fnreg2, HIST_VERBOSE|HIST_ADDITION),
        "cannot add the second %s\n", fnreg2);
  }

  if (fndist2 == NULL && fnreg2 == NULL) {
    if (!ctmd) { /* initialize weights for tempering */
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
    if (fnlog) {
      if ((fplog = fopen(fnlog, "w")) == NULL) {
        fprintf(stderr, "cannot open log %s\n", fnlog);
        return -1;
      }
    }
    if (fnbracket) {
      if (0 != loadbras(fnbracket))
        fnbracket = NULL;
    }
    for (i = 1; i <= ndata; i++) {
      sprintf(fxtc, "data%d/%s", i, fnxtc);
      if (!fexists(fxtc)) { /* test if the file exists */
        die_if(!ctmd, "cannot read %s\n", fnxtc);
        break;
      }
      if (0 != dotrj(bb, bbcnt,
            hs_fs, hs_reg,
            trj, fxtc, nat, box, &satmax, &helmax, &wtot)) {
        fprintf(stderr, "error doing %s\n", fnxtc);
        return -1;
      }
    }
  }
  dihwrite(bb, bbcnt, hs_reg, fnout);
  hs_save(hs_fs, fndist, HIST_ADDAHALF|HIST_VERBOSE);
  hs_close(hs_fs);
  hs_save(hs_reg, fnreg, HIST_VERBOSE);
  hs_close(hs_reg);
  sfree(mtop);
  if (!ctmd) {
    trj_free(trj);
    ze_free(ze);
  }
  if (fnlog) { fclose(fplog); fplog = NULL; }
  return 0;
}

int main(int argc, char **argv)
{
  (void) argtp;
  if (0 != doargs(argc, argv))
    return -1;
  run();
  return 0;
}


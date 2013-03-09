/* example of using multiple-histogram
 * we compute among alpha carbon atoms,
 * a histogram of radius of gyration, and
 * radial distribution function */
#include "zeutil.h" /* multiple-histogram module */

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnout = "out.txt";
const char *fnrdf = "rdf.txt";
double tequil = 1000; /* number of frames as equilibration */
double Tref = 300; /* reference temperature */
double delbet = 0.02;

/* print usage and then quit */
static void usage(const char *prog)
{
  fprintf(stderr, "%s [OPTIONS] file.xtc\n", prog);
  fprintf(stderr,
      " -s: followed by .tpr\n"
      " -f: followed by trajectory .xtc\n"
      " -z: followed by ZE file\n"
      " -t: followed by TRACE file\n"
      " -o: followed by output file\n"
      " -r: followed by rdf output\n"
      " -q: followed by time of equilibration\n"
      " -d: followed by beta cutoff\n"
      " -T: followed by reference temperature\n"
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
    if (strchr("fstzorqdT", ch) != NULL) {
      param = argv[i] + 2;
      if (param[0] == '\0') {
        if (i >= argc - 1)
          usage(argv[0]);
        param = argv[++i];
      }
      if (ch == 'f') fnxtc = param;
      else if (ch == 's') fntop = param;
      else if (ch == 't') fntr = param;
      else if (ch == 'z') fnze = param;
      else if (ch == 'o') fnout = param;
      else if (ch == 'r') fnrdf = param;
      else if (ch == 'q') tequil = atof(param);
      else if (ch == 'd') delbet = atof(param);
      else if (ch == 'T') Tref = atof(param);
      continue;
    }
    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'h': default: usage(argv[0]); break;
      }
  }
  return 0;
}

typedef struct { int idx, resid; } ca_t;
#define XMAX 5.0
#define XDEL 0.005
#define XCNT (int)(XMAX/XDEL + .5)
double grhist[XCNT+1];
double rdfhist[XCNT+1];

/* compute radius of gyration */
static double calcgr(ca_t *ca, int n, rvec x[])
{
  int i, id, j;
  double xc[3] = {0., 0., 0.}, x2, dx;

  /* compute the center */
  for (i = 0; i < n; i++) {
    id = ca[i].idx;
    for (j = 0; j < 3; j++)
      xc[j] += x[id][j];
  }
  for (j = 0; j < 3; j++) xc[j] /= n;
  /* compute deviations */
  for (x2 = 0., i = 0; i < n; i++) {
    id = ca[i].idx;
    for (j = 0; j < 3; j++) {
      dx = x[id][j] - xc[j];
      x2 += dx*dx;
    }
  }
  return sqrt(x2/n);
}

/* accumulate data for radial distribution function */
static void addrdf(ca_t *ca, int n, rvec x[], double w)
{
  int i, j, k, id, jd;
  double x2, dx;

  for (i = 0; i < n; i++) {
    id = ca[i].idx;
    for (j = i+1; j < n; j++) {
      jd = ca[j].idx;
      for (x2 = 0., k = 0; k < 3; k++) {
        dx = x[id][k] - x[jd][k];
        x2 += dx*dx;
      }
      dx = sqrt(x2);
      k = (int)(dx/XDEL);
      if (k < XCNT) rdfhist[k] += w;
    }
  }
}

/* add backbone index-set from moltype imt */
static ca_t *addca(ca_t *ca, int *cacnt, gmx_mtop_t *mtop, int imt, int idx,
    int resid)
{
  int ag, im, imb, nmb = mtop->nmolblock, nca;
  gmx_molblock_t *mb;

  nca = *cacnt;
  /* search molblock for those whose moltype is imt */
  for (ag = 0, imb = 0; imb < nmb; imb++) {
    mb = mtop->molblock + imb;
    if (mb->type == imt) {
      for (im = 0; im < mb->nmol; im++) {
        xrenew(ca, nca+1);
        ca[nca].idx = ag + im*mb->natoms_mol + idx;
        ca[nca].resid = resid;
        nca++;
      }
    }
    ag += mb->nmol * mb->natoms_mol;
  }
  *cacnt = nca;
  return ca;
}

/* obtain C-alpha indices */
static ca_t *ca_index(int *cacnt, gmx_mtop_t *mtop)
{
  gmx_moltype_t *mt;
  int imt, i, resid;
  ca_t *ca = NULL;

  xnew(ca, 1);
  *cacnt = 0;
  /* loop over moltypes to search a interaction that matches `allatoms' */
  for (imt = 0; imt < mtop->nmoltype; imt++) {
    mt = mtop->moltype + imt;
    for (i = 0; i < mt->atoms.nr; i++) {
      if (strcmp("CA", mt->atoms.atomname[i][0]) != 0)
        continue;
      resid = mt->atoms.atom[i].resnr;
      ca = addca(ca, cacnt, mtop, imt, i, resid);
    }
  } /* loop over moltypes */
  printf("C-alpha atoms:\n");
  for (i = 0; i < *cacnt; i++)
    printf("%d ", ca[i].idx+1);
  printf("\nTotal: %d atoms\n", *cacnt);
  return ca;
}

/* scan trajectory and accumulate data */
static int dotrj(ca_t *ca, int cacnt,
    trj_t *trj, const char *fntrj, int nat, matrix box)
{
  int nfr, i, match, natoms;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL;
  double t1, w, gr;
  real t = 0.f;

  /* loop over xtc frames */
  for (nfr = 0; ; nfr++) {
    if (nfr == 0) {
      natoms = read_first_x(&fpxtc, fntrj, &t, &x, box);
      if (natoms != nat) {
        fprintf(stderr, "natoms mismatch! top %d, xtc %d\n", nat, natoms);
        return -1;
      }
    } else {
      if (!read_next_x(fpxtc, &t, natoms, x, box))
        break;
    }
    if (t < tequil) continue;

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
    if (!trj->fr[trj->pos].in) continue;
    w = exp(trj->fr[trj->pos].lnw);

    gr = calcgr(ca, cacnt, x); /* compute radius of gyration */
    i = (int)(gr/XDEL);
    if (i < XCNT) grhist[i] += w;

    addrdf(ca, cacnt, x, w);
  }
  close_trj(fpxtc);
  return 0;
}

static int writehist(const char *fn, double arr[])
{
  FILE *fp;
  int i, i0, i1;
  double fac;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open file %s.\n", fn);
    return -1;
  }
  for (i0 = 0; i0 < XCNT; i0++) if (arr[i0] > 0.) break;
  for (i1 = XCNT-1; i1 >= i0; i1--) if (arr[i1] > 0.) break;
  for (fac = 0., i = i0; i <= i1; i++) fac += arr[i];
  fac = 1.0/(fac*XDEL);
  for (i = i0; i <= i1; i++) {
    fprintf(fp, "%g %g %g\n", XDEL*(i+.5), arr[i], arr[i]*fac);
  }
  fclose(fp);
  return 0;
}

int run(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat, cacnt;
  double dt;
  char fxtc[FILENAME_MAX];
  matrix box;
  rvec *xref = NULL;
  zedata_t *ze;
  trj_t *trj;
  ca_t *ca;
  double bref = 1.0/(BOLTZ*Tref);

  /* read topology */
  mtop = read_tpx_conf(fntop, &xref, box, &nat, &dt);
  ca = ca_index(&cacnt, mtop); /* get indices for backbone dihedrals */

  if ((ze = ze_load(fnze, 10.0)) == NULL) {
    fprintf(stderr, "cannot load %s\n", fnze);
    return -1;
  }
  /* we first load TRACE, compute weight of each frame */
  if ((trj = trj_loadseq(fntr, &ndata, dt, tequil,
          ze, bref, delbet)) == NULL) {
    fprintf(stderr, "failed to load trace\n");
    return -1;
  }
  /* we now successively load xtc, try to match the TRACE frames */
  for (trj->pos = 0, i = 1; i <= ndata; i++) {
    sprintf(fxtc, "data%d/%s", i, fnxtc);
    if (0 != dotrj(ca, cacnt, trj, fxtc, nat, box)) {
      fprintf(stderr, "error doing %s\n", fnxtc);
      return -1;
    }
  }
  writehist(fnout, grhist);
  writehist(fnrdf, rdfhist);
  sfree(mtop);
  trj_free(trj);
  ze_free(ze);
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


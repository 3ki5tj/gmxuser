/* compute contact number
 * Copyright (C) 2010, 2011 Cheng Zhang */
#include "zeutil.h"
#include "ztoputil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RV3
#define ZCOM_HIST
#include "zcom.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fnze = "ZE";
const char *fntr = "TRACE";
const char *fnidx = "core.idx";
const char *fncont = "cont.out";
const char *fncont2 = NULL;
const char *fnlog = "cont.log";
const char *fnrdf = "core.rdf";
double tequil = 10000;
double rcutoff = 1.2; /* cutoff distance to be deemed as a contact, in terms of nm */
double Tref = 300; /* reference temperature in Kelvin */
double delbet = 1.0; /* we use the entire range, safer: 0.02 */
int ctmd = 0; /* constant temperature md */
int logfreq = 1;
int CONTMAX = 2000;
double RDFMAX = 20.0;
double RDFDEL = 0.02;

/* print usage and then quit */
static void usage(const char *prog)
{
  fprintf(stderr, "%s [OPTIONS] file.xtc\n", prog);
  fprintf(stderr, "To combine two database:\n  %s -o cont1 -a cont2\n", prog);
  fprintf(stderr,
      " -f: followed by trajectory .xtc\n"
      " -s: followed by .tpr\n"
      " -i: followed by .idx\n"
      " -z: followed by ZE file\n"
      " -t: followed by TRACE file\n"
      " -o: followed by cont.out file\n"
      " -q: followed by time of equilibration\n"
      " -c: followed by contact cutoff in nm\n"
      " -r: followed by rdf file\n"
      " -a: followed by the second database\n"
      " -T: followed by the reference temperature\n"
      " -M: constant temperature simulation\n"
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
    if (strchr("ifsztoaqdTcgr", ch) != NULL) {
      param = argv[i] + 2;
      if (param[0] == '\0') {
        if (i >= argc - 1)
          usage(argv[0]);
        param = argv[++i];
      }
      if (ch == 'f') fnxtc = param;
      else if (ch == 's') fntop = param;
      else if (ch == 'i') fnidx = param;
      else if (ch == 'z') fnze = param;
      else if (ch == 't') fntr = param;
      else if (ch == 'o') fncont = param;
      else if (ch == 'a') fncont2 = param;
      else if (ch == 'r') fnrdf = param;
      else if (ch == 'q') tequil = atof(param);
      else if (ch == 'T') Tref = atof(param);
      else if (ch == 'd') delbet = atof(param);
      else if (ch == 'c') rcutoff = atof(param);
      else if (ch == 'g') fnlog = param;
      else if (ch == 'F') {
        logfreq = atoi(param);
        if (logfreq < 1) logfreq = 1;
        fprintf(stderr, "frequency of writing log %d\n", logfreq);
      }
      continue;
    }
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      switch (ch) {
      case 'M': ctmd = 1; break;
      case 'h': default: usage(argv[0]); break;
      }
    }
  }
  if (fnlog) delbet = 1e3;

  /* check if it is constant temperature md */
  if (!ctmd && !fexists("data1/ZE")) {
    fprintf(stderr, "no ZE, assume constant temperature md\n");
    ctmd = 1;
  }
  return 0;
}

typedef struct {
  int n;
  double vol, wtot;
} den_t;

#define sqr(x) ((x)*(x))
static __inline double dist(rvec a, rvec b)
{
  return sqrt(sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]));
}

/* count contact number for each group pair, save to
 * ct[0..ngg-1], where ngg = ng*(ng-1)/2 */
static int countcont(ca_t *ca, int cacnt, int ng,
    rvec *x, int *ct, hist_t *hs_rdf, double w)
{
  int i, j, ia, ja, ig, jg, igg, ngg, ctall;
  double dis;

  ngg = ng*(ng-1)/2;
  for (igg = 0; igg < ngg; igg++) ct[igg] = 0;
  for (i = 0; i < cacnt; i++) {
    for (j = i+1; j < cacnt; j++) {
      ia = ca[i].idx;
      ja = ca[j].idx;
      ig = ca[i].gid;
      jg = ca[j].gid;
      if (ig < 0 || jg < 0 || ig == jg) continue;
      dis = rv3_dist(x[ia], x[ja]); /* assume protein is whole */
      hs_add(hs_rdf, &dis, w, HIST_VERBOSE);
      /* compute distance between ia and ja */
      if (dis > rcutoff) continue;
      igg = getpairindex(ig, jg, ng);
      ct[igg] += 1;
    }
  }
  for (ctall = 0, igg = 0; igg < ngg; igg++) ctall += ct[igg];
  return ctall;
}

/* do trajectory and related trace */
static int dotrj(protop_t *pro, hist_t *hs_ct, hist_t *hs_rdf, den_t *den,
    trj_t *trj, const char *fxtc, int nat, matrix box,
    FILE *fplog)
{
  int nfr, match, natoms;
  int fpxtc = 0; /* gromacs file handle for trj */
  double t1, w;
  real t = 0.f;
  rvec *x = NULL;
  int igg, *ct, ctall;
  double *ctdbl;

  xnew(ct, pro->ngg+1);
  xnew(ctdbl, pro->ngg+1);

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
    if (t < tequil) continue;

    /* search trj for a matching frame */
    if (ctmd) {
      w = 1.;
    } else {
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
      if (!trj->fr[trj->pos].in) continue; /* skip irrelavant frame */
      w = exp(trj->fr[trj->pos].lnw);
    }

    ctall = countcont(pro->ca, pro->nca, pro->ngrp, x, ct+1, hs_rdf, w);
    ct[0] = ctall; /* ct 0 reserved for total contact */
    for (igg = 0; igg <= pro->ngg; igg++) ctdbl[igg] = ct[igg];
    hs_add(hs_ct, ctdbl, w, HIST_VERBOSE);
    den->wtot += w;

    if (fplog && (((int)(t+.5)) % logfreq == 0)) {
      fprintf(fplog, "%.0f ", t);
      for (igg = 0; igg <= pro->ngg; igg++)
        fprintf(fplog, "%d ", ct[igg]);
      fprintf(fplog, "%g\n", w);
    }
  }
  close_xtc(fpxtc);
  free(ct);
  free(ctdbl);
  return 0;
}

static int rdfw(FILE *fp, void *pdata)
{
  den_t *den = (den_t *)pdata;
  fprintf(fp, "RDF %d %g %g | ", den->n, den->wtot, den->vol);
  return 0;
}

static int rdfr(const char *s, void *pdata)
{
  den_t *den = (den_t *)pdata;
  int ret = sscanf(s, " RDF %d%lf%lf | ", &(den->n), &(den->wtot), &(den->vol));
  return (ret == 3) ? 0 : 1;
}

static double rdfnorm(int j, int i, double xmin, double dx, void *pdata)
{
  den_t *den = (den_t *)pdata;
  int npr;
  double x, fac, vsph;

  (void)j;
  x = xmin + i*dx;
  vsph = (4.*M_PI/3)*dx*(3*x*(x+dx) + dx*dx);
  npr = den->n*(den->n-1)/2;
  fac = den->vol/vsph/npr/den->wtot;
  return fac;
}


int run(void)
{
  gmx_mtop_t *mtop;
  char *ftop;
  int i, nat, ndata;
  char fxtc[FILENAME_MAX];
  double dt;
  double bref = 1.0/(BOLTZ*300);
  rvec *xref = NULL; /* coordinates of .tpr */
  matrix box;
  zedata_t *ze;
  trj_t *trj;

  protop_t *pro;
  hist_t *hs_ct, *hs_rdf;
  den_t den[1];
  FILE *fplog;

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  if (fnidx) fnidx = nfexists(fnidx, 7);

  if ((pro = pro_init(fnidx, mtop, xref,
          PROTOP_SEPHELIX|PROTOP_HYDROPHOBIC)) == NULL) {
    fprintf(stderr, "cannot initialize groups\n");
    return -1;
  }
  printf("%d groups, %d pairs, %d group atoms\n", pro->ngrp, pro->ngg, pro->ngat);
  hs_ct = hs_open(1+pro->ngg, 0, CONTMAX, 1);
  den->n = pro->ngat;
  den->wtot = 0.;
  den->vol = box[0][0]*box[1][1]*box[2][2];
  hs_rdf = hs_openx(1, 0., RDFMAX, RDFDEL, rdfw, rdfr, rdfnorm);

  /* do database addition */
  if (fncont2 != NULL) {
    die_if (0 != hs_load(hs_ct, fncont, HIST_VERBOSE),
      "cannot load master: %s\n", fncont);
    die_if (0 != hs_load(hs_ct, fncont2, HIST_VERBOSE|HIST_ADDITION),
       "cannot add additional: %s\n", fncont2);
  } else {
    if (!ctmd) {
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
    } else {
      ndata = 10000;
    }

    if (fnlog && (fplog = fopen(fnlog, "w")) == NULL) {
      fprintf(stderr, "cannot open log file %s\n", fnlog);
      return -1;
    }
    /* we now successively load xtc, try to match the TRACE frames */
    if (!ctmd) trj->pos = 0;
    for (i = 1; i <= ndata; i++) {
      sprintf(fxtc, "data%d/%s", i, fnxtc);
      if (!fexists(fxtc)) {
        die_if(!ctmd, "cannot open %s\n", fnxtc);
        break;
      }
      if (0 != dotrj(pro, hs_ct, hs_rdf, den,
            trj, fxtc, nat, box, fplog)) {
        fprintf(stderr, "error doing %s\n", fnxtc);
        return -1;
      }
    }
    if (!ctmd) {
      trj_free(trj);
      ze_free(ze);
    }
    if (fplog) fclose(fplog);
  }
  hs_save(hs_ct, fncont, HIST_KEEPHIST|HIST_VERBOSE);
  hs_close(hs_ct);
  hs_savex(hs_rdf, fnrdf, den, HIST_ADDAHALF|HIST_KEEPHIST);
  hs_close(hs_rdf);
  sfree(mtop);
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


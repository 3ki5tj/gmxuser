/* make multiple histogram index
 * Copyright (C) 2011 Cheng Zhang */
#include "zeutil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#include "zcom1.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnout = "mh.idx";
double tequil = 0;/* handle arguments */

double Tref = 300; /* reference temperature */
double delbet = 1000.0; /* we don't skip frames by default */
int every = 1;
double weightmin = 0.05; /* minimal WHAM weight to be indexed */

static int doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "generate multiple histogram index at 300K";
  argopt_add(ao, "-s", NULL,  &fntop,  "topology file");
  argopt_add(ao, "-f", NULL,  &fnxtc,  "trajectory");
  argopt_add(ao, "-z", NULL,  &fnze,   "ZE file");
  argopt_add(ao, "-t", NULL,  &fntr,   "TRACE file");
  argopt_add(ao, "-o", NULL,  &fnout,  "output file");
  argopt_add(ao, "-T", "%lf", &Tref,   "reference temperature");
  argopt_add(ao, "-e", "%d",  &every,  "skip interval");
  argopt_add(ao, "-w", "%lf", &weightmin, "minimal weight");
  argopt_add(ao, "-q", "%lf", &tequil, "time of equilibration");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);

  /* check if it is constant temperature md */
  die_if (!fexists("data1/ZE"), "cannot found ZE");
  return 0;
}

static int writeidx(trj_t *trj)
{
  int i, id = 0;
  double lnw, lnwmin = log(weightmin);
  FILE *fp;

  xfopen(fp, fnout, "w", exit(1));
  for (id = 0, i = 0; i < trj->n; i++) {
    if (!trj->fr[i].in) continue;
    lnw = trj->fr[i].lnw;
    if (lnw < lnwmin) continue;
    fprintf(fp, "%8d %9.0f %g %d\n", id++, trj->fr[i].t, exp(lnw), trj->fr[i].dataid);
  }
  fclose(fp);
  return 0;
}

static int domoltrj(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat;
  char fxtc[FILENAME_MAX], *ftop;
  matrix box;
  rvec *xref = NULL;
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  double dt, bref = 1.0/(BOLTZ*Tref);

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);

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

  writeidx(trj);
  sfree(mtop);
  if (trj) trj_free(trj);
  if (ze) ze_free(ze);
  return 0;
}

int main(int argc, char **argv)
{
  (void) argtp;
  doargs(argc, argv);
  domoltrj();
  return 0;
}


/* vdist.c
 * 1. compute vertical distance, from one helix to another two,
 *    a quantity that represents helix packing chirality
 * 2. compute helix distance
 * 3. compute helix angle*/
#include "zeutil.h" /* multiple-histogram module */
#include "ztoputil.h"

#define ZCOM_PICK
#define ZCOM_SS
#define ZCOM_HIST
#include "zcom.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnvdis = "v.dist";
const char *fnpdis = "p.dist";
const char *fngdis = "p.ang";
const char *fnidx = "helix.idx";
const char *fnvadd = NULL;
const char *fnpadd = NULL;
const char *fngadd = NULL;
char *fnvdis2 = NULL;  /* correlation, append fnvdis by fnvdis by "2" */
char *fnvadd2 = NULL;
double tequil = 1000; /* number of frames as equilibration */
double Tref = 300; /* reference temperature */
double delbet =  1000.0; /* 0.10; */
const char *fnlog = "vdist.log";  /* log file name */
int logfreq = 1;
int ctmd = 0;

/* in angstrom */
double VDMAX = 35.0;
double VDDEL = 0.1;
double PDMAX = 80.0;
double PDDEL = 0.05;

/* print usage and then quit */
static void usage(const char *prog)
{
  fprintf(stderr, "Usage:\n  %s [OPTIONS] file.xtc\n", prog);
  fprintf(stderr, "To combine data:\n  %s master.dist -a add.dist\n",
      prog);
  fprintf(stderr,
      " -s: followed by .tpr\n"
      " -f: followed by trajectory .xtc\n"
      " -z: followed by ZE file\n"
      " -t: followed by TRACE file\n"
      " -o: followed by vdist (chiral) distribution output file\n"
      " -p: followed by pdist (helix-helix distance) distribution output file\n"
      " -n: followed by pang (helix-helix angle) distribution output file\n"
      " -q: followed by time of equilibration\n"
      " -d: followed by beta cutoff\n"
      " -T: followed by reference temperature\n"
      " -a: adding the second database for vdist\n"
      " -b: adding the second database for pdist\n"
      " -c: adding the second database for pang\n"
      " -g: followed by the log file\n"
      " -i: followed by index file\n"
      " -F: followed by the frequence of writing log file\n"
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
    if (strchr("fsitzopnrqdTabcgF", ch) != NULL) {
      param = argv[i] + 2;
      if (param[0] == '\0') {
        if (i >= argc - 1)
          usage(argv[0]);
        param = argv[++i];
      }
      if (ch == 'f') fnxtc = param;
      else if (ch == 's') fntop = param;
      else if (ch == 'i') fnidx = param;
      else if (ch == 't') fntr = param;
      else if (ch == 'z') fnze = param;
      else if (ch == 'o') fnvdis = param;
      else if (ch == 'p') fnpdis = param;
      else if (ch == 'n') fngdis = param;
      else if (ch == 'q') tequil = atof(param);
      else if (ch == 'd') delbet = atof(param);
      else if (ch == 'T') Tref = atof(param);
      else if (ch == 'a') fnvadd = param;
      else if (ch == 'b') fnpadd = param;
      else if (ch == 'c') fngadd = param;
      else if (ch == 'g') fnlog = param;
      else if (ch == 'F') {
        logfreq = atoi(param);
        if (logfreq < 1) logfreq = 1;
        fprintf(stderr, "frequency of writing log %d\n", logfreq);
      }
      continue;
    }
    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'M': ctmd = 1; break;
      case 'h': default: usage(argv[0]); break;
      }
  }
  if (fnlog) delbet = 1e3; /* no saving of delbet */

  /* name transformations */
  fnvdis2 = ssdup(fnvdis); sscat(fnvdis2, "2");
  fnvadd2 = ssdup(fnvadd); sscat(fnvadd2, "2");
  /* check if it is constant temperature md */
  if (!ctmd && !fexists("data1/ZE")) {
    fprintf(stderr, "no ZE, assume constant temperature md\n");
    ctmd = 1;
  }
  return 0;
}

/* model helix as a rod from ps to pt
 * coordinates collected from C-alpha atom is..it
 * save modeled endpoints to ps and pt */
static void getrod(const protop_t *pro, int ig, rvec x[],
    rvec ps, rvec pt)
{
  int i, id, n, j, k, is = 9999, it = -1;
  real xc[3], dx[3], mat[3][3], val[3], vecs[3][3], fac;

  /* get center of mass */
  rv3_zero(xc);
  for (n = 0, i = 0; i < pro->nca; i++) {
    if (pro->ca[i].gid != ig) continue;
    id = pro->ca[i].idx;
    rv3_inc(xc, x[id]);
    if (i < is) is = i;
    if (i > it) it = i;
    n++;
  }
  die_if (n < 3, "helix %d is too short, n = %d\n", ig, n);
  rv3_smul(xc, 1.f/n);
  /* get moments of inertia */
  for (j = 0; j < 3; j++)
    for (k = 0; k < 3; k++)
      mat[j][k] = 0.f;
  for (i = 0; i < pro->nca; i++) {
    if (pro->ca[i].gid != ig) continue;
    id = pro->ca[i].idx;
    rv3_diff(dx, x[id], xc);
    for (j = 0; j < 3; j++)
      for (k = j; k < 3; k++)
        mat[j][k] += dx[j]*dx[k];
  }
  for (j = 0; j < 3; j++) {
    for (k = j+1; k < 3; k++)
      mat[k][j] = mat[j][k] = mat[j][k]/n;
    mat[j][j] *= 1.f/n;
  }
  /* eigenvector of the largest eigenvalue */
  rm3_eigsys(val, vecs, mat, 1);
  /* determine the sign */
  fac = (real) sqrt(3.f*val[0]);
  rv3_diff(dx, x[ pro->ca[it].idx ], x[ pro->ca[is].idx ]);
  if (rv3_dot(dx, vecs[0]) < 0.) fac = -fac;
  /* construct end points */
  rv3_lincomb2(ps, xc, vecs[0], 1.f, -fac);
  rv3_lincomb2(pt, xc, vecs[0], 1.f, fac);
/*
  printf("val = {%g,%g,%g}, fac = %g\n", val[0],val[1],val[2], fac);
  printf("xc %g,%g,%g, ps %g,%g,%g pt %g,%g,%g\n",
      xc[0],xc[1],xc[2], ps[0],ps[1],ps[2], pt[0],pt[1],pt[2]);
  printf("%g,%g,%g --> %g,%g,%g\n",
      x[pro->ca[is].idx][0], x[pro->ca[is].idx][1], x[pro->ca[is].idx][2],
      x[pro->ca[it].idx][0], x[pro->ca[it].idx][1], x[pro->ca[it].idx][2]);
  getchar();
*/
}

/* model helices as rods
 * save end points of helix i in rods[2*i] and rods[2*i+1] */
static void getrods(const protop_t *pro,
    rvec x[], rvec *rods)
{
  int ip;

  for (ip = 0; ip < pro->ngrp; ip++) {
    getrod(pro, ip, x, rods[2*ip], rods[2*ip+1]);
  }
}

/* compute vertical distance for three successive helices L M R
 * for each triplet, we have a distance from L to the plane of M and R
 * and that from R to plane of L and M.
 * thus we have 2*(nse - 2) quantities to output, nse is the number of helices
 * */
static void calcvd(double ap[], double am[],
    const ca_t *ca, int n,
    int nse, rvec *rods, rvec x[])
{
  int ip, cnt, i1, i4, i1s, i4s;
  double r1, r2;
  real *ps, *pt;

  /* loop over pivots */
  for (ip = 1; ip < nse - 1; ip++) {
    ps = rods[ip*2];
    pt = rods[ip*2+1];
    cnt = 0;
    ap[ip-1] = 0.;
    am[ip-1] = 0.;
    for (i1 = 0; i1 < n; i1++) {
      if (ca[i1].gid != ip-1) continue;
      i1s = ca[i1].idx;
      for (i4 = 0; i4 < n; i4++) {
        if (ca[i4].gid != ip+1) continue;
        i4s = ca[i4].idx;
        r1 = rv3_vpdist(x[i4s],   x[i1s], ps, pt);
        r2 = rv3_vpdist(x[i1s],   x[i4s], pt, ps);
        ap[ip-1] += r1*10; /* x 10 from nm to angstrom */
        am[ip-1] += r2*10;
        cnt++;
      }
    }
    ap[ip-1] /= cnt;
    am[ip-1] /= cnt;
  }
}

/* compute pair distance
 * hd[0..nse*(nse-1)/2-1] pair distance between helices
 *  ca[0..n-1]: index map from residue index to atom index
 *  se[0..2*nse]: index map from helix index to residue index */
static void calcpd(double hd[], double hg[],
    const ca_t *ca, int n,
    int nse, rvec *rods, rvec x[])
{
  int h1, h2, ipr, i, is, j, js, cnt1, cnt2;
  double dis, dmin, dsm, ang;
  rvec v1, v2;

  /* loop over pivots */
  for (ipr = 0, h1 = 0; h1 < nse - 1; h1++) { /* first helix */
    rv3_normalize(rv3_diff(v1, rods[h1*2], rods[h1*2+1]));

    for (h2 = h1+1; h2 < nse; h2++, ipr++) { /* second helix */
      /* A. compute distance between helices h1 and h2 */
      dsm = 0.;

      /* A1. for each point on h1 */
      for (cnt1 = 0, i = 0; i < n; i++) {
        if (ca[i].gid != h1) continue;
        is = ca[i].idx;
        /* minimal distance to any point h2 */
        for (dmin = 1e9, j = 0; j < n; j++) {
          if (ca[j].gid != h2) continue;
          js = ca[j].idx;
          dis = rv3_dist(x[is], x[js]);
          if (dis < dmin) dmin = dis;
        }
        dsm += dmin;
        cnt1++;
      }

      /* A2. for each point on h2 */
      for (cnt2 = 0, i = 0; i < n; i++) {
        if (ca[i].gid != h2) continue;
        is = ca[i].idx;
        for (dmin = 1e9, j = 0; j < n; j++) {
          if (ca[j].gid != h1) continue;
          js = ca[j].idx;
          dis = rv3_dist(x[is], x[js]);
          if (dis < dmin) dmin = dis;
        }
        dsm += dmin;
        cnt2++;
      }

      /* averaged distance between h1 and h2 */
      dsm /= (cnt1+cnt2);
      hd[ipr] = dsm*10.0; /* to angstrom */

      /* B. compute angle between h1 and h2 */
      rv3_normalize(rv3_diff(v2, rods[h2*2], rods[h2*2+1]));
      ang = acos(rv3_dot(v1, v2));
      hg[ipr] = ang*180.0/M_PI;
    }
  }
}

/* scan trajectory and accumulate data */
static int dotrj(protop_t *pro,
    hist_t *hs_vd, hist2_t *hs2_vd,
    hist_t *hs_pd, hist_t *hs_pg,
    trj_t *trj, const char *fntrj, int nat, matrix box,
    FILE *fplog)
{
  int nfr, i, match, natoms, cnt;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL;
  rvec *rods;
  double t1, w, *vp, *vm, *pd, *pg;
  real t = 0.f;

  cnt = pro->ngrp - 2;
  xnew(vp, cnt);
  xnew(vm, cnt);
  pro->ngg = pro->ngg;
  xnew(pd, pro->ngg);
  xnew(pg, pro->ngg);
  xnew(rods, pro->ngrp*2);

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

    /* compute helix rod models */
    getrods(pro, x, rods);

    /* A. vertical distance */
    calcvd(vp, vm, pro->ca, pro->nca, pro->ngrp, rods, x);
    /* 1D vdist histogram */
    for (i = 0; i < cnt; i++) {
      hs_add1(hs_vd, 2*i,   vp[i], w, HIST_VERBOSE);
      hs_add1(hs_vd, 2*i+1, vm[i], w, HIST_VERBOSE);
    }
    /* 2D vdist histogram, only vp to save space */
    if (hs2_vd != NULL)
      hs2_add(hs2_vd, vp, vp+1, 1, w, HIST_VERBOSE);

    /* B. pair distance */
    calcpd(pd, pg, pro->ca, pro->nca, pro->ngrp, rods, x);
    /* 1D pdist histogram */
    hs_add(hs_pd, pd, w, HIST_VERBOSE);
    hs_add(hs_pg, pg, w, HIST_VERBOSE);

    /* log */
    if (fplog && (((int)(t+.5)) % logfreq == 0)) {
      fprintf(fplog, "%.0f ", t);
      for (i = 0; i < pro->ngrp - 2; i++) /* vdist */
        fprintf(fplog, "%7.3f %7.3f ", vp[i], vm[i]);
      for (i = 0; i < pro->ngg; i++)  /* pdist */
        fprintf(fplog, "%7.3f ", pd[i]);
      for (i = 0; i < pro->ngg; i++)  /* pang */
        fprintf(fplog, "%6.2f ", pg[i]);
      fprintf(fplog, "%g\n", w);
    }
  }
  close_trj(fpxtc);
  free(pg);
  free(pd);
  free(vp);
  free(vm);
  free(rods);
  return 0;
}

int run(void)
{
  gmx_mtop_t *mtop;
  char *ftop, *fidx = NULL;
  int i, ndata, nat;
  double dt;
  double *vpref, *vmref;
  char fxtc[FILENAME_MAX];
  matrix box;
  rvec *xref = NULL;
  double bref = 1.0/(BOLTZ*Tref);
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  protop_t *pro;
  rvec *rods;
  int nvd;
  hist_t *hs_vd = NULL, *hs_pd = NULL, *hs_pg = NULL;
  hist2_t *hs2_vd = NULL;
  FILE *fplog = NULL;

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  if (fnidx) fidx = nfexists(fnidx, 7);

  pro = pro_init(fidx, mtop, xref, PROTOP_SEPHELIX);

  nvd = pro->ngrp-2;
  xnew(vpref, nvd);
  xnew(vmref, nvd);
  xnew(rods, pro->ngrp*2);
  getrods(pro, xref, rods);
  calcvd(vpref, vmref, pro->ca, pro->nca, pro->ngrp, rods, xref);
  for (i = 0; i < pro->ngrp -2; i++) {
    printf("ref. vd[%d-%d-%d] = %g, %g\n",
        i+1, i+2, i+3, vpref[i], vmref[i]);
  }
  /* histograms */
  hs_vd = hs_open(2*nvd, -VDMAX, VDMAX, VDDEL); /* vdist */
  if (pro->ngrp >= 4)
    hs2_vd = hs2_open(pro->ngrp-3, -VDMAX, VDMAX, VDDEL, -VDMAX, VDMAX, VDDEL);
  hs_pd = hs_open(pro->ngg, 0., PDMAX, PDDEL); /* pair dist */
  hs_pg = hs_open(pro->ngg, 0., 180.0, 1.0); /* pair angle */

  if (fnvadd != NULL) { /* database manipulation */
    die_if (0 != hs_load(hs_vd, fnvdis, HIST_VERBOSE),
      "cannot load master: %s\n", fnvdis);
    die_if (0 != hs_load(hs_vd, fnvadd, HIST_VERBOSE|HIST_ADDITION),
       "cannot add additional: %s\n", fnvadd);
    if (hs2_vd != NULL) {
      die_if (0 != hs2_load(hs2_vd, fnvdis2, HIST_VERBOSE),
        "cannot load master: %s\n", fnvdis2);
      die_if (0 != hs2_load(hs2_vd, fnvadd2, HIST_VERBOSE|HIST_ADDITION),
        "cannot add additional: %s\n", fnvadd2);
    }
  }

  if (fnpadd != NULL) { /* database manipulation */
    die_if (0 != hs_load(hs_pd, fnpdis, HIST_VERBOSE),
      "cannot load master: %s\n", fnpdis);
    die_if (0 != hs_load(hs_pd, fnpadd, HIST_VERBOSE|HIST_ADDITION),
       "cannot add additional: %s\n", fnpadd);
  }

  if (fngadd != NULL) { /* database manipulation */
    die_if (0 != hs_load(hs_pg, fngdis, HIST_VERBOSE),
      "cannot load master: %s\n", fnpdis);
    die_if (0 != hs_load(hs_pg, fngadd, HIST_VERBOSE|HIST_ADDITION),
       "cannot add additional: %s\n", fngadd);
  }
  if (fnpadd == NULL && fnvadd == NULL && fngadd == NULL) {
    /* do trajectories */

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
      if (0 != dotrj(pro,
            hs_vd, hs2_vd, hs_pd, hs_pg,
            trj, fxtc, nat, box, fplog)) {
        fprintf(stderr, "error doing %s\n", fnxtc);
        return -1;
      }
    }
    if (!ctmd) {
      ze_free(ze);
      trj_free(trj);
    }
    if (fplog) fclose(fplog);
  }

  /* write vdist data */
  hs_save(hs_vd, fnvdis, HIST_ADDAHALF|HIST_KEEPHIST|HIST_VERBOSE);
  hs_close(hs_vd);
  if (hs2_vd != NULL) {
    hs2_save(hs2_vd, fnvdis2, HIST_ADDAHALF|HIST_KEEPHIST|HIST_VERBOSE);
    hs2_close(hs2_vd);
  }

  /* write pdist data */
  hs_save(hs_pd, fnpdis, HIST_ADDAHALF|HIST_KEEPHIST|HIST_VERBOSE);
  hs_close(hs_pd);

  /* write angle data */
  hs_save(hs_pg, fngdis, HIST_ADDAHALF|HIST_KEEPHIST|HIST_VERBOSE);
  hs_close(hs_pg);

  sfree(mtop);
  free(rods);
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


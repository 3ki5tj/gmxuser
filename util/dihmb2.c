/* display dihedral information for an .xtc file
 * Copyright (C) 2011 Cheng Zhang */
#include "zeutil.h"
#include "ztoputil.h"

/* real is defined in GROMACS */
#define ZCOM_PICK
#define ZCOM_RV3
#define ZCOM_HIST
#define ZCOM_UTIL
#include "zcom.h"

#include "rotdb.h"

const char *fnxtc = "fit.xtc";
const char *fntop = "naked.tpr";
const char *fntr = "TRACE";
const char *fnze = "ZE";
const char *fnidx = "core.idx";
const char *fndist = "chi1.dist";
const char *fndist2 = NULL;
const char *fnout = "chi1.out";
const char *fndb  = "chi1.db";
const char *fndb2 = NULL;
const char *fnbracket = "fold.bra";
double tequil = 10000;
double Tref = 300; /* reference temperature */
double delbet = 0.2;
int ctmd = 0;
const char *fnlog = NULL;
FILE *fplog = NULL;

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
      " -i: followed by index file\n"
      " -D: followed by output distribution\n"
      " -a: followed by second distribution\n"
      " -o: followed by chi1 rotamer\n"
      " -g: followed by output log file\n"
      " -q: followed by time of equilibration\n"
      " -T: followed by reference temperature\n"
      " -d: followed by delta beta\n"
      " -b: followed by the bracket file\n"
      " -x: followed by rotamer database\n"
      " -y: followed by second rotamer database\n"
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
    if (strchr("fstziDoqdThgbaxy", ch) != NULL) {
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
      else if (ch == 'D') fndist = param;
      else if (ch == 'a') fndist2 = param; /* second database */
      else if (ch == 'o') fnout = param;
      else if (ch == 'i') fnidx = param;
      else if (ch == 'q') tequil = atof(param);
      else if (ch == 'd') delbet = atof(param);
      else if (ch == 'T') Tref = atof(param);
      else if (ch == 'g') fnlog = param;
      else if (ch == 'b') fnbracket = param;
      else if (ch == 'x') fndb = param;
      else if (ch == 'y') fndb2 = param; /* second database */
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

/* note we use degree instead of rad */
#define XMAX 180
#define XMIN (-XMAX)
#define XDEL 1.0

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

  if (!fexists(fn)) return -1;
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
    if (fgets(buf, sizeof buf, fp) == NULL) {
      fprintf(stderr, "%s: error of %dth brackets\n", fn, i+1);
      fclose(fp);
      return -1;
    }
    strip(buf);
    die_if (2 != sscanf(buf, "%lf%lf", &(bras[i][0]), &(bras[i][1])),
      "cannot scan bracket from line [%s]\n", buf);
  }
  fclose(fp);
  printf("brackets are loaded successfully\n");
  return 0;
}

/* angle to region */
static unsigned ang2int(double ang)
{
  if (fabs(ang) > 120.0) return 0;
  else if (ang < 0) return 1;
  else return 2;
}

/* code rotemar configurations into an interger */
static void encode_rotamers(unsigned *code, protop_t *pro, double *chi1)
{
  int ica, cnt, ib;
  unsigned base = 1u, icode = 0u, reg;
  /* encode to integer */
  die_if(pro->ngat > RBMAX*RBSIZ, "too many %d atoms to encode!\n", pro->ngat);
  for (ib = -1, cnt = 0, ica = 0; ica < pro->nca; ica++) {
    if (pro->ca[ica].gid < 0) continue;
    if (cnt % RBSIZ == 0) {
      ib++;
      icode = 0u;
      base = 1u;
    }
    //die_if (cnt > 20, "too many %d atoms to encode\n", cnt);
    reg = ang2int(chi1[ica]);
    icode += base*reg;
    code[ib] = icode;
    base *= 3;
    cnt++;
    //printf("%3d %u %u\n", cnt, reg, icode);
  }
  //printf("%3d %u %u %u\n", ib, icode, code[0], code[1]); getchar();
}

/* compute chi1 */
static void getchi1(unsigned *code, protop_t *pro, rvec *x, double *chi1)
{
  int ica;
  ca_t *ca;
  double ang;
  const double x180 = 180.0 - 1e-12;

  for (ica = 0; ica < pro->nca; ica++) {
    ca = pro->ca + ica;
    if (ca->idcb < 0 || ca->idcg < 0) {
      chi1[ica] = 0.;
      continue;
    }
    ang = (180./M_PI) * rv3_calcdih(NULL,
      x[ca->idc], x[ca->idx], x[ca->idcb], x[ca->idcg], 0);
    if (ang < -x180) ang = -x180; else if (ang > x180) ang = x180;
    //if (ang < 0) ang += 180; else ang -= 180;
    chi1[ica] = ang;
  }
  encode_rotamers(code, pro, chi1);
}

/* scan trajectory and accumulate data */
static int dotrj(protop_t *pro, rotdb_t *rot,
    hist_t *hs_chi1,
    trj_t *trj, const char *fxtc,
    int nat, matrix box)
{
  int nfr, i, match, natoms;
  int fpxtc = 0; /* gromacs file handle for trj */
  rvec *x = NULL;
  double t1, w, *chi1;
  unsigned rotcode[RBMAX] = {0};
  real t = 0.f;

  xnew(chi1, pro->nca);

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

    /* compute chi1 angle */
    getchi1(rotcode, pro, x, chi1);
    if (w > 0.001) /* set a cutoff */
      rotdb_add(rot, rotcode, w);
    for (i = 0; i < pro->nca; i++) {
      if (pro->ca[i].idcg < 0 || pro->ca[i].idcb < 0) continue;
      hs_add1(hs_chi1, i+1, chi1[i], w, HIST_VERBOSE);
      hs_add1(hs_chi1, 0,   chi1[i], w, HIST_VERBOSE);
    }
  }
  free(chi1);
  return 0;
}

/* calculate entropy */
int calc_entropy(hist_t *hs, protop_t *pro, rotdb_t *rot)
{
  int ica, i, len;
  double x, y, sm, p[3], ent0 = 0., ent1 = 0., *arr;

  /* entropy as each residue as independent */
  for (ica = 0; ica < pro->nca; ica++) { /* dihedrals */
    if (pro->ca[ica].idcb < 0 || pro->ca[ica].idcg < 0)
      continue;
    if (pro->ca[ica].gid < 0) continue;
    p[0] = p[1] = p[2] = 0.;
    arr = hs->arr + (ica+1)*hs->n;
    for (i = 0; i < hs->n; i++) { /* histogram bins: angles */
      x = XMIN + (i+.5)*XDEL;
      y = arr[i];
      p[ang2int(x)] += y;
    }
    sm = p[0] + p[1] + p[2];
    for (i = 0; i < 3; i++) {
      p[i] /= sm;
      if (p[i] < 1e-6) continue;
      ent0 += -log(p[i])*p[i];
    }
  }
  /* entropy from rotamer enumeration */
  for (sm = ent1 = 0., i = 0; i < rot->n; i++) {
    x = rot->arr[i].w;
    if (x < 1e-16) continue;
    sm += x;
    ent1 += - x*log(x);
  }
  ent1 = log(sm) + ent1/sm;
  len = rot->len;
  printf("%d, ent0 = %g, %g (%g, %g), ent1 = %g, %g (%g, %g), free = %g\n",
      len, ent0, exp(ent0), ent0/len, exp(ent0/len),
      ent1, exp(ent1), ent1/len, exp(ent1/len), len*log(3.));
  return 0;
}



/* print a summary information */
int print_summary(hist_t *hs, protop_t *pro, double *xref, const char *fn)
{
  int ica, i, iref, imax;
  double x, y, dif, difsum, absum, p[3], ent, *arr;
  double cfrac = 0., cfrac1 = 0., cnt = 0;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  for (ica = 0; ica < pro->nca; ica++) { /* dihedrals */
    fprintf(fp, "%4d %c ", ica+1, aa2letter(pro->ca[ica].resnm));
    if (pro->ca[ica].idcb < 0 || pro->ca[ica].idcg < 0) {
      fprintf(fp, "\n");
      continue;
    }
    p[0] = p[1] = p[2] = 0.;
    arr = hs->arr + (ica+1)*hs->n;
    for (difsum = absum = 0, i = 0; i < hs->n; i++) { /* histogram bins: angles */
      x = XMIN + (i+.5)*XDEL;
      y = arr[i];
      p[ang2int(x)] += y;

      dif = x - xref[ica];
      if (dif > 180) dif -= 360; else if (dif < -180) dif += 360;
      difsum += dif*y;
      absum += fabs(dif)*y;
    }
    x = p[0] + p[1] + p[2];
    difsum /= x;
    absum /= x;
    for (i = 0; i < 3; i++)
      p[i] /= x;
    for (ent = 0, i = 0; i < 3; i++) {
      if (p[i] < 1e-12) continue;
      ent += -p[i]*log(p[i]);
    }
    x = xref[ica];
    iref = ang2int(x);
    fprintf(fp, "%6.4f %6.4f %6.4f %10.6f %d %10.6f %10.6f %10.6f\n",
        p[0], p[1], p[2], xref[ica], iref, difsum, absum, ent);
    for (imax = 0, i = 1; i < 3; i++)
      if (p[i] > p[imax])
        imax = i;
    if (iref == imax) cfrac1 += 1.;
    cfrac += p[iref];
    cnt += 1.;
  }
  printf("correct fraction %g, %g, cnt %g\n", cfrac/cnt, cfrac1/cnt, cnt);
  fprintf(fp, "# %g %g %g\n", cfrac/cnt, cfrac1/cnt, cnt);
  fclose(fp);
  return 0;
}

int run(void)
{
  gmx_mtop_t *mtop;
  int i, ndata, nat;
  char fxtc[FILENAME_MAX];
  char *ftop, *fidx = NULL;
  matrix box;
  rvec *xref = NULL;
  zedata_t *ze = NULL;
  trj_t *trj = NULL;
  double dt, bref = 1.0/(BOLTZ*Tref);
  hist_t *hs_chi1 = NULL;
  double *chi1ref;
  protop_t *pro;
  rotdb_t *rot;
  unsigned rotref[RBMAX];

  /* read topology */
  die_if ((ftop = nfexists(fntop, 7)) == NULL, "no topology file %s\n", fntop);
  mtop = read_tpx_conf(ftop, &xref, box, &nat, &dt);
  if (fnidx) fidx = nfexists(fnidx, 7);

  pro = pro_init(fidx, mtop, xref, PROTOP_SEPHELIX|PROTOP_HYDROPHOBIC|PROTOP_ONLYCHI1);
  xnew(chi1ref, pro->nca);
  getchi1(rotref, pro, xref, chi1ref);
  //rotcode_print(stderr, rotref, pro->ngat); getchar();
  rot = rotdb_init(pro->ngat, rotref);
  hs_chi1 = hs_open(pro->nca+1, XMIN, XMAX, XDEL);


  if (fndist2 != NULL) {
    die_if (0 != hs_load(hs_chi1, fndist, HIST_VERBOSE),
        "cannot load master %s\n", fndist);
    die_if (0 != hs_load(hs_chi1, fndist2, HIST_VERBOSE|HIST_ADDITION),
        "cannot add second %s\n", fndist2);
  }

  if (fndb2 != NULL) {
    die_if (0 != rotdb_load(rot, fndb),
        "cannot load master db %s\n", fndb);
    die_if (0 != rotdb_load(rot, fndb2),
        "cannot add second db %s\n", fndb2);
  }

  if (fndist2 == NULL && fndb2 == NULL) {
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
      if (0 != dotrj(pro, rot,
            hs_chi1, trj, fxtc, nat, box)) {
        fprintf(stderr, "error doing %s\n", fnxtc);
        return -1;
      }
    }
  }
  hs_save(hs_chi1, fndist, HIST_ADDAHALF|HIST_VERBOSE);
  print_summary(hs_chi1, pro, chi1ref, fnout);
  rotdb_trim(rot);
  rotdb_write(rot, fndb);
  calc_entropy(hs_chi1, pro, rot);
  rotdb_free(rot);
  hs_close(hs_chi1);
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


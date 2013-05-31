#ifndef ZTOPUTIL_H__
#define ZTOPUTIL_H__

/* common routines of reading GROMACS topology file
 * pro = pro_init(fidx, mtop, xref, flags);
 *   pro = pro_index(mtop);
 * id = getpairindex(i, j, n);
 * parsepairindex(id, n, &i, &j);
 * */

#define ZCOM_PICK
#define ZCOM_RV3
#define ZCOM_UTIL
#include "zcom.h"

typedef struct {
  int idx; /* gromacs atom index of CA */
  int idn; /* index of N */
  int idc; /* index of C in C=O */
  int idcb; /* index of CB */
  int idcg; /* index of CG */
  int resid; /* residue index */
  char *resnm; /* residue name */
  int gid; /* group index */
} ca_t;

typedef struct {
  ca_t *ca;
  int nca;
  int ngrp, ngg;
  int ngat; /* number of atoms included in groups */
} protop_t;

/* add backbone index-set from moltype imt
 * also search for coorresponding N and C */
static int pro_addca(protop_t *pro, gmx_mtop_t *mtop, int imt, int idx,
    int resid, char *resnm)
{
  int ag, im, imb, nmb = mtop->nmolblock;
  gmx_molblock_t *mb;
  int i, rid2, idn = -1, idc = -1, idcb = -1, idcg = -1;
  gmx_moltype_t *mt = mtop->moltype + imt;
  char *atnm;

  /* search for corresponding N and C */
  for (i = 0; i < mt->atoms.nr; i++) {
    rid2 = mt->atoms.atom[i].resnr;
    if (rid2 != resid) continue;
    atnm = mt->atoms.atomname[i][0];
    if (strcmp("N", atnm) == 0)
      idn = i;
    else if (strcmp("C", atnm) == 0)
      idc = i;
    else if (strcmp("CB", atnm) == 0)
      idcb = i;
    else if (strcmp("CG", atnm) == 0 || strcmp("CG1", atnm) == 0
          || strcmp("OG", atnm) == 0 || strcmp("OG1", atnm) == 0
          || strcmp("SG", atnm) == 0)
      idcg = i;
  }
  if (idn < 0 || idc < 0) {
    fprintf(stderr, "cannot find N(%d) and C(%d) for residue %d, %d\n",
        idn, idc, resid, idx);
    return -1;
  }

  /* search molblock for those whose moltype is imt */
  for (ag = 0, imb = 0; imb < nmb; imb++) {
    mb = mtop->molblock + imb;
    if (mb->type == imt) {
      for (im = 0; im < mb->nmol; im++) {
        xrenew(pro->ca, pro->nca+1);
        pro->ca[pro->nca].idx  = ag + im*mb->natoms_mol + idx;
        pro->ca[pro->nca].idn  = ag + im*mb->natoms_mol + idn;
        pro->ca[pro->nca].idc  = ag + im*mb->natoms_mol + idc;
        pro->ca[pro->nca].idcb = ag + im*mb->natoms_mol + idcb;
        pro->ca[pro->nca].idcg = ag + im*mb->natoms_mol + idcg;
        pro->ca[pro->nca].resid = resid;
        pro->ca[pro->nca].resnm = resnm;
        pro->ca[pro->nca].gid = -1; /* invalid gid */
        pro->nca++;
      }
    }
    ag += mb->nmol * mb->natoms_mol;
  }
  return 0;
}

/* obtain C-alpha indices */
static protop_t *pro_index(gmx_mtop_t *mtop)
{
  gmx_moltype_t *mt;
  int imt, i, resid;
  protop_t *pro;
  char *resnm;

  xnew(pro, 1);
  pro->nca = 0;
  pro->ngrp = 0;
  xnew(pro->ca, 1);
  /* loop over moltypes to search a interaction that matches `allatoms' */
  for (imt = 0; imt < mtop->nmoltype; imt++) {
    mt = mtop->moltype + imt;
    for (i = 0; i < mt->atoms.nr; i++) {
      if (strcmp("CA", mt->atoms.atomname[i][0]) != 0)
        continue;
      resid = mt->atoms.atom[i].resnr;
      resnm = mt->atoms.resname[resid][0];
      if (0 != pro_addca(pro, mtop, imt, i, resid, resnm)) {
        free(pro->ca); free(pro);
        return NULL;
      }
    }
  } /* loop over moltypes */
  printf("C-alpha atoms:\n");
  for (i = 0; i < pro->nca; i++)
    printf("%d ", pro->ca[i].idx+1);
  printf("\nTotal: %d atoms\n", pro->nca);
  return pro;
}

/* find an ca for a residue resid */
static __inline int pro_findres(protop_t *pro, int resid)
{
  int i;
  for (i = 0; i < pro->nca; i++)
    if (pro->ca[i].resid == resid)
      return i;
  return -1;
}

/* recompute ngat */
static __inline void pro_recount(protop_t *pro)
{
  int i, cnt = 0;
  for (i = 0; i < pro->nca; i++)
    if (pro->ca[i].gid >= 0) cnt++;
  pro->ngat = cnt;
}

/* set gid of nonhydrophobic residues to -1 */
static void pro_hydrophobic(protop_t *pro)
{
#define HBCNT 5
  int i, j, n;
  const char *reshb[HBCNT] = {"LEU", "VAL", "ILE", "TRP", "PHE"};
  char *rnm;

  for (i = 0; i < pro->nca; i++) {
    rnm = pro->ca[i].resnm;
    n = strlen(rnm);
    if (n == 4) rnm++;
    for (j = 0; j < HBCNT; j++) {
      if (strcmp(rnm, reshb[j]) == 0) break;
    }
    if (j == HBCNT) /* not hydrophobic */
      pro->ca[i].gid = -1;
  }
  pro_recount(pro);
}

/* set gid of ala and gly residues to -1 */
static void pro_onlychi1(protop_t *pro)
{
#define SMLCNT 5
  int i, n;
  char *rnm;

  for (i = 0; i < pro->nca; i++) {
    rnm = pro->ca[i].resnm;
    n = strlen(rnm);
    if (n == 4) rnm++;
    if (strcmp(rnm, "GLY") == 0 || strcmp(rnm, "ALA") == 0)
      pro->ca[i].gid = -1;
  }
  pro_recount(pro);
}

/* load group index from file fn, initialize atom index
 * format:
 * # number_of_groups
 * gid0 : resid1 resid2
 * gid1 : resid3 resid4-resid5
 * gid starts from 0
 * resid starts from 1, can also has a prefix like LEU, or L
 * */
static int pro_loadgroups(protop_t *pro, const char *fn)
{
  FILE *fp;
  char s[512], item[8], *p, *q, *pivot;
  int ig, j, cid, x0, x, x1, next, na = 0, ng;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  /* read the first line to get the number of atom groups */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "no first line %s\n", fn);
    fclose(fp);
    return -1;
  }
  if (s[0] != '#') {
    fprintf(stderr, "%s: missing first line #\n", fn);
    goto ERR;
  }
  p = s+1;
  next = 0;
  if (1 != sscanf(p, "%d%n", &ng, &next)) {
    fprintf(stderr, "%s: cannot get the # of atom groups\n", fn);
    goto ERR;
  }
  pro->ngrp = ng;
  pro->ngg = ng*(ng-1)/2;

  /* read residue indices */
  for (ig = 0; ig < ng; ig++) {
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "no atom index for group %d\n", ig);
      goto ERR;
    }
    next = 0;
    p = s;
    if (1 != sscanf(p, "%d : %n", &x, &next)) {
      fprintf(stderr, "no id in group %d\n", ig);
      goto ERR;
    }
    if (x != ig) {
      fprintf(stderr, "index mismatch got %d, should be %d\n", x, ig);
      goto ERR;
    }
    p += next;
    for (j = 0; ; p += next) {
      if (1 != sscanf(p, "%s%n", item, &next))
        break; /* out of numbers  */
      /* determine if it is range or a number */
      if ((pivot = strchr(item, '-')) != NULL) { /* range */
        *pivot = '\0';
        for (q = item; isalpha(*q); q++) ;
        x0 = x1 = -1;
        sscanf(q, "%d", &x0);
        for (q = pivot+1; isalpha(*q); q++) ;
        sscanf(q, "%d", &x1);
      } else { /* single number */
        for (q = item; isalpha(*q); q++) ;
        x0 = -1;
        sscanf(q, "%d", &x0); /* residue id */
        x1 = x0;
      }
      for (x = x0; x <= x1; x++) {
        cid = pro_findres(pro, x - 1); /* index - 1 convention */
        if (cid < 0) {
          fprintf(stderr, "cannot find ca for residue %d\n", x);
          goto ERR;
        }
        pro->ca[cid].gid = ig;
        //printf("%s%d ", pro->ca[cid].resnm, x);
      }
      j += x1 - x0 + 1;
    }
    //printf("group %d has %d\n", ig, j);
    na += j;
  }
  fclose(fp);
  pro->ngat = na;
  printf("%s loaded, %d atoms, %d groups\n", fn, na, ng);
  return 0;
ERR:
  fclose(fp);
  return -1;
}

/* make groups according helical conformation in xref */
static void pro_counthelix(protop_t *pro, rvec *xref)
{
  int i, j, nse = 0, is, it, nres = pro->nca;
  int *se, *ishx, quin[5];
  double phi, psi;

  /* A. make an array of nres, identify if each residue is helix */
  xnew(ishx, nres);
  ishx[0] = ishx[nres-1] = 0;
  for (i = 1; i < nres-1; i++) {
    /* make local quintet */
    quin[0] = pro->ca[i-1].idc;
    quin[1] = pro->ca[i].idn;
    quin[2] = pro->ca[i].idx;
    quin[3] = pro->ca[i].idc;
    quin[4] = pro->ca[i+1].idn;
    phi = rv3_calcdihv(NULL, xref, quin, 0);
    psi = rv3_calcdihv(NULL, xref, quin+1, 0);
    ishx[i] = (phi < 0 && psi > -100*M_PI/180 && psi < 80*M_PI/180);
  }

  /* B. searching for segments
   * make 2*pro->ngrp for start/end of each segment
   * range of segment k is se[2*k] <= id < se[2*k+1] */
  xnew(se, 2);
  for (i = 0, is = 0; i < nres; ) { /* try to find the helices */
    while (i < nres && !ishx[i]) i++;
    if (i >= nres) break; /* no more helices */
    is = i;
    while (ishx[i] && i < nres) i++;
    it = i;
    for (; is < it; is++) { /* skip terminal GLY and PRO */
      if (strcmp(pro->ca[is].resnm, "GLY") != 0
       && strcmp(pro->ca[is].resnm, "PRO") != 0) break;
    }
    for (; it > is; it--) {
      if (strcmp(pro->ca[it-1].resnm, "GLY") != 0
       && strcmp(pro->ca[it-1].resnm, "PRO") != 0) break;
    }
    if (it - is >= 4) { /* successfully find a helical segment */
      xrenew(se, 2*(nse+1));
      se[2*nse] = is;
      se[2*nse+1] = it;
      nse++;
      printf("found helix %d: from %d to %d\n", nse, is+1, it);
    } else { } /* just let go, don't increment nse */
  }
  pro->ngrp = nse;
  pro->ngg = pro->ngrp*(pro->ngrp-1)/2;

  /* C. fill ca[x].gid from se/nse */
  for (pro->ngat = 0, i = 0; i < nse; i++) {
    for (j = se[2*i]; j < se[2*i+1]; j++) {
      if (j >= pro->nca) {
        fprintf(stderr, "j %d >= cacnt %d\n", j, pro->nca);
        exit(1);
      }
      pro->ca[j].gid = i;
      pro->ngat++;
    }
  }
  free(se);
  free(ishx);
}

/* print group information */
static void pro_printgroup(const protop_t *pro)
{
  int ig, i, cnti = 0, cnt = 0;

  for (ig = 0; ig < pro->ngrp; ig++) {
    printf("group %2d: ", ig);
    cnti = 0;
    for (i = 0; i < pro->nca; i++) {
      if (pro->ca[i].gid == ig) {
        printf("%s%-3d ", pro->ca[i].resnm, i+1);
        cnti++;
        if (cnti % 10 == 0) printf("\n          ");
      }
    }
    cnt += cnti;
    printf("\n");
  }
  printf("%d groups, %d atoms\n", pro->ngrp, cnt);
}

#define PROTOP_SEPHELIX    0x1000
#define PROTOP_HYDROPHOBIC 0x2000
#define PROTOP_ONLYCHI1    0x4000

/* load ca index and partition
 * load groups from fidx if it is not NULL
 * or init groups if PROTOP_SEPHELIX is given
 * */
static protop_t *pro_init(const char *fidx, gmx_mtop_t *mtop, rvec *xref,
    unsigned flags)
{
  protop_t *pro;

  /* compute ca index */
  pro = pro_index(mtop);
  if (fidx != NULL) {
    pro_loadgroups(pro, fidx);
  } else if (flags & PROTOP_SEPHELIX) {
    pro_counthelix(pro, xref);
    if (flags & PROTOP_HYDROPHOBIC) {
      pro_hydrophobic(pro); /* set gid of non-hydrophobic residues to -1 */
    }
    printf("Parse helices, %d atoms; %d helices\n", pro->nca, pro->ngrp);
  }
  if (flags & PROTOP_ONLYCHI1) {
    pro_onlychi1(pro); /* set gid of ALA & GLY to -1 */
  }
  pro_recount(pro);
  pro_printgroup(pro);
  return pro;
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

#endif


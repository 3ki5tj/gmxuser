/*
  analyzer for backbone dihedral angles and their distributions

  Copyright (C) 2010  Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#define ZCOM_PICK
#define ZCOM_LU
#define ZCOM_RV3
#define ZCOM_SS  /* safe string */
#include "zcom.h"

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)
#define BOLTZ   0.008314511212
#define ROOMT   300.0
#define KT (BOLTZ*ROOMT)

/* global parameters */
int analysis_level=0;
int skip_pro=1,skip_gly=1;
int verbose=0;
char *listname=NULL;
char *fname_input=NULL;
char *fname_ref=NULL;
int average_bins=2;
int makewei_bins=2;
double epsilon=1e-4;
double demon_factor=2.0;
double exclude_last_sine=0;
int grouping_neis=1;
#define CSORDER 32
int cs_order=6;
int cartoon_style=0;

static void draw_cartoon0(int phi, int psi){
  const char h1[3][8] = {"**    ", "||    ", "~~~~~~"};
  const char h2[3][8] = {"    **", "    ||", "~~~~~~"};
  printf("%s%s", h1[phi], h2[psi]);
}

static void draw_cartoon1(int phi, int psi){
  if((phi!=psi) || phi==0){
    printf("     **     ");
  }else if(phi==2){
    printf("~~~~~~~~~~~~");
  }else{
    printf("||        ||");
  }
}

static void draw_cartoon(int phi, int psi)
{
  if (cartoon_style == 0)
    draw_cartoon0( phi, psi );
  else
    draw_cartoon1( phi, psi );
}

// polynominal expansion of the reference force field
// AMBER03
double cp_phi03[2*CSORDER]={8.08349,-1.41503,-2.88780,-3.78066};
double cp_psi03[2*CSORDER]={16.95692,-2.93131,-12.16456,7.72366};
// AMBER99SB
double cp_phisb[2*CSORDER]={1.75724,5.27184,2.25936,-7.02912};
double cp_psisb[2*CSORDER]={17.40544,-5.0208,-13.22144,9.20480};


double *cp_phiref=cp_phi03;
double *cp_psiref=cp_psi03;
double cs_phiref[2*CSORDER]; // Fourier coefficients
double cs_psiref[2*CSORDER];

// standard cosine and sine values
double c36,s36, c12,s12, c72,s72, c30,s30;

// a normal vector for the area A-B-C
static double *areavec(double w[], double a[], double b[], double c[])
{
  double u[3], v[3];
  rv3_diff(u, b, a);
  rv3_diff(v, c, b);
  return rv3_cross(w, u, v);
}

#define normareavec(w, a, b, c) rv3_normalize( areavec(w, a, b, c) )

// remove leading/tailing spaces
static char *trim(char *s)
{
  char *p=s+strlen(s)-1, *q;
  // trim the tail
  while (isspace(*p) && p >= s)
    *p-- = '\0';

  // trim the leading spaces
  // first nonspace place
  for(q=s; *q && isspace(*q); q++) ;
  // shifting
  for(p=s; *q; ) *p++=*q++;
  if(*q=='\0') *p='\0';
  return s;
}

enum AATYPE{ ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE,
  LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL};
char aaletter[64]="ARNDCEQGHILKMFPSTWYV";

// return the integer index
static int resname2itype(const char *name){
  const char *p=name;

  // be careful with 4-letter residues
  if(strcmp(name, "CYS2")==0){
    return CYS;
  }else if(name[3] != '\0' && (name[0]=='C'||name[0]=='N'))
    p++; // skip the terminal letter 'C' or 'N'

  if      (strcmp (p, "ALA")  == 0)  return ALA;
  else if (strcmp (p, "ARG")  == 0)  return ARG;
  else if (strcmp (p, "ASN")  == 0)  return ASN;
  else if (strcmp (p, "ASP")  == 0)  return ASP;
  else if (strncmp(p, "CY",2) == 0)  return CYS;
  else if (strcmp (p, "GLU")  == 0)  return GLU;
  else if (strcmp (p, "GLN")  == 0)  return GLN;
  else if (strcmp (p, "GLY")  == 0)  return GLY;
  else if (strncmp(p, "HI",2) == 0)  return HIS;
  else if (strcmp (p, "ILE")  == 0)  return ILE;
  else if (strcmp (p, "LEU")  == 0)  return LEU;
  else if (strncmp(p, "LY",2) == 0)  return LYS;
  else if (strcmp (p, "MET")  == 0)  return MET;
  else if (strcmp (p, "PHE")  == 0)  return PHE;
  else if (strcmp (p, "PRO")  == 0)  return PRO;
  else if (strcmp (p, "SER")  == 0)  return SER;
  else if (strcmp (p, "THR")  == 0)  return THR;
  else if (strcmp (p, "TRP")  == 0)  return TRP;
  else if (strcmp (p, "TYR")  == 0)  return TYR;
  else if (strcmp (p, "VAL")  == 0)  return VAL;
  return -1;
}


#define NM 20480

typedef struct tag_residue{
  double rn[3];
  double rca[3];
  double rc[3];
  char resname[8];
  int resid;
  /* additional information */
  double phi, psi;
  double aphi, apsi;
  int phi_type, psi_type;
}residue_t;

// atom name should always have 4 letters, the next character is the alternative location
#define GET_ATNAME(name, line) alter_atname( trim( memcpy(name, line+12, 4) ), 0)
// allow residue name to have 4 letters
#define GET_RESNAME(name, line) trim( memcpy(name, line+17, 4) )
// residue id
#define GET_RESID(id, line) id=-1; sscanf(line+22, "%d", &id)
// get coordinates
#define GET_XYZ(r, line) sscanf(buf+30, "%lf%lf%lf",&(r[0]),&(r[1]),&(r[2]))

// handle alternative form of atom names, like 1HD1, 2HD1, 3HD1, instead of HD11, HD12 and HD13, etc
// flag=0 to normal form, flag=1 to alternative form
static char *alter_atname(char *atname, int flag){
  char atname2[8]=" ";
  if(flag==0){ // normal form
    if(isdigit(atname[0])){
      sprintf(atname2, "%s%c", atname+1, atname[0]);
      strcpy(atname, atname2); // copy back
    }
  }else{ // to alternative form
    if(!isdigit(atname[0])){
      int size=strlen(atname);
      if(!isdigit(atname[size-1])) return atname;
      atname2[0]=atname[size-1];
      atname[size-1]='\0';
      strcpy(atname2+1, atname);
      strcpy(atname, atname2); // copy back
      //printf("alternative form is [%s]\n", atname);
    }
  }
  return atname;
}

#define DANG 1
/* #define ANGCNT ((int)(360.0/DANG-1e-6)+1) */
#define ANGCNT 360


int Hphi[ANGCNT];
int Hpsi[ANGCNT];
int Hdihs[ANGCNT][ANGCNT];
double dist_phi[ANGCNT], wei_phi[ANGCNT];
double dist_psi[ANGCNT], wei_psi[ANGCNT];
double phiA[ANGCNT], phiB[ANGCNT], wei_phiA[ANGCNT], wei_phiB[ANGCNT];
double psiA[ANGCNT], psiB[ANGCNT], wei_psiA[ANGCNT], wei_psiB[ANGCNT];

static int load_backbone(residue_t *seq, const char *fname)
{
  FILE *fp;
  char *buf=NULL, resname[8], atname[8];
  int i, resid;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s\n", fname);
    return -1;
  }

  resname[4] = atname[4] = '\0';
  for (i = -1; ssfgets(buf, NULL, fp); ) {
    // we don't need multiple models
    if(strncmp(buf, "ENDMDL", 6) == 0) break;

    // we are not interested in HEMATM, ANISOU, etc ...
    if(strncmp(buf, "ATOM  ", 6) != 0) continue;
    // get rid of alternative locations
    if(buf[16] != ' ' && buf[16]!='A') continue;
    GET_ATNAME(atname, buf);
    GET_RESNAME(resname, buf);
    GET_RESID(resid, buf);
    if (resname2itype(resname) < 0)
      continue;

    // math the previous residue
    if( i<0 || seq[i].resid != resid ){
      i++;
      if(i>=NM){
        fprintf(stderr, "too many residues from %s\n", fname);
        i--;
        break;
      }
      seq[i].resid =resid;
      strcpy(seq[i].resname, resname);
    }

    if (strcmp(atname, "N") == 0)
      GET_XYZ(seq[i].rn, buf);
    else if (strcmp(atname, "CA") == 0)
      GET_XYZ(seq[i].rca, buf);
    else if (strcmp(atname, "C") == 0)
      GET_XYZ(seq[i].rc, buf);
  }

  fclose(fp);
  return i+1;
}

#define LOOP  0
#define ALPHA 1
#define BETA  2
static int computedihs(residue_t seq[], int i, int nres)
{
  double dot, w[3], p[3], q[3];
  double phi, psi;
  int phi_type, psi_type;

  /* calculate phi/psi angle, no pbc */
  phi = psi = 0.0;
  phi_type = psi_type = LOOP;
  normareavec(q, seq[i].rn, seq[i].rca, seq[i].rc);

  /* calculate phi angle, if possible, trans is zero
     phi = C(previous)-N-CA-C */
  if(i > 0 && seq[i-1].resid + 1 == seq[i].resid) {
              /* make sure the sequence continuity */
    normareavec(p, seq[i].rca, seq[i].rn, seq[i-1].rc);
    dot = rv3_dot(p, q);
    phi = acos(dot);
    rv3_diff(w, seq[i-1].rc, seq[i].rn);
    dot = rv3_dot(w, q);
    if(dot < 0)
      phi = -phi;
    if(phi > M_PI*80/180)
      phi_type = ALPHA;
    else if(phi > -M_PI*50/180)
      phi_type = BETA;
    else
      phi_type = LOOP;
  }

  /* calculate psi angle,
     psi = N-CA-C-N(next) */
  if (i < nres-1 && seq[i].resid + 1 == seq[i + 1].resid) {
    normareavec(p, seq[i+1].rn, seq[i].rc, seq[i].rca);
    dot = rv3_dot(p, q);
    psi = acos(dot);
    rv3_diff(w, seq[i].rc, seq[i+1].rn);
    dot = rv3_dot(w, q);
    if(dot < 0)
      psi = -psi;
    if (psi > M_PI/3 || psi < -M_PI*160/180)
      psi_type = ALPHA;
    else if (psi >-M_PI*2.0/3)
      psi_type = BETA;
    else
      psi_type = LOOP;
  }

/*
    printf("%3d: RES=[%4s] N=%8.3f%8.3f%8.3f, CA=%8.3f%8.3f%8.3f, C=%8.3f%8.3f%8.3f phi=%8.3f psi=%8.3f\n", i, seq[i].resname,
        seq[i].rn[0], seq[i].rn[1], seq[i].rn[2],
        seq[i].rca[0],seq[i].rca[1],seq[i].rca[2],
        seq[i].rc[0], seq[i].rc[1], seq[i].rc[2],
        phi*RAD2DEG, psi*RAD2DEG);
*/

  /* we don't do it in style 0, which is supposed to show
   * individual phi/psi angles */
  if (cartoon_style == 1) {
    /* since for phi, the division between alpha and beta is unclear
     * we try to improve from psi */
    if (psi_type == ALPHA && phi_type == BETA && phi > M_PI*20/180)
      phi_type = ALPHA;
    else if (psi_type == BETA && phi_type == ALPHA && phi < M_PI*140/180)
      phi_type = BETA;
  }

  seq[i].phi = phi;
  seq[i].psi = psi;
  seq[i].aphi = phi*RAD2DEG;
  seq[i].apsi = psi*RAD2DEG;
  seq[i].phi_type = phi_type;
  seq[i].psi_type = psi_type;
  return 0;
}

/* from a single file (possible with a reference)
 * load coordinates of backbone atoms and compute dihedrals
 * results are saved to `seq'
 * */
static int loadxyz(const char *fname, const char *fname2)
{
  static residue_t seq[NM], seq2[NM];
  char lab[8] = "LAB", *p;
  char buf[64], title[8];
  int nres, nres2;
  int i, iphi, ipsi;
  int itype, histable;

  /* we only keep four letters for the identification */
  strncpy(title, fname, 4);
  title[4]='\0';
  if ((p=strchr(title, '.')) != NULL)
    *p = '\0';

  /* load backbone coordinates */
  nres = load_backbone(seq, fname);
  nres2 = (fname2 != NULL) ? load_backbone(seq2, fname2) : nres;

  if (nres <= 0 || nres2 <= 0) {
    fprintf(stderr, "error during loading residues, nres=%d, nres2=%d.\n",
        nres, (fname2 != NULL) ? 0: nres2);
    return 1;
  } else if (nres != nres2) {
    fprintf(stderr, "reference # of residues: %d, should be %d.\n", nres2, nres);
    return 1;
  }

  for (i = 0; i < nres; i++) {
    itype = resname2itype(seq[i].resname);

    /* skip proline and glycine */
    if ( (itype == PRO && skip_pro)
      || (itype == GLY && skip_gly)) {
      histable = 0;
    } else {
      histable = 1;
    }

    computedihs(seq, i, nres);
    if (fname2 != NULL)
      computedihs(seq2, i, nres2);

    // accumulate histograms
    if (seq[i].phi != 0.0 && seq[i].psi != 0.0 && histable) {
      iphi = (int)((seq[i].aphi + 180.0)/DANG);
      ipsi = (int)((seq[i].apsi + 180.0)/DANG);
      if (iphi >= 0 && iphi < ANGCNT && ipsi >= 0 && ipsi < ANGCNT)
        Hdihs[iphi][ipsi]++;
    }

    if (verbose > 0) {
      sprintf(buf, "[%s]", seq[i].resname);
      printf("%-4s %3d(%3d): RES=%-6s %c",
        title, i+1, seq[i].resid, buf, histable ? ' ' : seq[i].resname[0]);

      printf(" phi=%8.3f[%c] psi=%8.3f[%c] ",
        seq[i].aphi, lab[ seq[i].phi_type], seq[i].apsi, lab[ seq[i].psi_type ]);
      draw_cartoon(seq[i].phi_type, seq[i].psi_type);

      if (fname2 != NULL) {
        printf("   -&-   ");
        draw_cartoon(seq2[i].phi_type, seq2[i].psi_type);
        printf(" phi=%8.3f[%c] psi=%8.3f[%c] ",
          seq2[i].aphi, lab[ seq2[i].phi_type], seq2[i].apsi, lab[ seq2[i].psi_type ]);
      }
      printf("\n");
    }
  }

  return 0;
}

/*
// fine tuning the distribution,
// further average for rare data
static double *disttune(double dist[], int del, double low){
  int i,j,jwrap,cnt;
  static double tmp[ANGCNT+1];
  double max=0.0,x,sum;
  for(i=0; i<ANGCNT; i++){
    x=dist[i];
    if(x>max) max=x;
  }
  low *= max;
  for(sum=0.0, i=0; i<ANGCNT; i++){
    if(dist[i]>low){
      tmp[i]=dist[i];
    }else{
      cnt=0; x=0.0;
      for(j=i-del; j<=i+del; j++){
        jwrap=(j+ANGCNT*100)%ANGCNT;
        x+=dist[jwrap];
        cnt++;
      }
      tmp[i]=x/cnt;
      if(i<30){
        fprintf(stderr, "i=%d old value %g, new value %g\n", i, dist[i], tmp[i]);
      }
    }
    sum+=tmp[i];
  }
  sum=1.0/(sum*DANG);
  for(i=0; i<ANGCNT; i++) dist[i]=tmp[i]*sum;
  return dist;
}
*/
// smoothen and convert histogram to distribution
static double *hist2dist(double dist[], int hist[], int del, int pbc, double *eps)
{
  int i, j, jwrap,cnt;
  double x, tot, max;

  for (tot = 0.0, i = 0; i < ANGCNT; i++) {
    /* average over neighboring bins */
    for(x = 0.0, cnt = 0, j = i-del; j < i+del; j++){
      if (pbc) {
        jwrap = (j+ANGCNT*100)%ANGCNT;
        x += hist[jwrap];
        cnt++;
      }else if (j >= 0 && j < ANGCNT) {
        x += hist[j];
        cnt++;
      }
    }
    dist[i] = x/cnt;
    tot += dist[i];
  }
  tot = 1.0/tot/DANG;
  for (i = 0; i < ANGCNT; i++)
    dist[i] *= tot;

  /* normalize */
  for (max = 0.0, i = 0; i < ANGCNT; i++) {
    if ((x=fabs(dist[i])) > max) {
      max = x;
    }
  }
  for(i=0; i<ANGCNT; i++){
    // distribution must be positive
    x=dist[i];
    if(x<max*epsilon) dist[i]=max*epsilon;
  }
  if(eps) *eps=max*epsilon;

  return dist;
}



static void demonize(double w[]){
  int i;
  double max=0.0,eps=1e-2;
  for(i=0; i<ANGCNT; i++){
    w[i]=pow(w[i], demon_factor);
    if(w[i]>max)max=w[i];
  }
  for(i=0; i<ANGCNT; i++){
    w[i]/=max;
    if(w[i]<eps) w[i]=eps;
  }
}

// obtained by e.g., TrigExpand[Cos[5 x]], or TrigExpand[Sin[6 x]/Sin[x]] WolframAlpha
// which is just ChebyshevT[5, Cos[x]], for the cosine part
// currently to 12th order
static int cs2cp(double pc[], double fc[]){
  int j;
  double x;
  double *ps, *fs;
  for(j=0; j<2*cs_order; j++) pc[j]=0.0;
  if(cs_order>0){ pc[0]=fc[0]; }
  if(cs_order>1){ pc[1]=fc[1]; }
  if(cs_order>2){ x=fc[2]; pc[0]-=x;   pc[2]+=2*x;}
  if(cs_order>3){ x=fc[3]; pc[1]-=3*x; pc[3]+=4*x;}
  if(cs_order>4){ x=fc[4]; pc[0]+=x;   pc[2]-=8*x;   pc[4]+=8*x;}
  if(cs_order>5){ x=fc[5]; pc[1]+=5*x; pc[3]-=20*x;  pc[5]+=16*x;}
  if(cs_order>6){ x=fc[6]; pc[0]-=x;   pc[2]+=18*x;  pc[4]-=48*x;  pc[6]+=32*x; }
  if(cs_order>7){ x=fc[7]; pc[1]-=7*x; pc[3]+=56*x;  pc[5]-=112*x; pc[7]+=64*x; }
  if(cs_order>8){ x=fc[8]; pc[0]+=x;   pc[2]-=32*x;  pc[4]+=160*x; pc[6]-=256*x; pc[8]+=128*x;}
  if(cs_order>9){ x=fc[9]; pc[1]+=9*x; pc[3]-=120*x; pc[5]+=432*x; pc[7]-=576*x; pc[9]+=256*x;}
  if(cs_order>10){x=fc[10];pc[0]-=x;   pc[2]+=50*x;  pc[4]-=400*x; pc[6]+=1120*x;pc[8]-=1280*x;pc[10]+=512*x; }
  if(cs_order>11){x=fc[11];pc[1]-=11*x;pc[3]+=220*x; pc[5]-=1232*x;pc[7]+=2816*x;pc[9]-=2816*x;pc[11]+=1024*x;}

  ps=pc+cs_order;
  fs=fc+cs_order;
  if(cs_order>0){ ps[0]=fs[0];}
  if(cs_order>1){ ps[1]=2*fs[1];}
  if(cs_order>2){ x=fs[2]; ps[0]-=x;   ps[2]+=4*x;}
  if(cs_order>3){ x=fs[3]; ps[1]-=4*x; ps[3]+=8*x;}
  if(cs_order>4){ x=fs[4]; ps[0]+=x;   ps[2]-=12*x; ps[4]+=16*x;}
  if(cs_order>5){ x=fs[5]; ps[1]+=6*x; ps[3]-=32*x; ps[5]+=32*x;}
  if(cs_order>6){ x=fs[6]; ps[0]-=x;   ps[2]+=24*x; ps[4]-=80*x;  ps[6]+=64*x; }
  if(cs_order>7){ x=fs[7]; ps[1]-=8*x; ps[3]+=80*x; ps[5]-=192*x; ps[7]+=128*x; }
  if(cs_order>8){ x=fs[8]; ps[0]+=x;   ps[2]-=40*x; ps[4]+=240*x; ps[6]-=448*x;  ps[8]+=256*x; }
  if(cs_order>9){ x=fs[9]; ps[1]+=10*x;ps[3]-=160*x;ps[5]+=672*x; ps[7]-=1024*x; ps[9]+=512*x; }
  if(cs_order>10){x=fs[10];ps[0]-=x;   ps[2]+=60*x; ps[4]-=560*x; ps[6]+=1792*x; ps[8]-=2304*x; ps[10]+=1024*x;}
  if(cs_order>11){x=fs[11];ps[1]-=12*x;ps[3]+=280*x;ps[5]-=1792*x;ps[7]+=4608*x; ps[9]-=5120*x; ps[11]+=2048*x;}

  // stupid flip of SIGN due to convention used in md4!
  for(j=0; j<cs_order; j++) ps[j]=-ps[j];

  return (cs_order>12)?1:0;
}

/* obtained by e.g., TrigReduce[Cos[x]^2]
 * input pc
 * output fc
 * */
static int cp2cs(double fc[], double pc[]){
  int j;
  double x;
  double *ps, *fs;
  for(j=0; j<2*cs_order; j++) fc[j]=0.0;
  if(cs_order>0){ fc[0]=pc[0]; }
  if(cs_order>1){ fc[1]=pc[1]; }
  if(cs_order>2){ x=pc[2]; fc[0]+=x*1.0/2;     fc[2]+=x*1.0/2; }
  if(cs_order>3){ x=pc[3]; fc[1]+=x*3.0/4;     fc[3]+=x*1.0/4;}
  if(cs_order>4){ x=pc[4]; fc[0]+=x*3.0/8;     fc[2]+=x*4.0/8;     fc[4]+=x*1.0/8;}
  if(cs_order>5){ x=pc[5]; fc[1]+=x*10.0/16;   fc[3]+=x*5.0/16;    fc[5]+=x*1.0/16;}
  if(cs_order>6){ x=pc[6]; fc[0]+=x*10.0/32;   fc[2]+=x*15.0/32;   fc[4]+=x*6.0/32;    fc[6]+=x*1.0/32;}
  if(cs_order>7){ x=pc[7]; fc[1]+=x*35.0/64;   fc[3]+=x*21.0/64;   fc[5]+=x*7.0/64;    fc[7]+=x*1.0/64;}
  if(cs_order>8){ x=pc[8]; fc[0]+=x*35.0/128;  fc[2]+=x*56.0/128;  fc[4]+=x*28.0/128;  fc[6]+=x*1.0/128;  fc[8]+=x*1.0/128; }
  if(cs_order>9){ x=pc[9]; fc[1]+=x*126.0/256; fc[3]+=x*84.0/256;  fc[5]+=x*36.0/256;  fc[7]+=x*9.0/256;  fc[9]+=x*1.0/256;}
  if(cs_order>10){x=pc[10];fc[0]+=x*126.0/512; fc[2]+=x*210.0/512; fc[4]+=x*120.0/512; fc[6]+=x*45.0/512; fc[8]+=x*10.0/512; fc[10]+=x*1.0/512;}
  if(cs_order>11){x=pc[11];fc[1]+=x*462.0/1024;fc[3]+=x*330.0/1024;fc[5]+=x*165.0/1024;fc[7]+=x*55.0/1024;fc[9]+=x*11.0/1024;fc[11]+=x*1.0/1024;}

  ps=pc+cs_order;
  fs=fc+cs_order;
  if(cs_order>0) fs[0]=ps[0];
  if(cs_order>1) fs[1]=ps[1]*0.5;
  if(cs_order>2){ x=ps[2]; fs[0]+=x*.25;       fs[2]+=x*.25; }
  if(cs_order>3){ x=ps[3]; fs[1]+=x*.25;       fs[3]+=x*.125; }
  if(cs_order>4){ x=ps[4]; fs[0]+=x*.125;      fs[2]+=x*.1875;     fs[4]+=x*.0625;}
  if(cs_order>5){ x=ps[5]; fs[1]+=x*.15625;    fs[3]+=x*.125;      fs[5]+=x/32.0;}
  if(cs_order>6){ x=ps[6]; fs[0]+=x*5.0/64;    fs[2]+=x*9.0/64;    fs[4]+=x*5.0/64;    fs[6]+=x*1.0/64;}
  if(cs_order>7){ x=ps[7]; fs[1]+=x*14.0/128;  fs[3]+=x*14.0/128;  fs[5]+=x*6.0/128;   fs[7]+=x*1.0/128;}
  if(cs_order>8){ x=ps[8]; fs[0]+=x*14.0/256;  fs[2]+=x*28.0/256;  fs[4]+=x*20.0/256;  fs[6]+=x*7.0/256;   fs[8]+=x*1.0/256;}
  if(cs_order>9){ x=ps[9]; fs[1]+=x*42.0/512;  fs[3]+=x*48.0/512;  fs[5]+=x*27.0/512;  fs[7]+=x*8.0/512;   fs[9]+=x*1.0/512;}
  if(cs_order>10){x=ps[10];fs[0]+=x*42.0/1024; fs[2]+=x*90.0/1024; fs[4]+=x*75.0/1024; fs[6]+=x*35.0/1024; fs[8]+=x*9.0/1024; fs[10]+=x*1.0/1024;}
  if(cs_order>11){x=ps[11];fs[1]+=x*132.0/2048;fs[3]+=x*165.0/2048;fs[5]+=x*110.0/2048;fs[7]+=x*44.0/2048; fs[9]+=x*10.0/2048;fs[11]+=x*1.0/2048;}

  // stupid flip of sign
  for(j=0; j<cs_order; j++) fs[j]=-fs[j];
  return (cs_order>12)?1:0;
}

#define GET_CS(bj, j) bj=((j==0)?1.0:((j<cs_order)?cos(j*ang):sin((j-cs_order+1)*ang)))

static double *cs_combine(double dest[], double a[], double b[], double fa, double fb, double offsetb){
  int i;
  double x;
  for(i=0; i<cs_order*2; i++){
    x=b[i];
    if(i==0)x-=offsetb;
    dest[i]=a[i]*fa+x*fb;
  }
  return dest;
}

// fix baseline, input
static double fixbaseline(double a[]){
  int j;
  double min=1e9,ang,y,bj;
  for(ang=-M_PI; ang<M_PI; ang+=0.00001){
    for(y=0,j=0; j<2*cs_order; j++){
      GET_CS(bj, j);
      y+=bj*a[j];
    }
    if(y<min) min=y;
  }
  if(verbose>=0) printf("fixing baseline according to the minimal %g\n", min);
  a[0]-=min;
  return min;
}

#define CSFMT  "%+.7f"
#define CPFMT  "%.7f"
#define CS_FACTOROUT  0x1
#define CS_NOLOG      0x10
#define CS_ARRAY      0x100
#define CS_EXPR       0x200
#define CS_COSPOLY    0x1000

static int cs_print(double a[], int flag, char *title){
  int j;
  int factorout=flag&CS_FACTOROUT;
  int arrayform=flag&CS_ARRAY;
  int exprform=flag&CS_EXPR;
  int cospoly=flag&CS_COSPOLY;
  double max=1.0,x;
  double p[CSORDER*2]={0.0};

  //printf("COS=%d, EXPR=%d, ARRAY=%d, FACT=%d\n", cospoly, exprform, arrayform, factorout);

  // convert cos(n x) to a polynomial of ( cos(x) )^k
  if(cospoly){
    if( cs2cp(p, a) != 0){
      fprintf(stderr, "Limited intelligence: unable to do a full conversion!\n");
    }
  }

  if(factorout){
    for(max=0.0, j=0; j<2*cs_order; j++){
      x=a[j];
      if((x=fabs(x)) > max) max=x;
    }
  }else{
    max=1.0;
  }

  if(exprform){
    if(title) printf("%s(x)=",title);
    if(cospoly){
      // cosine part
      printf(CPFMT, p[0]);
      for(j=1; j<cs_order; j++)
        printf("+cos(x)*("CPFMT, p[j]);
      for(j=1; j<cs_order; j++) printf(")");

      // sine part
      printf("-sin(x)*("CPFMT, p[cs_order]);
      for(j=1; j<cs_order; j++)
        printf("+cos(x)*("CPFMT, p[j+cs_order]);
      for(j=0; j<cs_order; j++) printf(")");
    }else{
      if(factorout) printf("%g*(  ", max);
      for(j=0; j<2*cs_order; j++){
        x=a[j];
        if(factorout) x/=max;
        if(j==0){
          if(x==1.0) printf("1");
          else printf(CPFMT, x);
        }else{
          printf(CSFMT, x);
          if(j<cs_order) printf("*cos(");
          else printf("*sin(");

          if(j==1 || j==cs_order){
            printf("x)");
          }else if(j<cs_order){
            printf("%d*x)", j);
          }else{
            printf("%d*x)", j-cs_order+1);
          }
        }
      }
      if(factorout) printf("  )");
    }
    printf("\n");
  }

  if(arrayform){
    if(title) printf("%s=",title);
    for(j=0; j<2*cs_order; j++){
      if(j>0) printf(","); else printf("{");
      printf(CPFMT, (cospoly?p[j]:a[j]) );
    }
    printf("}\n");
  }


  return 0;
}

// get cosine/sine transform coefficients of the potential, pmf = -log(hist)
static int pcossin(double a[], double dist[], double weight[], int flag, char *title){
  int i, j, k, dim;
  double sum,ang, x, bj, bk, w;
  static double mat[4*CSORDER*CSORDER];
  int nolog;

  nolog=(flag&CS_NOLOG);

  dim=2*cs_order;
  for(j=0; j<dim; j++)a[j]=0.0;

  if(exclude_last_sine) dim--;

  // fill in the matrix
  for(j=0; j<dim; j++){
    for(k=j; k<dim; k++){
      sum=0.0;
      for(i=0; i<ANGCNT; i++){
        ang=(i*DANG-180.0)*DEG2RAD;
        GET_CS(bj, j);
        GET_CS(bk, k);
        sum+=bj*bk*weight[i];
      }
      mat[j*dim+k]=sum;
      if(k!=j) mat[k*dim+j]=sum;
    }
  }

  // fill in the solution column
  for(j=0; j<dim; j++){
    sum=0.0;
    for(i=0; i<ANGCNT; i++){
      ang=(i*DANG-180.0)*DEG2RAD;
      if(nolog){
        x=dist[i];
      }else{
        if(dist[i] == 0){
          fprintf(stderr, "how to handle zero? at i=%d, ang=%g\n", i, ang*RAD2DEG);
          return 1;
        }
        x=-log(dist[i]);
      }
      w=((weight==NULL)?1.0:weight[i]);

      GET_CS(bj, j);
      sum+=x*bj*w;
    }
    a[j]=sum;
  }

  // solve a by LU decomposition
  lusolve(mat, a, dim);

  // print out coefficients
  if (verbose >= 0)
    cs_print(a, flag&(~CS_COSPOLY), title); // fourier coefficents
  cs_print(a, flag|CS_COSPOLY, title); // cosine polynomial
  return 0;
}

static void psplit(double dist[], double partA[], double partB[], double a0, double a1, double eps){
  int i0, i1, i;
  double k;

  i0=(int)((a0+180)/DANG+1e-6);
  i1=(int)((a1+180)/DANG+1e-6);
  k=(dist[i1]-dist[i0])/(i1-i0);
  for(i=0; i<ANGCNT; i++){
    if(i>=i0 && i<i1){
      partA[i]=dist[i0]+k*(i-i0);
    }else{
      partA[i]=dist[i];
    }
    partB[i]=dist[i]-partA[i];
    if(partA[i]<eps) partA[i]=eps;
    if(partB[i]<eps) partB[i]=eps;
  }
}

// print portable format
static int pgentle(char fnm[], int del){
  FILE *fp;
  int i,j,cnt,iwrap;
  double x,sum1,sum2,max1=0.0,max2=0.0;
  static double arr1[ANGCNT], arr2[ANGCNT];
  if(del<=0) del=1;
  if(ANGCNT%del != 0){
    fprintf(stderr, "%d is not divisible by %d\n", ANGCNT, del);
    return 1;
  }
  printf("grouping %d bins\n", del);
  cnt=ANGCNT/del;
  for(i=0; i<cnt; i++){
    sum1=sum2=0.0;
    //for(j=0; j<=del; j++){
    for(j=-del; j<del; j++){
      iwrap=((i*del+j+ANGCNT*10)%ANGCNT);
      x=Hphi[iwrap];
      //if(j==0 || j==del)x*=0.5;
      sum1+=x;

      x=Hpsi[iwrap];
      //if(j==0 || j==del)x*=0.5;
      sum2+=x;
    }
    //printf("i=%d, %g %g\n", i, sum1,sum2);
    if((arr1[i]=sum1/(2*del))>max1) max1=arr1[i];
    if((arr2[i]=sum2/(2*del))>max2) max2=arr2[i];
  }
  if((fp=fopen(fnm,"w"))==NULL){
    fprintf(stderr, "cannot open %s\n",fnm);
    return 1;
  }
  for(i=0; i<cnt; i++){
    arr1[i]/=max1;
    arr2[i]/=max2;
    fprintf(fp, "%8.3f %g %g\n",
        -180.0+i*del*DANG, arr1[i], arr2[i]);
  }
  fclose(fp);
  printf("spb0_distref=%d %g %g | \\\n", cnt, -M_PI, M_PI);
  for(i=0; i<cnt; i++) {
    printf("%8.6f ", arr1[i]);
    if((i+1)%10==0 && i!=cnt-1) printf("\\\n");
  }
  printf("\n\n");
  printf("spb1_distref=%d %g %g | \\\n", cnt, -M_PI, M_PI);
  for(i=0; i<cnt; i++){
    printf("%8.6f ", arr2[i]);
    if((i+1)%10==0 && i!=cnt-1) printf("\\\n");
  }
  printf("\n\n");
  return 0;
}


// deal with histogram information
static void phistogram(char *fnm_distfs, char *fnm_dist2){
  int i, j, tot=0,ltot, tmp;
  double epsf, epss;
  FILE *fp;
  double phics[CSORDER*2]={0.0},dphi[CSORDER*2]={0.0};
  double psics[CSORDER*2]={0.0},dpsi[CSORDER*2]={0.0};
  static double tmpcs[CSORDER*2], tmpd[ANGCNT];
  int csflags=CS_ARRAY|CS_EXPR;

  for(j=0; j<ANGCNT; j++) Hpsi[j]=0;

  for(i=0; i<ANGCNT; i++){
    ltot=0;
    for(j=0; j<ANGCNT; j++){
      ltot+=(tmp=Hdihs[i][j]);
      Hpsi[j]+=tmp;
    }
    tot+=ltot;
    Hphi[i]=ltot;
  }

  // make averages to make distribution from histogram
  hist2dist(dist_phi, Hphi, average_bins, 1, &epsf);
  hist2dist(dist_psi, Hpsi, average_bins, 1, &epss);

  // construct weight function
  demonize( hist2dist(wei_phi, Hphi, makewei_bins, 0, NULL) );
  demonize( hist2dist(wei_psi, Hpsi, makewei_bins, 0, NULL) );

  // try to obtain FT for phi distributions
  pcossin(phics, dist_phi, wei_phi, csflags, "pf");
/*
  // test conversion
  {
    static double cs1[CSORDER*2], cs2[CSORDER*2];
    // NOTE here we do not allow CS_COSPOLY for funny conversion, so the result is raw data
    cs_print(phics, CS_ARRAY, "raw");
    cs2cp(cs1, phics);
    cs_print(cs1, CS_ARRAY, "cs1");
    cp2cs(cs2, cs1);
    cs_print(cs2, CS_ARRAY, "cs2");
  }
*/
  fixbaseline( cs_combine(dphi, phics, cs_phiref, KT, -1.0, 0.0) );
  cs_print(dphi, csflags|CS_COSPOLY, "dphi");
  fixbaseline( cs_combine(dphi, phics, cs_phiref, KT*2, -1.0, 0.0) );
  cs_print(dphi, csflags|CS_COSPOLY, "dphi1");
  printf("\n");

  // split and calculate coefficients
  psplit(dist_phi, phiA, phiB,-180,-20, epsf);
  demonize( memcpy(wei_phiA, phiA, ANGCNT*sizeof(double)) );
  pcossin(tmpcs, phiA, wei_phiA, csflags, "pfA");
  demonize( memcpy(wei_phiB, phiB, ANGCNT*sizeof(double)) );
  pcossin(tmpcs, phiB, wei_phiB, csflags, "pfB");

  printf("\n==================================================================\n\n");


  // deal with the psi angle
  pcossin(psics, dist_psi, wei_psi, csflags, "ps");
  fixbaseline( cs_combine(dpsi, psics, cs_psiref, KT, -1.0, 0.0) );
  cs_print(dpsi, csflags|CS_COSPOLY, "dpsi");
  fixbaseline( cs_combine(dpsi, psics, cs_psiref, KT*2, -1.0, 0.0) );
  cs_print(dpsi, csflags|CS_COSPOLY, "dpsi1");
  printf("\n");

  psplit(dist_psi, psiA, psiB, -120,40, epss);
  demonize( memcpy(wei_psiA, psiA, ANGCNT*sizeof(double)) );
  pcossin(tmpcs, psiA, wei_psiA, csflags, "psA");
  demonize( memcpy(wei_psiB, psiB, ANGCNT*sizeof(double)) );
  pcossin(tmpcs, psiB, wei_psiB, csflags, "psB");

  if((fp=fopen(fnm_distfs, "w"))==NULL){
    fprintf(stderr, "cannot write dihedral distributions\n");
  }else{
    for(i=0; i<ANGCNT; i++){
      fprintf(fp, "%8.3f %e %e %e %e %e\n",
          i*DANG-180.0, dist_phi[i], dist_psi[i], psiA[i], psiB[i], tmpd[i]);
    }
    fclose(fp);
  }

  if ((fp=fopen(fnm_dist2, "w")) == NULL) {
    fprintf(stderr, "cannot write joint dihedral distributions\n");
  } else {
    for (i = 0; i < ANGCNT; i++) {
      for (j = 0; j < ANGCNT; j++) {
        fprintf(fp, "%8.3f %8.3f %e\n",
            (i+0.5)*DANG-180, (j+0.5)*DANG-180, 1.0*Hdihs[i][j]/tot/(1.0*DANG*DANG));
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

void help_msg(const char *prog)
{
  static char options[]=
    "OPTIONS:\n"
    " -h: print this message\n"
    " -A: followed by analysis level\n"
    " -a: followed by number of bins to make averages\n"
    " -w: followed by number of bins to make a weighting function\n"
    " -o: order of cosine-sine expansion\n"
    " -f: followed by a factor to magnify the weighting function\n"
    " -l: followed by a list file, which contains the files to be handled\n"
    " -e: exclude the highest order sine term\n"
    " -r: followed reference force field 03(default) or 99sb\n"
    " -R: followed by reference structure\n"
    " -S: followed by cartoon style (0 or 1)\n"
    " -g: # of grouping neighoring bins in gentle form\n"
    " -v: verbose\n";
  static char input_desc[]=
    "  I print out dihedral angle information from the input\n"
    "  Use -l file.lst to handle a bunch of files.";

  fprintf(stderr, "%s  Copyright (C) 2010  Cheng Zhang\n"
  "This program comes with ABSOLUTELY NO WARRANTY.  "
  "It is free software, and you are welcome to redistribute "
  "it under certain conditions.\n\n", prog);

  fprintf(stderr, "USAGE:\n");
  fprintf(stderr, "%s [OPTIONS] filename\n\n", prog);
  fprintf(stderr, "%s\n", options);
  fprintf(stderr, "%s\n", input_desc);
  exit(1);
}

/* handle input parameters */
int handle_params(int argc, char *argv[])
{
  int i,itmp,err=0,ch;
  char *follower;

  // handling options and get the PDB file name
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      sscpy(fname_input, argv[i]);
      continue;
    }

    // locate the follower for trailing flags
    ch=argv[i][1];
    if(strchr("AawolfgrRS", ch) != NULL){
      if(argv[i][2]=='\0'){
        if(i==argc-1 || argv[i+1][0]=='-'){
          err=1;
          goto END_OPT;
        }
        follower = argv[++i];
      }else{
        follower=argv[i]+2;
      }
    }

    switch(ch){
    case 'A': // analysis level
      analysis_level=atoi(follower);
      break;
    case 'a': // average how many bins
      average_bins=atoi(follower);
      break;
    case 'w': // average for weights
      makewei_bins=atoi(follower);
      break;
    case 'o': // order
      itmp=atoi(follower);
      if(itmp>CSORDER || itmp<3){
        fprintf(stderr, "the given order %d is unacceptable.\n", itmp);
      }else{
        cs_order=itmp;
      }
      break;
    case 'f': // factor
      demon_factor=strtod(follower, NULL);
      break;

    case 'S': // cartoon style
      cartoon_style=(atoi(follower)!=0);
      break;

    case 'g':
      grouping_neis=atoi(follower);
      break;

    case 'l': // input file list
      sscpy(listname, follower);
      break;

    case 'R':
      sscpy(fname_ref, follower);
      break;

    case 'r': // reference force field
      if( strcmp(follower, "99sb")==0 ||
          strcmp(follower, "99SB")==0 ||
          strcmp(follower, "sb")==0   ||
          strcmp(follower, "SB")==0){
        cp_phiref=cp_phisb;
        cp_psiref=cp_psisb;
      }else if( strcmp(follower, "03")==0){
        cp_phiref=cp_phi03;
        cp_psiref=cp_psi03;
      }
      break;

    case 'v':
      verbose=1;
      if(argv[i][2]=='v') verbose++;
      break;

    case 'e':
      exclude_last_sine=1;
      break;

    case 'h':
    default:
      err=1;
      goto END_OPT;
    }
  }

END_OPT:
  if (err || (listname == NULL && fname_input == NULL))
    help_msg(argv[0]);

  /* promote the analysis level by 2, to allow deeper analysis
   * otherwise for a single file, the level is 0 by default. */
  if (listname != NULL) {
    analysis_level+=2;
    verbose--;
  }else{
    /* be more verbose if the input is just a single file */
    verbose++;
  }

  return 0;
}

int main(int argc, char *argv[])
{
  char *p;
  char *fnm_distfs=NULL;
  char *fnm_dist2=NULL;

  if(0 != handle_params(argc, argv) ) return 1;

  fnm_distfs = ssdup("distfs");
  fnm_dist2 = ssdup("dist2");

  /* load input file(s) */
  if (listname != NULL) { /* load files in the input list */
    FILE *fp;
    char *buf=NULL;
    size_t size;

    if ((fp=fopen(listname, "r")) == NULL) { /* open the file list */
      fprintf(stderr, "cannot open the file list %s\n", listname);
      return 1;
    }
    while (ssfgets(buf, &size, fp)) { /* each line is a file */
      trim(buf);
      if(buf[0] == '\0') break;
      loadxyz(buf, NULL);
      printf("succesfully loading data from %s   \r", buf);
    }
    printf("\n");
    fclose(fp);
    ssdel(buf);
  } else {
    if (loadxyz(fname_input, fname_ref) != 0)
      help_msg(argv[0]);
    printf("succesfully loading data from %s | %s\n",
      fname_input, (fname_ref != NULL) ? fname_ref : "");

    /* by default, we just analyze one structure and quit */
    if (analysis_level < 1)
      goto QUIT;

    /* modify the file names for distributions */
    sscpy(fnm_distfs, fname_input);
    if ((p=strrchr(fnm_distfs, '.')) != NULL) *p='\0';
    sscpy(fnm_dist2, fnm_distfs);
    sscat(fnm_distfs, "_distfs");
    sscat(fnm_dist2,  "_dist2");
  }

  /* convert cosine polynomial (cp) to fourier components (cs)
   * for the force field reference potential */
  cp2cs(cs_phiref, cp_phiref);
  cp2cs(cs_psiref, cp_psiref);

  // generate histograms
  phistogram(fnm_distfs, fnm_dist2);

  printf("\n==================================================================\n\n");

  pgentle("gentle.txt", grouping_neis);
  printf("abins=%d, wbins=%d, order=%d, factor=%g\n",
      average_bins, makewei_bins, cs_order, demon_factor);

QUIT:
  ssdelall(); /* free all strings */
  return 0;
}


/*
  To convert a pdb file to a format that is suitable for simulation using GROMACS+ffamber ports.

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

   The output name is am_original.pdb
   It does the following things:
   * Strip away non-coordinate PDB entries.
   * Remove multiple models/instances (-m to keep multiple models)
   * Optionally remove strange (nonprotein) molecules
   * Remove alternative locations (only keep location A)
   * Convert terminal residues, e.g. first MET ==> NMET
     also rename the C-caps atoms O/OXT to OC1/OC2 respectively
     also rename HN1 and HN2 in NH2 to H1 and H2 respectively
   * Change residue names according to protonation states, e.g. HIS --> HID/HIE/HIP
     prompt a change if the protonation state cannot be determined.
   * Renumbering hydrogen atoms (to avoid -ignh during pdb2gmx)

   It can also reverse the AMBER names in a pdb to normal ones by -R (to de-amberize)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define LMAX 64000
#define LLEN  96
char lines[LMAX][LLEN];
int resids[LMAX];

int remove_strange_residues = 0; // -g to change to 1
int always_lyp = 0; // assume LYS is LYP
int always_hid = 0; // assume HIS is HID
int multiple_models = 0; // allow multiple models
int remove_chainid = 0; // remove chain identifier COL 22
int verbose = 0;
int alter_atname_output = 0; // [1HB ] instead of [ HB1]
int no_ter = 0; // 1 to eliminate TER in PDB
int add_oxt = 0; /* automatically add OXT */

int reverse_to_normal = 0; // reverse an amberized pdb to normal names

const char *fnin = NULL;

#define STDRES_MAX 200
// we don't include ions and water molecules
char stdres[STDRES_MAX][8]={"",
"ALA","ARG","ASN","ASP", "GLU","GLN", "GLY", "ILE","LEU",
"MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
"NALA","NARG","NASN","NASP", "NGLU","NGLN", "NGLY", "NILE","NLEU",
"NMET","NPHE","NPRO","NSER","NTHR","NTRP","NTYR","NVAL",
"CALA","CARG","CASN","CASP", "CGLU","CGLN", "CGLY", "CILE","CLEU",
"CMET","CPHE","CPRO","CSER","CTHR","CTRP","CTYR","CVAL",

"HIS","HID","HIE","HIP","NHID","CHID",

"CYS","CYN","CYM","CYX","CYS2",
"NCYN","NCYM","NCYX",
"CCYN","CCYM","CCYX",

"LYS","LYP","LYN",
"NLYP","NLYN",
"CLYP","CLYN",

"DAB","ORN",  // amino acid stops here
"ACE","NME","NHE","NH2", // caps

"HEM","CMO","OH", // heme
"RCN","RC3","RC5","RC","RGN","RG3","RG5","RG","DCN","DC3","DC5","DC","DGN","DG3","DG5","DG","RUN","RU3","RU5","RU",
"DTN","DT3","DT5","DT","DAN","DA3","DA5","DA","RAN","RA3","RA5","RA",
""};

// names for reverse mapping
char remap_names[128][8]={
"HID","HIS", "HIE","HIS", "HIP","HIS",
"CYN","CYS", "CYM","CYS", "CYX","CYS", "YS2","CYS",
"LYP","LYS", "LYN","LYS",
// the following are for partial recovering, broken names
// this is usually a result of saving coordinates from VMD
"!@#","!@#",
"NAL","ALA", "CAL","ALA",
"NAR","ARG", "CAR","ARG",
"NIL","ILE", "CIL","ILE",
"NLE","LEU", "CLE","LEU",
"NME","MET", "CME","MET",
"NPH","PHE", "CPH","PHE",
"NPR","PRO", "CPR","PRO",
"NSE","SER", "CSE","SER",
"NTH","THR", "CTH","THR",
"NTR","TRP", "CTR","TRP",
"NTY","TYR", "CTY","TYR",
"NVA","VAL", "CVA","VAL",
"NHI","HIS", "CHI","HIS",
"NCY","CYS", "CCY","CYS",
"NLY","LYS", "CLY","LYS",
// NAS,CAS and NGL,CGL cannot be determined
""};


int AALAST; // the index of ORN of the above array
char waterres[32][8]={"","HOH","DOD","T4H","T4P","T3P","TIP","T5P","SOL","WAT",""};

// multiple hydrogen data base for H's that can lead to confusion
// adapted from ffamber03
// modify and run gen_mhdb.c to generate the list
typedef struct tag_multih_t{
  char resname[8];
  char hname[8];
  int count;
}multih_t;
multih_t mhdb[300]={{"","", 0},
{"ACE", "HH3", 3},
{"NME", "HH3", 3},
{"NHE", "H", 2},
{"NH2", "H", 2},
{"ALA", "HB", 3},
{"GLY", "HA", 2},
{"SER", "HB", 2},
{"THR", "HG2", 3},
{"LEU", "HB", 2},
{"LEU", "HD1", 3},
{"LEU", "HD2", 3},
{"ILE", "HG2", 3},
{"ILE", "HG1", 2},
{"ILE", "HD", 3},
{"VAL", "HG1", 3},
{"VAL", "HG2", 3},
{"ASN", "HB", 2},
{"ASN", "HD2", 2},
{"GLN", "HB", 2},
{"GLN", "HG", 2},
{"GLN", "HE2", 2},
{"ARG", "HB", 2},
{"ARG", "HG", 2},
{"ARG", "HD", 2},
{"ARG", "HH1", 2},
{"ARG", "HH2", 2},
{"HID", "HB", 2},
{"HIE", "HB", 2},
{"HIP", "HB", 2},
{"TRP", "HB", 2},
{"PHE", "HB", 2},
{"TYR", "HB", 2},
{"GLU", "HB", 2},
{"GLU", "HG", 2},
{"ASP", "HB", 2},
{"LYP", "HB", 2},
{"LYP", "HG", 2},
{"LYP", "HD", 2},
{"LYP", "HE", 2},
{"LYP", "HZ", 3},
{"ORN", "HB", 2},
{"ORN", "HG", 2},
{"ORN", "HD", 2},
{"ORN", "HE", 3},
{"DAB", "HB", 2},
{"DAB", "HG", 2},
{"DAB", "HD", 3},
{"LYN", "HB", 2},
{"LYN", "HG", 2},
{"LYN", "HD", 2},
{"LYN", "HE", 2},
{"LYN", "HZ", 2},
{"PRO", "HD", 2},
{"PRO", "HG", 2},
{"PRO", "HB", 2},
{"HYP", "HD2", 2},
{"HYP", "HB", 2},
{"CYN", "HB", 2},
{"CYM", "HB", 2},
{"CYX", "HB", 2},
{"CYS2", "HB", 2},
{"MET", "HB", 2},
{"MET", "HG", 2},
{"MET", "HE", 3},
{"ASH", "HB", 2},
{"GLH", "HB", 2},
{"GLH", "HG", 2},
{"CALA", "HB", 3},
{"CGLY", "HA", 2},
{"CSER", "HB", 2},
{"CTHR", "HG2", 3},
{"CLEU", "HB", 2},
{"CLEU", "HD1", 3},
{"CLEU", "HD2", 3},
{"CILE", "HG2", 3},
{"CILE", "HG1", 2},
{"CILE", "HD", 3},
{"CVAL", "HG1", 3},
{"CVAL", "HG2", 3},
{"CASN", "HB", 2},
{"CASN", "HD2", 2},
{"CGLN", "HB", 2},
{"CGLN", "HG", 2},
{"CGLN", "HE2", 2},
{"CARG", "HB", 2},
{"CARG", "HG", 2},
{"CARG", "HD", 2},
{"CARG", "HH1", 2},
{"CARG", "HH2", 2},
{"CHID", "HB", 2},
{"CHIE", "HB", 2},
{"CHIP", "HB", 2},
{"CTRP", "HB", 2},
{"CPHE", "HB", 2},
{"CTYR", "HB", 2},
{"CGLU", "HB", 2},
{"CGLU", "HG", 2},
{"CASP", "HB", 2},
{"CLYP", "HB", 2},
{"CLYP", "HG", 2},
{"CLYP", "HD", 2},
{"CLYP", "HE", 2},
{"CLYP", "HZ", 3},
{"CPRO", "HD", 2},
{"CPRO", "HG", 2},
{"CPRO", "HB", 2},
{"CCYN", "HB", 2},
{"CCYX", "HB", 2},
{"CMET", "HB", 2},
{"CMET", "HG", 2},
{"CMET", "HE", 3},
{"NALA", "H", 3},
{"NALA", "HB", 3},
{"NGLY", "H", 3},
{"NGLY", "HA", 2},
{"NSER", "H", 3},
{"NSER", "HB", 2},
{"NTHR", "H", 3},
{"NTHR", "HG2", 3},
{"NLEU", "H", 3},
{"NLEU", "HB", 2},
{"NLEU", "HD1", 3},
{"NLEU", "HD2", 3},
{"NILE", "H", 3},
{"NILE", "HG2", 3},
{"NILE", "HG1", 2},
{"NILE", "HD", 3},
{"NVAL", "H", 3},
{"NVAL", "HG1", 3},
{"NVAL", "HG2", 3},
{"NASN", "H", 3},
{"NASN", "HB", 2},
{"NASN", "HD2", 2},
{"NGLN", "H", 3},
{"NGLN", "HB", 2},
{"NGLN", "HG", 2},
{"NGLN", "HE2", 2},
{"NARG", "H", 3},
{"NARG", "HB", 2},
{"NARG", "HG", 2},
{"NARG", "HD", 2},
{"NARG", "HH1", 2},
{"NARG", "HH2", 2},
{"NHID", "H", 3},
{"NHID", "HB", 2},
{"NHIE", "H", 3},
{"NHIE", "HB", 2},
{"NHIP", "H", 3},
{"NHIP", "HB", 2},
{"NTRP", "H", 3},
{"NTRP", "HB", 2},
{"NPHE", "H", 3},
{"NPHE", "HB", 2},
{"NTYR", "H", 3},
{"NTYR", "HB", 2},
{"NGLU", "H", 3},
{"NGLU", "HB", 2},
{"NGLU", "HG", 2},
{"NASP", "H", 3},
{"NASP", "HB", 2},
{"NLYP", "H", 3},
{"NLYP", "HB", 2},
{"NLYP", "HG", 2},
{"NLYP", "HD", 2},
{"NLYP", "HE", 2},
{"NLYP", "HZ", 3},
{"NPRO", "H", 2},
{"NPRO", "HD", 2},
{"NPRO", "HG", 2},
{"NPRO", "HB", 2},
{"NCYN", "H", 3},
{"NCYN", "HB", 2},
{"NCYX", "H", 3},
{"NCYX", "HB", 2},
{"NMET", "H", 3},
{"NMET", "HB", 2},
{"NMET", "HG", 2},
{"NMET", "HE", 3},
{"DA5", "H5'", 2},
{"DA5", "H6", 2},
{"DA5", "H2'", 2},
{"DA", "H5'", 2},
{"DA", "H6", 2},
{"DA", "H2'", 2},
{"DA3", "H5'", 2},
{"DA3", "H6", 2},
{"DA3", "H2'", 2},
{"DAN", "H5'", 2},
{"DAN", "H6", 2},
{"DAN", "H2'", 2},
{"RA5", "H5'", 2},
{"RA5", "H6", 2},
{"RA", "H5'", 2},
{"RA", "H6", 2},
{"RA3", "H5'", 2},
{"RA3", "H6", 2},
{"RAN", "H5'", 2},
{"RAN", "H6", 2},
{"DT5", "H5'", 2},
{"DT5", "H7", 3},
{"DT5", "H2'", 2},
{"DT", "H5'", 2},
{"DT", "H7", 3},
{"DT", "H2'", 2},
{"DT3", "H5'", 2},
{"DT3", "H7", 3},
{"DT3", "H2'", 2},
{"DTN", "H5'", 2},
{"DTN", "H7", 3},
{"DTN", "H2'", 2},
{"RU5", "H5'", 2},
{"RU", "H5'", 2},
{"RU3", "H5'", 2},
{"RUN", "H5'", 2},
{"DG5", "H5'", 2},
{"DG5", "H2", 2},
{"DG5", "H2'", 2},
{"DG", "H5'", 2},
{"DG", "H2", 2},
{"DG", "H2'", 2},
{"DG3", "H5'", 2},
{"DG3", "H2", 2},
{"DG3", "H2'", 2},
{"DGN", "H5'", 2},
{"DGN", "H2", 2},
{"DGN", "H2'", 2},
{"RG5", "H5'", 2},
{"RG5", "H2", 2},
{"RG", "H5'", 2},
{"RG", "H2", 2},
{"RG3", "H5'", 2},
{"RG3", "H2", 2},
{"RGN", "H5'", 2},
{"RGN", "H2", 2},
{"DC5", "H5'", 2},
{"DC5", "H4", 2},
{"DC5", "H2'", 2},
{"DC", "H5'", 2},
{"DC", "H4", 2},
{"DC", "H2'", 2},
{"DC3", "H5'", 2},
{"DC3", "H4", 2},
{"DC3", "H2'", 2},
{"DCN", "H5'", 2},
{"DCN", "H4", 2},
{"DCN", "H2'", 2},
{"RC5", "H5'", 2},
{"RC5", "H4", 2},
{"RC", "H5'", 2},
{"RC", "H4", 2},
{"RC3", "H5'", 2},
{"RC3", "H4", 2},
{"RCN", "H5'", 2},
{"RCN", "H4", 2},
{"HEME", "HAA", 2},
{"HEME", "HBA", 2},
{"HEME", "HMA", 3},
{"HEME", "HMB", 3},
{"HEME", "HBB", 2},
{"HEME", "HMC", 3},
{"HEME", "HBC", 2},
{"HEME", "HMD", 3},
{"HEME", "HAD", 2},
{"HEME", "HBD", 2},
{"","", 0}};


// atom name should always have 4 letters, the next character is the alternative location
#define GET_ATNAME(name, line) alter_atname( trim( memcpy(name, line+12, 4) ), 0)
// allow residue name to have 4 letters
#define GET_RESNAME(name, line) trim( memcpy(name, line+17, 4) )


static int mh_lookup(char *resname, char *atname)
{
  int i;
  for (i = 1; mhdb[i].resname[0] != '\0'; i++) {
    if (strcmp(resname, mhdb[i].resname) == 0 &&
       strcmp(atname, mhdb[i].hname) == 0) {
      return mhdb[i].count;
    }
  }
  return -1;
}

// determine if it is a common protein/ligand residue
static int is_stdres(char *name)
{
  int i;
  for (i = 1; stdres[i][0]; i++) {
    if (strcmp(name, stdres[i]) == 0) return i;
  }
  return 0;
}

// determine if it is water residue
static int is_waterres(char *name)
{
  int i;
  for (i = 1; waterres[i][0]; i++) {
    if (strcmp(name, waterres[i]) == 0) return i;
  }
  return 0;
}

// remove leading/tailing spaces
static char *trim(char *s)
{
  char *p = s+strlen(s)-1, *q;
  // trim the tail
  while( isspace(*p) && p >= s)  *p--='\0';

  // trim the leading spaces
  // first nonspace place
  for (q = s; *q && isspace(*q); q++) ;
  // shifting
  for (p = s; *q; ) *p++=*q++;
  if (*q=='\0') *p='\0';
  return s;
}

/* add a prefix 'C' or 'N' to the residue name */
static char *addcn(char *name, char cn)
{
  char name2[8];

  name2[0] = cn; name2[1]='\0';
  strcat(name2, name);
  strcpy(name, name2);
  return name;
}

// handle alternative form of atom names, like 1HD1, 2HD1, 3HD1, instead of HD11, HD12 and HD13, etc
// flag=0 to normal form, flag=1 to alternative form
static char *alter_atname(char *atname, int flag)
{
  char atname2[8]=" ";
  if (flag == 0) { // normal form
    if (isdigit(atname[0])) {
      sprintf(atname2, "%s%c", atname+1, atname[0]);
      strcpy(atname, atname2); // copy back
    }
  } else { // to alternative form
    if (!isdigit(atname[0])) {
      int size = strlen(atname);
      if (!isdigit(atname[size-1])) return atname;
      atname2[0] = atname[size-1];
      atname[size-1]='\0';
      strcpy(atname2+1, atname);
      strcpy(atname, atname2); // copy back
      //printf("alternative form is [%s]\n", atname);
    }
  }
  return atname;
}

static char *patch_atname(char *atname, char *line)
{
  char atname2[8]=" ";
  int i;

  if (strlen(atname) > 4) {
    fprintf(stderr, "Atom name is too long! %s\n", atname);
    return line;
  }

  if (alter_atname_output) {
    strcpy(atname2, atname);
    alter_atname(atname2, 1);
  } else {
    // add a leading space for short atom name
    if (strlen(atname) <= 3) strcat(atname2, atname);
    else strcpy(atname2, atname);
  }

  // append spaces
  for (i = strlen(atname2); i < 4; i++) atname2[i]=' ';
  atname2[4]='\0';

  memcpy(line+12, atname2, 4);
  return line;
}

static int get_inputch(void)
{
  char inputbuf[32]="";
  int ch;
  fgets(inputbuf, sizeof inputbuf, stdin);
  ch = inputbuf[0];
  if (verbose) {
    if (isspace(ch)) {
      printf("You gave me a space.\n");
    } else {
      printf("Your input charactor is %c\n", ch);
    }
  }
  return ch;
}


// this is a mini independent program to reverse an amber file to normal convention
static int deamberize(void)
{
  FILE *fp;
  char fname[FILENAME_MAX];
  char buf[LLEN];
  char atname[8]="", atname2[8]="", resname[8]="", resname2[8]="";
  int i, j, nlin, cap, mod, resid = 0, broken = 0;

  if ((fp = fopen(fnin, "r")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fnin);
    return 1;
  }
  for (i = 0; i < LMAX; ) {
    if (fgets(buf, sizeof buf, fp) == NULL) break;

    if (strncmp(buf, "ATOM", 4) != 0 &&
       strncmp(buf, "HETATM", 6) != 0 &&
       strncmp(buf, "TER", 3) != 0) {
      // just copy the line;
      resids[i] = -1;
    } else {
      sscanf(buf+23, "%d", &resid);
      resids[i] = resid;
    }
    strcpy(lines[i++], buf);
  }
  nlin = i;
  fclose(fp);

  for (i = 0; i < nlin; ) {
    // don't do any modification
    if (resids[i] < 0) { i++; continue; }

    GET_ATNAME(atname, lines[i]);
    GET_RESNAME(resname, lines[i]);
    strcpy(resname2, resname);

    mod = 0;
    cap = 0;
    // decap
    if (strlen(resname) == 4) {
      if (strcmp(resname, "CYS2") == 0) {
        strcpy(resname2, "CYS");
        mod = 1;
      } else if (strchr("CN", (cap = resname[0])) != NULL) {
        strcpy(resname2, resname+1);
        mod = 1;
      }
    } else if (strlen(resname) == 3 && strchr("CN", resname[0]) != NULL) {
      int who = 0;
      // investigate deeper
      if ( strcmp(resname+1, "AS") == 0 ) { // ASP or ASN
        strcpy(resname2, "ASP");
        for (j = i; resids[j] == resids[i]; j++) {
          GET_ATNAME(atname2, lines[j]);
          //printf("j=%d, atname=%s\n", j, atname2);
          if (strcmp(atname2, "OD2") == 0) {
            who = 0; break;
          } else if (strcmp(atname2, "ND2") == 0) {
            who = 1; break;
          }
        }
        if (who == 1) strcpy(resname2, "ASN");
        fprintf(stderr, "assume residue %d %s as %s\n",
          resids[i], resname, resname2);
        cap = resname[0];
        mod = 1;
      } else if (strcmp(resname+1, "GL") == 0 ) { // GLY or GLN, GLU
        strcpy(resname2, "GLY");
        for (j = i; resids[j] == resids[i]; j++) {
          GET_ATNAME(atname2, lines[j]);
          if (strcmp(atname2, "OE2") == 0) {
            who = 1; break;
          } else if (strcmp(atname2, "NE2") == 0) {
            who = 2; break;
          }
        }
        if (who == 1) strcpy(resname2, "GLU");
        else if (who == 2) strcpy(resname2, "GLN");
        fprintf(stderr, "will assume residue %d %s as %s\n",
          resids[i], resname, resname2);
        cap = resname[0];
        mod = 1;
      }// other fall through normal mapping substitution
    }

    // name conversions for LYS, HIS, and CYS
    for (broken = 0, j = 0; remap_names[2*j][0]!='\0'; j++) {
      if (strcmp(remap_names[2*j], "!@#") == 0) {
        broken = 1; // start with broken names
        continue;
      }
      //printf("%2d:%s,%s.\n", j, remap_names[2*j],resname2);
      if (strcmp(remap_names[2*j], resname2) == 0) {
        strcpy(resname2, remap_names[2*j+1]);
        //printf("renaming %s to %s\n", resname, resname2);
        mod = 1;
        if (broken) cap = remap_names[2*j][0];
        break;
      }
    }


    if (mod) {
      printf("renaming %s to %s\n", resname, resname2);

      // append spaces to 4
      for (j = strlen(resname2); j < 4; j++) resname2[j]=' ';
      resname2[j] = '\0';
      for (j = i; resids[j] == resids[i]; j++) {
        // make caps modification
        if (cap == 'C') {
          GET_ATNAME(atname2, lines[j]);
          if (strcmp(atname2, "OC1") == 0) {
            printf("fixing C-terminal atom: OC1-->O\n");
            patch_atname("O", lines[j]);
          } else if (strcmp(atname2, "OC2") == 0) {
            printf("fixing C-terminal atom: OC2-->OXT\n");
            patch_atname("OXT", lines[j]);
          }
        }
        memcpy(lines[j]+17, resname2, 4);
      }
      i = j; // goto the next residue
    } else {
      i++;
    }
  }
  strcat(strcpy(fname, fnin),".bak");
  rename(fnin, fname);

  // output the file
  if ((fp = fopen(fnin, "w")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fnin);
    return 1;
  }
  for (i = 0; i < nlin; i++) {
    fprintf(fp, "%s", lines[i]);
  }
  fclose(fp);
  return 0;
}

const char *prog = "amberize";

/* print help message and die */
static void help(void)
{
  const char options[]=
  "OPTIONS:\n"
  " -h: print this message\n"
  " -g: remove all strange (non-protein) residues (like water)\n"
  "     be careful if you want to keep ligands\n"
  " -l: when uncertain, assume every lysine has 3 HZs, or LYP\n"
  " -d: when uncertain, assume every histodine only has delta H or HID(risky)\n"
  " -m: allow multiple models (for NMR structures)\n"
  " -c: remove chain identifier\n"
  " -q: alternative form outputf for hydrogen, i.e., use [1HD1], [2HD1] instead of [HD11], [HD12]\n"
  " -t: remove TER commands in PDB\n"
  " -v: verbose, -vv to be more\n"
  " -u: (also -D or -R) reverse a PDB to normal (nonamber) names\n"
  " -x: add OXT automatically\n"
  "\n";
  fprintf(stderr, "%s  Copyright (C) 2010  Cheng Zhang\n"
  "This program comes with ABSOLUTELY NO WARRANTY. "
  "It is free software, and you are welcome to redistribute "
  "it under certain conditions.\n\n", prog);

  fprintf(stderr, "USAGE:\n");
  fprintf(stderr, "%s [OPTIONS] your.pdb\n\n", prog);
  fprintf(stderr, "%s", options);
  exit(1);
}

/* handling options and get the PDB file name */
static int doargs(int argc, char **argv)
{
  int i, j, ch;

  prog = argv[0];
  verbose = 0;
  fnin = NULL;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnin = argv[i];
      continue;
    }
    for (j = 1; (ch = argv[i][j]) != '\0';  j++)
      switch (ch) {
      case 'D': case 'R': case 'u':
        reverse_to_normal = 1; break;
      case 'g': remove_strange_residues = 1; break;
      case 'l': always_lyp = 1; break;
      case 'd': always_hid = 1; break;
      case 'm': multiple_models = 1; break;
      case 'c': remove_chainid = 1; break;
      case 'q': alter_atname_output = 1; break;
      case 't': no_ter = 1; break;
      case 'x': add_oxt = 1; break;
      case 'v': verbose++; break;
      case 'h': default: help();
      }
  }
  if (fnin == NULL) help();
  return 0;
}

void diffnorm(double r[], double o[])
{
  int i;
  double x = 0.;

  for (i = 0; i < 3; i++) {
    r[i] -= o[i];
    x += r[i]*r[i];
  }
  x = 1.0/sqrt(x);
  for (i = 0; i < 3; i++) r[i] *= x;
}

/* add missing OXT  */
int addoxt(int start, int n)
{
  int i, end = -1, flags = 0, id;
  double ca[3], c[3], o[3], oxt[3];
  const char *p;
  char buf[80];

  if (n >= LMAX) {
    fprintf(stderr, "too many lines, cannot add more\n");
    exit(1);
  }

  for (i = start; i < n; i++)
    if (strncmp(lines[i], "ATOM", 4) != 0) break;
  end = i;

  /* search for coordinate for CA, C, O */
  #define GETXYZ(r, s) sscanf(s+30, "%lf%lf%lf", r, r+1, r+2)
  for (i = start; i < end; i++) {
    p = lines[i] + 13;
    if (strncmp(p, "CA", 2) == 0) {
      GETXYZ(ca, lines[i]);
      flags |= 1;
    } else if (strncmp(p, "C ", 2) == 0) {
      GETXYZ(c, lines[i]);
      flags |= 2;
    } else if (strncmp(p, "O ", 2) == 0) {
      GETXYZ(o, lines[i]);
      flags |= 4;
    }
  }
  if (flags != 7) {
    fprintf(stderr, "cannot find all atoms %d\n", flags);
    exit(1);
  }

  sscanf(lines[end-1] + 6, "%d", &id); id++;
  strcpy(buf, lines[end-1] + 17); buf[14] = '\0';

  diffnorm(ca, c);
  diffnorm(o, c);
  for (i = 0; i < 3; i++) {
    oxt[i] = c[i] - (ca[i] + o[i])*1.23;
  }

  /* nudge lines after */
  for (i = n; i > end; i--)
   strcpy(lines[i], lines[i-1]);

  sprintf(lines[end], "ATOM  %5d  OXT %s%7.3f %7.3f %7.3f  1.00  0.00           O  \n",
      id, buf, oxt[0], oxt[1], oxt[2]);
  //printf("%s", lines[end]); getchar();
  return n + 1;
}

int main(int argc, char **argv)
{
  FILE *fp;
  char fname[FILENAME_MAX];
  char buf[LLEN], *p, *line;
  char atname[8]="", atname2[8]="", resname[8]="", resname2[8]="";
  int i, j, nlin;
  int resid, atomid, resntype;
  int resfirst = 1000000, reslast = -1;
  int hasncap = 0, hasccap = 0;

  if (doargs(argc, argv) != 0)
    return 1;

  if (reverse_to_normal) {
    return deamberize();
  }

  if ((fp = fopen(fnin, "r")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fnin);
    return 1;
  }
  atname[4] = '\0';
  resname[4] = '\0';
  AALAST = is_stdres("ORN");

  // First round, read pdb line by line,
  // do simple modifications.
  fprintf(stderr, "First round...\n");
  for (i = 0; i < LMAX; ) {
    if (fgets(buf, sizeof buf, fp) == NULL) break;

    // we only need one model
    if (multiple_models) {
      if (strncmp(buf, "MODEL", 5) == 0 ||
         strncmp(buf, "ENDMDL", 6) == 0) {
        resids[i] = -1;
        strcpy(lines[i++], buf);
        //printf("line: %4d, buf=[%s]\n", i, trim(buf));
        continue;
      }
    } else {
      if (strncmp(buf, "ENDMDL", 6) == 0) break;
    }
    if (strncmp(buf, "ATOM", 4) != 0 &&
       strncmp(buf, "HETATM", 6) != 0) {
      if (no_ter) continue;
      if (strncmp(buf, "TER", 3) != 0) continue;
    }

    // change D to H (12-15) is atom name
    for (p = buf+12; *p && (isspace(*p) || isdigit(*p)); p++) ;
    if (*p == 'D') *p = 'H';

    // get atom id, atom name, residue id, residue name
    sscanf(buf+6, "%d", &atomid); // 6-10
    GET_ATNAME(atname, buf);
    GET_RESNAME(resname, buf);
    resntype = is_stdres(resname);
    sscanf(buf+23, "%d", &resid);
    resids[i] = resid; // save the residue id

    if (remove_chainid && buf[21]!=' ') {
      if (verbose) {
        fprintf(stderr, "remove chain identifier [%c]\n", buf[21]);
      }
      buf[21]=' ';
    }

    // remove alternative locations
    if (buf[16] != ' ') {
      if (buf[16] != 'A') {
        if (verbose)
          fprintf(stderr, "discard alternative location %c for atom %d, residue %d\n",
            buf[16], atomid, resid);
        continue;
      } else {
        buf[16]=' ';
      }
    }

    // remove water etc...
    if (remove_strange_residues &&
      resntype == 0) {
      // don't prompt if it's just water
      if (!is_waterres(resname) || verbose) {
        fprintf(stderr, "removing strange atom %s,%d residue %s,%d\n",
          atname, atomid, resname, resid);
      }
      continue;
    }

    if (resntype > 0 && resntype <= AALAST) { // amino acid residue
      if (resid < resfirst) resfirst = resid;
      if (resid > reslast)  reslast =resid;
    }
    if (strcmp(resname, "ACE") == 0) hasncap = 1;
    if (strcmp(resname, "NME") == 0 ||
       strcmp(resname, "NHE") == 0 ||
       strcmp(resname, "NH2") == 0) hasccap = 1;

    if (strcmp(resname, "NH2") == 0) {
      if (strcmp(atname, "HN1") == 0) {
        patch_atname("H1", buf);
      } else if (strcmp(atname, "HN2") == 0) {
        patch_atname("H2", buf);
      }
    }

    strcpy(lines[i], buf);
    i++;
    if (i >= LMAX) {
      fprintf(stderr, "the file is too large to handle\n");
      break;
    }
  }
  nlin = i;
  fclose(fp);

  fprintf(stderr, "Lines:%d, protein residues: %d to %d, ncap=%d, ccap=%d\n",
      nlin, resfirst, reslast, hasncap, hasccap);
  if (hasccap) reslast = 100000000;
  if (hasncap) resfirst = -100000000;

  // Second round: advanced modifications (residue names, cap atoms)
  fprintf(stderr, "Second round: residue-renaming...\n");
  for (i = 0; i < nlin; ) {
    int hash, has1, has2, resmod = 0, ch;
    if (resids[i] < 0) { i++; continue; }
    line = lines[i];
    GET_RESNAME(resname, line);
    resntype = is_stdres(resname);
    //printf("resname=%s, line=%4d, type=%d\n", resname, i, resntype);

#define ADDCAP(resnm,id) \
    if (id == resfirst) { resmod = 2; addcn(resnm, 'N'); } \
    else if (id == reslast) { resmod = 3; addcn(resnm, 'C'); }

    /* modify LYS to LYP (+) or LYN (neutral) */
    if (strcmp(resname, "LYS") == 0) {
      has1 = hash = 0;
      for (j = i; j < nlin; j++) {
        if (resids[j] != resids[i]) break;
        GET_ATNAME(atname, lines[j]);
        if (atname[0] == 'H') hash = 1;
        if (strcmp(atname, "HZ3") == 0) has1 = 1;
      }
      if (!hash) { // X-ray structure, prompt change
        strcpy(resname2,"LYP");
        // only LYP can be at the caps
        if (resids[i] != resfirst && resids[i] != reslast && !always_lyp) {
          fprintf(stderr, "I am about to change residue %d LYS.  ", resids[i]);
          fprintf(stderr, "If it only has two HZs, enter n, "
            "otherwise just press anything else (default): ");
          ch = get_inputch();
          if (ch=='n'||ch=='N') {
            strcpy(resname2,"LYN");
          }
        }
        ADDCAP(resname2, resids[i]);
      } else if (has1) {
        strcpy(resname2, "LYP");
        ADDCAP(resname2, resids[i]);
      } else {
        strcpy(resname2, "LYN");
        ADDCAP(resname2, resids[i]);
      }
      //fprintf(stderr, "modify the lysine residue\n");
      if (resmod == 0) resmod = 1;

    // modify CYS
    } else if (strcmp(resname, "CYS") == 0) {
      has1 = hash = 0;
      for (j = i; j < nlin; j++) {
        if (resids[j] != resids[i]) break;
        GET_ATNAME(atname, lines[j]);
        if (atname[0] == 'H') hash = 1;
        if (strcmp(atname, "HG") == 0) has1 = 1;
      }
      if (!hash || has1) {
        strcpy(resname2,"CYN");
        if (!hash) { // X-ray structure, be careful
          fprintf(stderr, "I am about to change residue %d CYS.  ", resids[i]);
          fprintf(stderr, "If it is involved in a disulfide bond, enter x or y "
              "(otherwise just press anything else): ");
          ch = get_inputch();
          if (ch == 'x' || ch == 'X' || ch == 'y' || ch == 'Y') {
            strcpy(resname2, "CYX");
            fprintf(stderr, "Now make sure you use the modified specbond.dat.\n");
          }
        }
        ADDCAP(resname2, resids[i]);
      } else {
        strcpy(resname2,"CYM");
        fprintf(stderr, "Warning: you might want to use CYX (charge neutral) instead of CYM\n");
        if (resids[i] == resfirst || resids[i] == reslast) {
          fprintf(stderr, "Active cysteine terminal!\n"
              "I'll change it to CYX, but then you need to modify specbond.dat if you want "
              "pdb2gmx to create disulfide bond for you.  "
              "But you can always manually add disulfide bond.\n");
          strcpy(resname2,"CYX");
          ADDCAP(resname2, resids[i]);
        }
      }
      if (resmod == 0) resmod = 1;

    // modify HIS
    } else if (strcmp(resname, "HIS") == 0) {
      has1 = has2 = hash = 0;
      for (j = i; j < nlin; j++) {
        if (resids[j] != resids[i]) break;
        GET_ATNAME(atname, lines[j]);
        if (atname[0]=='H') hash = 1;
        if (strcmp(atname, "HD1") == 0) has1 = 1;
        if (strcmp(atname, "HE2") == 0) has2 = 1;
      }
      if (!hash) { // X-ray
        strcpy(resname2,"HID");
        if (!always_hid) {
          fprintf(stderr, "I am about to change residue %d HIS.  ", resids[i]);
          fprintf(stderr, "Press d (has HD1, default), e (has HE2) or p (has both) "
            "to specify the protonation state:");
          ch = get_inputch();
          if (ch=='e'||ch=='E') {
            strcpy(resname2,"HIE");
          } else if (ch=='p'||ch=='P') {
            strcpy(resname2,"HIP");
          }
        }
      } else {
        if (has1 && has2) {
          strcpy(resname2, "HIP");
        } else if (has2&&!has1) {
          strcpy(resname2, "HIE");
        } else { // default
          strcpy(resname2, "HID");
        }
      }
      // terminal caps
      if (resids[i] == resfirst || resids[i] == reslast) {
        if (resname2[2] != 'D') {
          fprintf(stderr, "we do not have %s for terminals, revert to HID\n", resname2);
          strcpy(resname2, "HID");
        }
        ADDCAP(resname2, resids[i]);
      }
      if (resmod == 0) resmod = 1;
    // regular residues
    } else if (strlen(resname) == 3) { // not modified names
      strcpy(resname2, resname);
      ADDCAP(resname2, resids[i]);
    }

    // residue name is modified,
    if (resmod) {
      fprintf(stderr, "[%d] change residue %3d: from [%s] to [%s]\n",
        resmod, resids[i], resname, resname2);
      // append spaces after residue name, to facilitate copying
      for (j = strlen(resname2); j < 4; j++) resname2[j]=' ';
      resname2[j]='\0';

      // check C-terminal OXT
      if (resmod == 3) {
        for (has1 = 0, j = i; j < nlin; j++) {
          if (resids[j] != resids[i]) break;
          GET_ATNAME(atname, lines[j]);
          if (strcmp(atname, "OXT") == 0) {
            has1 = 1;
            break;
          }
        }
        if (!has1) {
          if (add_oxt) { /* automatically add OXT coordinates */
            fprintf(stderr, "automatically add OXT terminal...\n");
            nlin = addoxt(i, nlin);
            resids[nlin-1] = resids[nlin-2];
            has1 = 1;
          } else {
            fprintf(stderr, "cannot find OXT, C-terminal is incomplete!\n");
            return 1;
          }
        }
        for (j = i; j < nlin; j++) {
          if (resids[j] != resids[i]) break;
          GET_ATNAME(atname, lines[j]);
          if (strcmp(atname, "O") == 0) {
            memcpy(lines[j]+13, "OC1", 3);
            //printf("line:%s", lines[j]);
          } else if (strcmp(atname, "OXT") == 0) {
            memcpy(lines[j]+13, "OC2", 3);
            //printf("line:%s", lines[j]);
          }
        }
      }

      // copy the modified resiude name, resname2
      for (j = i; j < nlin; j++) {
        if (resids[j] != resids[i]) break;
        memcpy(lines[j]+17, resname2, 4);
      }
      i = j; // goto the next residue
    } else {// end of modifing the residue name
      i++; // next line
    }
  }

  // Third round: re-numbering hydrogen atoms (assume they are continuous)
  fprintf(stderr, "Third round: H-renumbering...\n");
  for (i = 0; i < nlin; i++) {
    int mhcnt, size, size2, mhid, id0;

    if (resids[i] < 0) continue;
    line = lines[i];
    GET_ATNAME(atname, line);
    GET_RESNAME(resname, line);

    size = strlen(atname);
    if (atname[0]!='H' || !isdigit(atname[size-1]) ) continue;
    id0 = atname[--size];
    atname[size]='\0';
    mhcnt = mh_lookup(resname, atname);
    if (mhcnt <= 0) continue;

    for (j = i; j < i+mhcnt && j < nlin; j++) {
      GET_ATNAME(atname2, lines[j]);
      GET_RESNAME(resname2, lines[j]);

      size2 = strlen(atname2)-1;
      if (size2 != size || !isdigit(atname2[size2]) ||
          (strncmp(atname, atname2, size2) != 0) || (strcmp(resname, resname2) != 0) ) {
        fprintf(stderr, "a possible PDB corruption, atom:%s(%c), res:%s; i=%d,j=%d,mhcnt=%d, size=%d,%d;%c;%s,%s;%s,%s\nline:%s\n",
            atname, id0, resname, i, j, mhcnt, size, size2, atname2[size2], atname, atname2, resname, resname2, lines[j]);
        break;
      }
      mhid = j-i+1;

      // we do it even if atname2[size2]-'0' == mhid, because user may want to use the alternative form
      sprintf(atname2, "%s%d", atname, mhid);
      if (verbose) {
        if (mhid == 1 || verbose > 1) fprintf(stderr, "patching hydrogen %s, residue %s(%d)\n", atname2, resname, resids[i]);
        if (verbose > 2) fprintf(stderr, "  Before: %s", lines[j]);
      }
      patch_atname(atname2, lines[j]);
      if (verbose > 2) {
        fprintf(stderr, "  After : %s", lines[j]);
      }
    }
    i += (mhcnt-1);
    //printf("find a confusing H, at:%s, res:%s.\n", atname, resname);
  }

  sprintf(fname, "am_%s", fnin);
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fname);
    return 1;
  }
  for (i = 0; i < nlin; i++) {
    fprintf(fp, "%s", lines[i]);
  }
  fclose(fp);
  return 0;
}


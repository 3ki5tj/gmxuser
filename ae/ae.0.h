#define HAVEREAL 1
#define ZCOM_PICK
#define ZCOM_LOG
#define ZCOM_CFG
#define ZCOM_RNG
#include "zcom.h"
#define USE_MPI GMX_MPI

/* BOLTZ is defined in GROMACS */
#ifndef BOLTZ
#define BOLTZ  8.314511212e-3
#endif

typedef struct {
  int mode;         /* command-line input
                       $usr:cfg */
  double tmstep;    /* MD integration step, for convenience;
                       $def: 0.002; $usr:cfg; $io:; */
  double T0;        /* thermostat temperature; $key: T0; $def: 300.0; */
  double beta;      /* current inverse temperature; $io:none; */
  real   scale;     /* force scaling factor; $def: 1.0f; $io:none; */
  char   *rngfile;  /* file name of random number state;
                       $def: NULL; $closecall: if ($ismaster) mtsave(@@); */
  int    nstlog;    /* interval of writing trace file; -1: only at ns, 0: disable; $def: -1;  */
  char   *logfile;  /* name of trace file; $def: "mdae.log";
                       $closecall: if ($ismaster) log_close(@log); */
  logfile_t *log;   /* logfile $cfg:none; */
} ae_t; /* $cfgopen; $rb:0; $wb:0; $reduce:0; $bcast:0; */

static void ae_updscl(ae_t *ae)
{
  ae->scale = (real) (ae->beta / (BOLTZ * ae->T0));
}

static int ae_loaddata(ae_t *ae, int iscont)
{
  if (!iscont) return 0; /* initial run */
  /* read in data */
}

/* dump the object data, for an array, only the first and last
 * `arrmax' elements are printed */
static void ae_dump(ae_t *ae, const char *fn, int arrmax)
{
  FILE *fp = stderr;

  if ((fp = fopen(fn, "w")) == NULL)
    fprintf(stderr, "cannot write %s\n", fn);
  ae_manifest(ae, fp, arrmax); /* objgen function */
  if (fp && fp != stderr) fclose(fp);
}

/* master only basic initialization */
static ae_t *ae_open(const char *fncfg, int iscont, double tmstep, int mode)
{
  ae_t *ae;

  if (fncfg == NULL) fncfg = "ae.cfg";
  die_if ((ae = ae_cfgopen(fncfg, mode, tmstep)) == NULL,
    "cannot open cfg file %s\n", fncfg);
  ae->log = log_open(ae->logfile);
  ae->beta = 1.0/(BOLTZ * ae->T0);

  die_if (ae_loaddata(ae, iscont) != 0,
      "cannot load previous data, %d", iscont);
  ae_updscl(ae);
  ae_dump(ae, "ae.manifest", 3); /* objgen function */
  return ae;
}

static int ae_move(ae_t *ae, gmx_large_int_t step,
   int bfirst, int blast, int br)
{
  
  /* grab the local temperature of the potential energy */
  ae->beta = alge_getk(ae->epot);
  ae_updscl(ae);
  return 0;
}

/* $objgen_funcs; */
/* $OBJ_FUNCS */

/* $BONDFREE_C */
/* $FORCE_C */
/* $SIMUTIL_C */
/* $MD_C */
/* $RUNNER_C */
/* $MDRUN_C */

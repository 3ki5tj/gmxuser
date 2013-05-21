''' write auxiliary code to handle input '''

import re, getopt, os, sys
from ccgmx import CCGMX

class CCAUX(CCGMX):

  def __init__(c, obj):
    ''' override the constructor '''
    # since CCAUX is a CCGMX, gromacs version is computed
    CCGMX.__init__(c, None, obj, {})


  struct_src = r'''
typedef struct {
  int mode;
} $OBJ_t;
'''

  init_src = r'''
/* initialize $OBJ, every node calls this */
$OBJ_t *$PFX_init(const char *fncfg, unsigned fromcpt,
    gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr, int mode)
{
  $OBJ_t *$OBJ = NULL;
  int err = 0;

  if (mtop == NULL) {
    fprintf(stderr, "no topology\n");
    return NULL;
  }

  if (!(cr->duty & DUTY_PP)) /* return NULL for a PME-only node */
    return NULL;

  if (SIMMASTER(cr)) {
    /* call gromacs independent routine to initalize */
    if (($OBJ = calloc(1, sizeof(*$OBJ))) == NULL) {
      err = -1;
      goto PAR1;
    }
    $OBJ->mode = mode;
  } else { /* for non-master, simply allocate space */
    if ( ($OBJ = calloc(1, sizeof(*$OBJ)) ) == NULL) {
      fprintf(stderr, "no memory for $OBJ_t.\n");
      return NULL;
    }
  }

PAR1: /* check error in the master branch */
  if (PAR(cr)) gmx_bcast(sizeof(err), &err, cr);
  if (err) return NULL;

#ifdef GMX_MPI
  /* tell everyone settings on the master
   * valid only for PP only node, maybe we need to
   * consider using mpi_comm_mysim for more advanced versions
   * we pass MPI_COMM_NULL to avoid the case of one-node-mpi */
#endif

  return $OBJ;
}'''


  finish_src = '''
void $PFX_finish($OBJ_t *$OBJ)
{
  free($OBJ);
}
'''

  move_src = r'''
/* only a PP-processor calls this function
 * also assume global_stat() has been called
 * SIMMASTER(cr) should be equivalent to MASTER(cr) */
int $PFX_move($OBJ_t *$OBJ, gmx_enerdata_t *enerd,
    gmx_large_int_t step, int bFirstStep, int bLastStep,
    int bGStat, int bXTC, int bNS,
    t_commrec *cr)
{
  if ($OBJ->mode == 0) {
  }

  return 0;
}'''


  opt2fn_src = r'''
char *$PFX_opt2fn(const char *opt, int nfile, const t_filenm fnm[])
{
  int i;
  for (i = 0; i < nfile; i++) {
    if (strcmp(opt, fnm[i].opt) == 0) { /* a match is found */
      char *fname = fnm[i].fns[0], *p;
      if (fnm[i].ftp == efMDP) {
        /* modify the extension from .mdp to .cfg */
        if (strcmp(fname, "grompp.mdp") == 0) { /* replace the default name by NULL */
          return NULL; /* we do not have default name for .cfg files */
        } else if ((p = strstr(fname, ".cfg.mdp")) != NULL) {
          p[4] = '\0';
        } else if ((p = strstr(fname, ".mdp")) != NULL) {
          strcpy(p, ".cfg");
        }
      }
      return fname;
    }
  }
  return NULL;
}'''

  def getdecl(c):
    decl = c.temprepl(c.struct_src, True)
    if c.isgmx4:
      decl += [ "typedef long int gmx_large_int_t;" ]
    return decl

  def getfuncs(c):
    ''' get basic functions '''
    basic = (c.temprepl(c.init_src, True)
           + c.temprepl(c.finish_src, True)
           + c.temprepl(c.move_src, True)
           + c.temprepl(c.opt2fn_src, True))
    return basic




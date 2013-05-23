#ifndef GMXHMC_H__
#define GMXHMC_H__

/* GROMACS HMC trajectory manipulation
 * since this file is intrinsically GROMACS related
 * we will use GROMACS functions
 * we assume basic functions in zcom.h are included */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec.h"
#include "typedefs.h"
#include "smalloc.h"
#include "mvdata.h"
#include "domdec.h"
#include "partdec.h"
#include "vsite.h"
#include "mdrun.h"


/* assuming a single simulation, no replica exchange */

typedef struct {
  t_state *state;
} gmxhmc_t;

#define COPYRVN(s, t, n) if (s && t) copy_rvecn(s, t, 0, n);

/* copy relavent variables from the state `s' to state `t'
 * cf. exchange_state(), copy_state_nonatomdata() in kernel/repl_ex.c */
INLINE void gmxhmc_copystate(t_state *t, t_state *s)
{
  /* there are generally `s->ngtc' Nose-Hoover chains */
  int ngtc = s->ngtc * s->nhchainlength;
  int nnhpres = s->nnhpres * s->nhchainlength;
  int i;

  COPYRVN(s->box, t->box, DIM);
  COPYRVN(s->box_rel, t->box_rel, DIM);
  COPYRVN(s->boxv, t->boxv, DIM);
  t->veta = s->veta;
  t->vol0 = s->vol0;
  COPYRVN(s->svir_prev, t->svir_prev, DIM);
  COPYRVN(s->fvir_prev, t->fvir_prev, DIM);
  COPYRVN(s->pres_prev, t->pres_prev, DIM);
  for (i = 0; i < ngtc; i++) {
    t->nosehoover_xi[i] = s->nosehoover_xi[i];
    t->nosehoover_vxi[i] = s->nosehoover_vxi[i];
  }
  for (i = 0; i < nnhpres; i++) {
    t->nhpres_xi[i] = s->nhpres_xi[i];
    t->nhpres_vxi[i] = s->nhpres_vxi[i];
  }
  for (i = 0; i < s->ngtc; i++) {
    t->therm_integral[i] = s->therm_integral[i];
  }
  COPYRVN(s->x, t->x, s->natoms);
  COPYRVN(s->v, t->v, s->natoms);
  COPYRVN(s->sd_X, t->sd_X, s->natoms);
}



/* save the current state
 * cf. replica_exchange() in kernel/repl_ex.c */
INLINE void gmxhmc_push(gmxhmc_t *hmc, t_commrec *cr,
    t_state *state, t_state *state_local)
{
  /* collect state to the master node */
  if ( PAR(cr) ) {
    if ( DOMAINDECOMP(cr) ) {
      dd_collect_state(cr->dd, state_local, state);
    } else { /* see pd_collect_state() in kernel/repl_ex.c */
      int shift;
    
      shift = cr->nnodes - cr->npmenodes - 1;
      move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->x,NULL,shift,NULL);
      if (state->v)
        move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->v,NULL,shift,NULL);
      if (state->sd_X)
        move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->sd_X,NULL,shift,NULL);
    }
  }

  /* save the current state to the gmxhmx structure */
  if ( MASTER(cr) ) {
    gmxhmc_copystate(hmc->state, state); /* hmc->state = state */
  }
}


/* revert the current frame to the saved one
 * the prototype is similar to dd_partition_system() */
INLINE void gmxhmc_pop(gmxhmc_t *hmc, FILE *fplog, gmx_large_int_t step,
    t_state *state_global, t_state *state_local, rvec **f,
    gmx_mtop_t *top_global, gmx_localtop_t *top_local,
    t_inputrec *ir, t_commrec *cr, t_mdatoms *mdatoms, t_forcerec *fr,
    gmx_vsite_t *vsite, gmx_shellfc_t shellfc,
    gmx_constr_t constr, t_nrnb *nrnb, gmx_wallcycle_t wcycle)
{
  /* master recovers the global state first */
  if ( MASTER(cr) ) {
    gmxhmc_copystate(state_global, hmc->state);
  }

  /* distribute the changed global state to nonmasters */
  if ( DOMAINDECOMP(cr) ) {
    dd_partition_system(fplog, step, cr, TRUE, 1,
        state_global, top_global, ir, state_local, f, mdatoms,
        top_local, fr, vsite, shellfc, constr,
        nrnb, wcycle, FALSE);
  } else { /* cf. end of replica_exchange() in kernel/repl_ex.c */
    if ( PAR(cr) ) {
      bcast_state(cr, state_global, FALSE);
    }
  }

  /* remember to call neighbor searching in the next step */
}



/* initialize the structure
 * only MASTER(cr) should call this routine */
gmxhmc_t *gmxhmc_open(t_state *st, t_inputrec *ir, t_commrec *cr)
{
  gmxhmc_t *hmc;

  snew(hmc, 1);
  snew(hmc->state, 1);
  /* the following initialization calls are local */
  /* allocate state->x, ->v, ->nosehoover_xi, ->nosehoover_vxi,
   *   ->therm_integral, ->nhpres_xi, ->nhpres_vxi, ... */
  init_state(hmc->state, st->natoms, st->ngtc, st->nnhpres,
      st->nhchainlength);
  /* allocate state->sd_X, ->ld_rng, ->ld_rngi, ... */
  set_state_entries(hmc->state, ir, cr->nnodes);
  return hmc;
}

/* only MASTER(cr) should call this routine */
void gmxhmc_close(gmxhmc_t *hmc)
{
  /* free state->... */
  done_state(hmc->state);
  sfree(hmc->state);
  sfree(hmc);
}

#endif


#ifndef GMXHMC_H__
#define GMXHMC_H__

/* GROMACS HMC trajectory manipulation
 * since this file is intrinsically GROMACS related
 * we will use GROMACS functions
 * we assume basic functions in zcom.h are included */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "domdec.h"
#include "partdec.h"


/* assuming a single simulation, no replica exchange */

typedef struct {
  t_state *state;
} gmxhmc_t;



/* copy relavent variables from the state `s' to state `t'
 * cf. exchange_state() in kernel/repl_ex.c */
INLINE int gmxhmc_copystate(t_state *t, t_state *s)
{
  /* there are generally `s->ngtc' Nose-Hoover chains */
  int ngtc = s->ngtc * s->nhchainlength;
  int nnhpres = s->nnhpres * s->nhchainlength;
  int i;

  copy_mat(s->box, t->box);
  copy_mat(s->box, t->box_rel);
  copy_mat(s->box, t->box_rel);
  t->veta = s->veta;
  t->vol0 = s->vol0;
  copy_mat(s->svir_prev, t->svir_prev);
  copy_mat(s->fvir_prev, t->fvir_prev);
  copy_mat(s->pres_prev, t->pres_prev);
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
  copy_rvecn(s->x, t->x, 0, s->natoms);
  copy_rvecn(s->v, t->v, 0, s->natoms);
  copy_rvecn(s->sd_X, t->sd_X, 0, s->natoms);
}



/* save the current state
 * cf. replica_exchange() in kernel/repl_ex.c */
INLINE int gmxhmc_push(gmxhmc_t *gh, t_commrec *cr,
    t_state *state)
{
  /* collect state to the master node */
  if ( PAR(cr) ) {
    if ( DOMAINDECOMP(cr) ) {
      dd_collect_state(cr);
    } else {
      pd_collect_state(cr, state);
    }
  }

  /* save the current state to the gmxhmx structure */
  if ( MASTER(cr) ) {
    copy_t_state(gh->state, state); /* gh->state = state */
  }
}


/* revert the current frame to the saved one
 * the prototype is similar to dd_partition_system() */
INLINE void gmxhmc_pop(gmxhmc_t *gh, FILE *fplog,
    gmx_large_int_t, t_commrec *cr,
    t_state *state_global, t_inputrec *ir,
    t_state *state_local, rvec **f,
    t_mdatoms *mdatoms, gmx_localtop_t *top_local,
    t_forcerec *fr, gmxvsite_t *vsite, gmx_shellfc_t shellfc,
    gmx_constr_t *constr, t_nrnb *nrnb,
    gmx_wallcycle_t wcycle)
{
  /* master recovers the global state first */
  if ( MASTER(cr) ) {
    gmxhmc_copystate(state_global, gh->state);
  }

  /* distribute the changed global state to nonmasters */
  if ( DOMAINDECOMP(cr) ) {
    dd_partition_system(fplog, step, cr, TRUE, 1,
        state_global, top_global, state, f, mdatoms, top, fr,
        vsite, shellfc, constr, nrnb, wcycle, FALSE);
  } else {
    copy_state_nonatomdata();

    if ( PAR(cr) ) {
      bcast_state(cr, state, FALSE);
    }
  }
  
  /* remember to call neighbor searching in the next step */
}



/* initialize the structure
 * only MASTER(cr) should call this routine */
gmxhmc_t *gmxhmc_open(t_state *st)
{
  gmxhmc_t *gh;

  snew(gh, 1);
  snew(gh, gh->state);
  /* the following initialization calls are local */
  /* allocate state->x, ->v, ->nosehoover_xi, ->nosehoover_vxi,
   *   ->therm_integral, ->nhpres_xi, ->nhpres_vxi, ... */
  init_state(gh->state, st->natoms, st->ngtc, st->nnhpres,
      st->nhchainlength);
  /* allocate state->sd_X, ->ld_rng, ->ld_rngi, ... */
  set_state_entries(gh->state, input, cr->nnodes);
  return gh;
}

/* only MASTER(cr) should call this routine */
void gmxhmc_close(gmxhmc_t *gh)
{
  /* free state->... */
  done_state(gh->state);
  sfree(gh->state);
  sfree(gh);
}

#endif


/*
  Accelerated tempering using integral identities
  adapted from GROMACS package written by
  David van der Spoel, Erik Lindahl, Berk Hess, and others.
  Copyright (C) 2009-2010  Cheng Zhang

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
#include "disre.h"
#include "orires.h"
#include "dihre.h"
#include "pppm.h"
#include "pme.h"
#include "mdatoms.h"
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

#include "gmx_unused.h"

#define GMXVERSION 40007
#include "md_tutil.h"

/* The following two variables and the signal_handler function
 * is used from pme.c as well
 */
extern volatile bool bGotTermSignal, bGotUsr1Signal;

static RETSIGTYPE signal_handler(int n)
{
  switch (n) {
  case SIGTERM:
    bGotTermSignal = TRUE;
    break;
  case SIGUSR1:
    bGotUsr1Signal = TRUE;
    break;
  }
}

time_t md(FILE *fplog, t_commrec *cr, int nfile, t_filenm fnm[],
          bool bVerbose, bool bCompact,
          gmx_vsite_t *vsite, gmx_constr_t constr,
          int stepout, t_inputrec *ir,
          gmx_mtop_t *top_global,
          t_fcdata *fcd,
          t_state *state_global, rvec f[],
          rvec buf[], t_mdatoms *mdatoms,
          t_nrnb *nrnb, gmx_wallcycle_t wcycle,
          t_forcerec *fr,
          real cpt_period, real max_hours,
          unsigned long Flags,
          int *nsteps_done, at_t *at)
{
  int        fp_ene = 0, fp_trn = 0, fp_xtc = 0, step, step_rel, step_ene;
  char       *fn_cpt;
  FILE       *fp_dgdl = NULL, *fp_field = NULL;
  time_t     start_t;
  double     run_time;
  real       t, t0, lam0;
  bool       bGStatEveryStep, bGStat;
  bool       bNS, bSimAnn, bStopCM, bRerunMD, bNotLastFrame = FALSE,
             bFirstStep, bStateFromTPX, bLastStep;
  bool       bNEMD, do_ene, do_log, do_verbose, bRerunWarnNoV = TRUE,
                 bForceUpdate = FALSE, bX, bV, bF, bXTC, bCPT;
  bool       bMasterState;
  tensor     force_vir, shake_vir, total_vir, pres, ekin;
  int        i, m, status;
  rvec       mu_tot;
  t_vcm      *vcm;
  int        step_ns = 0, step_nscheck = 0, nns = 0, nabnsb = 0, ns_lt;
  double     ns_s1 = 0, ns_s2 = 0, ns_ab = 0, ns_lt_runav = 0, ns_lt_runav2 = 0;
  matrix     *scale_tot;
  t_trxframe rerun_fr;
  int        nchkpt = 1;
  /* Booleans (disguised as a reals) to checkpoint and terminate mdrun */
  real       chkpt = 0, terminate = 0, terminate_now = 0;

  gmx_localtop_t *top;
  t_mdebin *mdebin = NULL;
  t_state    *state = NULL;
  rvec       *f_global = NULL;
  gmx_enerdata_t *enerd;
  gmx_stochd_t sd = NULL;
  t_graph    *graph = NULL;

  gmx_groups_t *groups;
  gmx_ekindata_t *ekind;
  gmx_shellfc_t shellfc;
  int         count, nconverged = 0;
  double      tcount = 0;
  bool        bHaveConstr = FALSE;
  bool        bConverged = TRUE, bOK, bSumEkinhOld;
  bool        bAppend;
  real        temp0, dvdl;
  int         a0, a1, ii;
  matrix      lastbox;
  double      cycles;

  /* Check for special mdrun options */
  bRerunMD = (Flags & MD_RERUN);
  bGStatEveryStep = !(Flags & MD_NOGSTAT);
  bAppend  = (Flags & MD_APPENDFILES);

  if (!bGStatEveryStep && !EI_DYNAMICS(ir->eI)) {
    const char *warn = "\nWARNING:\nNo energy summing can only be used with dynamics, ignoring this option\n";
    fprintf(stderr,"%s\n", warn);
    if (fplog)
      fprintf(fplog,"%s\n", warn);
    bGStatEveryStep = TRUE;
  }
  if (!bGStatEveryStep) {
    if (fplog) {
      fprintf(fplog,"\nWill not sum the energies at every step,\n"
              "therefore the energy file does not contain exact averages and fluctuations.\n\n");
      if (ir->etc != etcNO || ir->epc != epcNO) {
        fprintf(fplog,"WARNING:\nThe temperature and/or pressure for scaling will only be updated every nstlist (%d) steps\n\n", ir->nstlist);
      }
      if (ir->nstlist == -1) {
        die_if (at->nsttemp < 0, "must specify a fixed ns frequency.\n");
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  < n.list life time > - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
      }
    }
    if (ir->comm_mode != ecmNO && ir->nstcomm == 1) {
      if (fplog) {
        fprintf(fplog,"WARNING:\nWe should not remove the COM motion every step with option -nosum,\n");
      }
      if (ir->nstlist > 0) {
        ir->nstcomm = ir->nstlist;
        if (fplog) {
          fprintf(fplog,"setting nstcomm to nstlist (%d)\n\n", ir->nstcomm);
        }
      } else {
        ir->nstcomm = 100;
        if (fplog) {
          fprintf(fplog,"setting nstcomm to %d\n\n", ir->nstcomm);
        }
      }
    }
  }

  if (bRerunMD)
    ir->nstxtcout = 0;

  groups = &top_global->groups;

  /* Initial values */
  init_md(fplog, cr, ir, &t, &t0, &state_global->lambda, &lam0,
          nrnb, top_global, &sd,
          nfile, fnm, &fp_trn, &fp_xtc, &fp_ene, &fn_cpt,
          &fp_dgdl, &fp_field, &mdebin,
          force_vir, shake_vir, mu_tot, &bNEMD, &bSimAnn, &vcm, Flags);

  /* Energy terms and groups */
  snew(enerd, 1);
  init_enerdata(fplog, top_global->groups.grps[egcENER].nr, enerd);
  /* Kinetic energy data */
  snew(ekind, 1);
  init_ekindata(fplog, top_global, &(ir->opts), ekind);
  /* Copy the cos acceleration to the groups struct */
  ekind->cosacc.cos_accel = ir->cos_accel;

  /* Check for polarizable models and flexible constraints */
  /* Normally, it does nothing, and directly returns null, when both shell and flexcon are 0  */
  shellfc = init_shell_flexcon(fplog,
                               top_global, n_flexible_constraints(constr),
                               (ir->bContinuation ||
                                (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                               NULL : state_global->x);

  {
    double io = compute_io(ir, top_global->natoms, groups, mdebin->ebin->nener, 1);
    if ((io > 2000) && MASTER(cr))
      fprintf(stderr,
              "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
              io);
  }

  if (DOMAINDECOMP(cr)) {
    top = dd_init_local_top(top_global);

    snew(state, 1);
    dd_init_local_state(cr->dd, state_global, state);

    if (DDMASTER(cr->dd) && ir->nstfout) {
      snew(f_global, state_global->natoms);
    }
  } else {
    if (PAR(cr)) {
      /* Initialize the particle decomposition and split the topology */
      top = split_system(fplog, top_global, ir, cr);

      pd_cg_range(cr, &fr->cg0, &fr->hcg);
      pd_at_range(cr, &a0, &a1);
    } else {
      top = gmx_mtop_generate_local_top(top_global, ir);

      a0 = 0;
      a1 = top_global->natoms;
    }

    state = partdec_init_local_state(cr, state_global);
    f_global = f;

    atoms2md(top_global, ir, 0, NULL, a0, a1-a0, mdatoms);

    if (vsite) {
      set_vsite_top(vsite, top, mdatoms, cr);
    }

    if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
      graph = mk_graph(fplog, &(top->idef), 0, top_global->natoms, FALSE, FALSE);
    }

    if (shellfc) {
      make_local_shells(cr, mdatoms, shellfc);
    }

    if (ir->pull && PAR(cr)) {
      dd_make_local_pull_groups(NULL, ir->pull, mdatoms);
    }
  }

  if (DOMAINDECOMP(cr)) {
    /* Distribute the charge groups over the nodes from the master node */
    dd_partition_system(fplog, ir->init_step, cr, TRUE,
                        state_global, top_global, ir,
                        state, &f, &buf, mdatoms, top, fr, vsite, shellfc, constr,
                        nrnb, wcycle, FALSE);
  }

  update_mdatoms(mdatoms, state->lambda);

  if (MASTER(cr))
  {
    /* Update mdebin with energy history if appending to output files */
    if ( Flags & MD_APPENDFILES )
    {
      restore_energyhistory_from_state(mdebin, &state_global->enerhist);
    }
    /* Set the initial energy history in state to zero by updating once */
    update_energyhistory(&state_global->enerhist, mdebin);
  }


  if (sd && (Flags & MD_READ_RNG)) {
    /* Set the random state if we read a checkpoint file */
    set_stochd_state(sd, state);
  }

  /* Initialize constraints */
  if (constr) {
    if (!DOMAINDECOMP(cr))
      set_constraints(constr, top, ir, mdatoms, NULL);
    bHaveConstr = TRUE;
  }

  if (!ir->bContinuation && !bRerunMD) {
    if (mdatoms->cFREEZE && (state->flags & (1<<estV))) {
      /* Set the velocities of frozen particles to zero */
      for (i = mdatoms->start; i < mdatoms->start+mdatoms->homenr; i++) {
        for (m = 0; m < DIM; m++) {
          if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m]) {
            state->v[i][m] = 0;
          }
        }
      }
    }
    if (bHaveConstr) {
      /* Constrain the initial coordinates and velocities */
      do_shakefirst(fplog, constr, ir, mdatoms, state, buf, f,
                    graph, cr, nrnb, fr, &top->idef);
    }
    if (vsite) {
      /* Construct the virtual sites for the initial configuration */
      construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, NULL,
                       top->idef.iparams, top->idef.il,
                       fr->ePBC, fr->bMolPBC, graph, cr, state->box);
    }
  }

  if (Flags & MD_READ_EKIN)
  {
      restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
  }
  else
  {
      /* Compute initial EKin for all.. */
      if (ekind->cosacc.cos_accel == 0) {
          calc_ke_part(state->v, &(ir->opts), mdatoms, ekind, nrnb, state->lambda);
      } else {
          calc_ke_part_visc(state->box, state->x, state->v, &(ir->opts),
                            mdatoms, ekind, nrnb, state->lambda);
      }
      if (PAR(cr))
      {
          GMX_MPE_LOG(ev_global_stat_start);

          global_stat(fplog, cr, enerd, force_vir, shake_vir, mu_tot,
                      ir, ekind, FALSE, constr, vcm, NULL, NULL, &terminate);

          GMX_MPE_LOG(ev_global_stat_finish);
      }
  }

  /* Calculate the initial half step temperature */
  temp0 = sum_ekin(TRUE, &(ir->opts), ekind, ekin, NULL);

  if (MASTER(cr)) {
    if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
      fprintf(fplog,
              "RMS relative constraint deviation after constraining: %.2e\n",
              constr_rmsd(constr, FALSE));
    fprintf(fplog,"Initial temperature: %g K\n", temp0);
    if (bRerunMD) {
      fprintf(stderr,"starting md rerun '%s', reading coordinates from"
              " input trajectory '%s'\n\n",
              *(top_global->name), opt2fn("-rerun", nfile, fnm));
      if (bVerbose)
        fprintf(stderr,"Calculated time to finish depends on nsteps from "
                "run input file,\nwhich may not correspond to the time "
                "needed to process input trajectory.\n\n");
    } else {
      if (ir->init_step > 0)
      {
         fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps (continuing from step %d, %8.1f ps).\n",
           *(top_global->name), ir->nsteps+ir->init_step, (ir->nsteps+ir->init_step)*ir->delta_t,
           ir->init_step, ir->init_step*ir->delta_t);
      }
      else
      {
         fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n",
           *(top_global->name), ir->nsteps, ir->nsteps*ir->delta_t);
      }
    }
    fprintf(fplog,"\n");
  }

  if (ir->nstlist == -1) {
    snew(scale_tot, 1);
  } else {
    scale_tot = NULL;
  }

  /* Write start time */
  start_t = print_date_and_time(fplog, cr->nodeid,"Started mdrun");
  wallcycle_start(wcycle, ewcRUN);
  if (fplog)
    fprintf(fplog,"\n");

  /* Set the node time counter to 0 after initialisation */
  start_time();
  /***********************************************************
   *
   *             Loop over MD steps
   *
   ************************************************************/

  /* if rerunMD then read coordinates and velocities from input trajectory */
  if (bRerunMD) {
    if (getenv("GMX_FORCE_UPDATE"))
      bForceUpdate = TRUE;

    bNotLastFrame = read_first_frame(&status, opt2fn("-rerun", nfile, fnm),
                                     &rerun_fr, TRX_NEED_X | TRX_READ_V);
    if (rerun_fr.natoms != top_global->natoms)
      gmx_fatal(FARGS,"Number of atoms in trajectory (%d) does not match the "
                "run input file (%d)\n", rerun_fr.natoms, top_global->natoms);
    if (ir->ePBC != epbcNONE) {
      if (!rerun_fr.bBox)
        gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used", rerun_fr.step, rerun_fr.time);
      if (max_cutoff2(ir->ePBC, rerun_fr.box) < sqr(fr->rlistlong))
        gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions", rerun_fr.step, rerun_fr.time);

      /* Set the shift vectors.
       * Necessary here when have a static box different from the tpr box.
       */
      calc_shifts(rerun_fr.box, fr->shift_vec);
    }
  }

  /* loop over MD steps or if rerunMD to end of input trajectory */
  bFirstStep = TRUE;
  /* Skip the first Nose-Hoover integration when we get the state from tpx */
  bStateFromTPX = !opt2bSet("-cpi", nfile, fnm);
  bLastStep = FALSE;
  bSumEkinhOld = FALSE,

  step = ir->init_step;
  step_rel = 0;
  step_ene = 0;

  bLastStep = (bRerunMD || step_rel > ir->nsteps);
  while (!bLastStep || (bRerunMD && bNotLastFrame)) {

    wallcycle_start(wcycle, ewcSTEP);

    GMX_MPE_LOG(ev_timestep1);

    if (bRerunMD) {
      if (rerun_fr.bStep) {
        step = rerun_fr.step;
        step_rel = step - ir->init_step;
      }
      if (rerun_fr.bTime)
        t = rerun_fr.time;
      else
        t = (real)step;
    } else {
      bLastStep = (step_rel == ir->nsteps);

      t = t0 + step*ir->delta_t;
    }
    if (Flags & MD_APPENDFILES) {
      step_ene = step;
    } else {
      step_ene = step_rel;
    }

    if (bRerunMD) {
      if (!(DOMAINDECOMP(cr) && !MASTER(cr))) {
        for (i = 0; i < state_global->natoms; i++) {
          copy_rvec(rerun_fr.x[i], state_global->x[i]);
        }
        if (rerun_fr.bV) {
          for (i = 0; i < state_global->natoms; i++) {
            copy_rvec(rerun_fr.v[i], state_global->v[i]);
          }
        } else {
          for (i = 0; i < state_global->natoms; i++) {
            clear_rvec(state_global->v[i]);
          }
          if (bRerunWarnNoV) {
            fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                    "         Ekin, temperature and pressure are incorrect,\n"
                    "         the virial will be incorrect when constraints are present.\n"
                    "\n");
            bRerunWarnNoV = FALSE;
          }
        }
      }
      copy_mat(rerun_fr.box, state_global->box);
      copy_mat(state_global->box, state->box);

      if (vsite && (Flags & MD_RERUN_VSITE)) {
        if (DOMAINDECOMP(cr)) {
          gmx_fatal(FARGS,"Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
        }
        if (graph) {
          /* Following is necessary because the graph may get out of sync
           * with the coordinates if we only have every N'th coordinate set
           */
          mk_mshift(fplog, graph, fr->ePBC, state->box, state->x);
          shift_self(graph, state->box, state->x);
        }
        construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
                         top->idef.iparams, top->idef.il,
                         fr->ePBC, fr->bMolPBC, graph, cr, state->box);
        if (graph)
          unshift_self(graph, state->box, state->x);
      }
    }

    /* Stop Center of Mass motion */
    bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

    bNS = bFirstStep;
    if (bRerunMD) {
      /* for rerun MD always do Neighbour Searching */
      if (ir->nstlist != 0) {
        bNS = TRUE;
      }
    } else {
      /* Determine whether or not to do Neighbour Searching */
      if (ir->nstlist > 0 && (step % ir->nstlist == 0)) {
        bNS = TRUE;
      } else if (ir->nstlist == -1) {
        bNS = (bFirstStep || nabnsb > 0);
        if (bNS) {
          if (bFirstStep) {
            ns_lt_runav = 0;
            step_nscheck = step;
          } else {
            /* Determine the neighbor list life time */
            ns_lt = step - step_ns;
            if (debug) {
              fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %d steps\n", nabnsb, ns_lt);
            }
            nns++;
            ns_s1 += ns_lt;
            ns_s2 += ns_lt*ns_lt;
            ns_ab += nabnsb;
            if (ns_lt_runav == 0) {
              ns_lt_runav  = ns_lt;
              /* Initialize the fluctuation average such that at startup
               * we check after 0 steps.
               */
              ns_lt_runav2 = dblsqr(ns_lt/2.0);
            }
            /* Running average with 0.9 gives an exp. history of 9.5 */
            ns_lt_runav2 = 0.9*ns_lt_runav2 + 0.1*dblsqr(ns_lt_runav - ns_lt);
            ns_lt_runav  = 0.9*ns_lt_runav  + 0.1*ns_lt;
            if (bGStatEveryStep) {
              /* Always check the nlist validity */
              step_nscheck = step;
            } else {
              /* We check after:  <life time> - 2*sigma
               * The factor 2 is quite conservative,
               * but we assume that with nstlist=-1 the user prefers
               * exact integration over performance.
               */
              step_nscheck = step
                + (int)(ns_lt_runav - 2.0*sqrt(ns_lt_runav2)) - 1;
            }
            if (debug) {
              fprintf(debug,"nlist life time %d run av. %4.1f sig %3.1f check %d check with -nosum %d\n",
                      ns_lt, ns_lt_runav, sqrt(ns_lt_runav2),
                      step_nscheck-step+1,
                      (int)(ns_lt_runav - 2.0*sqrt(ns_lt_runav2)));
            }
          }
          step_ns = step;
          /* Initialize the cumulative coordinate scaling matrix */
          clear_mat(*scale_tot);
          for (ii = 0; ii < DIM; ii++)
            (*scale_tot)[ii][ii] = 1.0;
        }
      }
    }

    if (terminate_now > 0 || (terminate_now < 0 && bNS)) {
      bLastStep = TRUE;
    }

    do_log = do_per_step(step, ir->nstlog) || bFirstStep || bLastStep;
    do_verbose = bVerbose && (step % stepout == 0 || bFirstStep || bLastStep);

    if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD)) {
      if (bRerunMD) {
        bMasterState = TRUE;
      } else {
        bMasterState = FALSE;
        /* Correct the new box if it is too skewed */
        if (DYNAMIC_BOX(*ir)) {
          if (correct_box(fplog, step, state->box, graph))
            bMasterState = TRUE;
        }
        if (DOMAINDECOMP(cr) && bMasterState)
          dd_collect_state(cr->dd, state, state_global);
      }

      if (DOMAINDECOMP(cr)) {
        /* Repartition the domain decomposition */
        wallcycle_start(wcycle, ewcDOMDEC);
        dd_partition_system(fplog, step, cr, bMasterState,
                            state_global, top_global, ir,
                            state, &f, &buf, mdatoms, top, fr, vsite, shellfc, constr,
                            nrnb, wcycle, do_verbose);
        wallcycle_stop(wcycle, ewcDOMDEC);
      }
    }

    if (MASTER(cr) && do_log)
      print_ebin_header(fplog, step, t, state->lambda);

    if (bRerunMD && rerun_fr.bV) {
      /* We need the kinetic energy at minus the half step for determining
       * the full step kinetic energy and possibly for T-coupling.
       */
      calc_ke_part(state->v, &(ir->opts), mdatoms, ekind, nrnb, state->lambda);
      if (PAR(cr)) {
        global_stat(fplog, cr, enerd, force_vir, shake_vir, mu_tot,
                    ir, ekind, FALSE, constr, vcm, NULL, NULL, &terminate);
      }
      sum_ekin(FALSE, &(ir->opts), ekind, ekin, NULL);
    }

    clear_mat(force_vir);

    GMX_MPE_LOG(ev_timestep2);

    if (shellfc) {
      /* Now is the time to relax the shells */
      count = relax_shell_flexcon(fplog, cr, bVerbose, step,
                                ir, bNS, bStopCM, top, constr, enerd, fcd,
                                state, f, buf, force_vir, mdatoms,
                                nrnb, wcycle, graph, groups,
                                shellfc, fr, t, mu_tot,
                                state->natoms, &bConverged, vsite,
                                fp_field);
      tcount += count;

      if (bConverged)
        nconverged++;
    }
    else {
      /* The coordinates (x) are shifted (to get whole molecules) in do_force
       * This is parallellized as well, and does communication too.
       * Check comments in sim_util.c
       */
      atgmx_doforce(fplog, cr, ir, step, nrnb, wcycle, top, groups,
               state->box, state->x, &state->hist,
               f, buf, force_vir, mdatoms, enerd, fcd,
               state->lambda, graph,
               fr, vsite, mu_tot, t, fp_field, NULL,
               GMX_FORCE_STATECHANGED | (bNS ? GMX_FORCE_NS : 0) |
	       GMX_FORCE_ALLFORCES, at);
    }

    GMX_BARRIER(cr->mpi_comm_mygroup);
    /* Now we have the energies and forces corresponding to the
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */

    GMX_MPE_LOG(ev_output_start);

    bX   = do_per_step(step, ir->nstxout);
    bV   = do_per_step(step, ir->nstvout);
    bF   = do_per_step(step, ir->nstfout);
    bXTC = do_per_step(step, ir->nstxtcout);
    if ((bNS || bLastStep) && (step > ir->init_step) && !bRerunMD) {
      bCPT = ((chkpt < 0 && do_per_step(step, ir->nstenergy)) || chkpt > 0 ||
               bLastStep);
      if (bCPT) {
        chkpt = 0;
      }
    } else {
      bCPT = FALSE;
    }

    if (bX || bV || bF || bXTC || bCPT) {
      wallcycle_start(wcycle, ewcTRAJ);
      if (bCPT) {
        if (sd) {
          get_stochd_state(sd, state);
        }
        if (MASTER(cr)) {
          if (bSumEkinhOld) {
            state_global->ekinstate.bUpToDate = FALSE;
          } else {
            update_ekinstate(&state_global->ekinstate, ekind);
            state_global->ekinstate.bUpToDate = TRUE;
          }
          update_energyhistory(&state_global->enerhist, mdebin);
        }
      }
      write_traj(fplog, cr, fp_trn, bX, bV, bF, fp_xtc, bXTC, (int)(ir->xtcprec), fn_cpt, bCPT,
                                 top_global, ir->eI, ir->simulation_part, step, t, state, state_global, f, f_global);
      if (bLastStep && step_rel == ir->nsteps &&
          (Flags & MD_CONFOUT) && MASTER(cr) &&
          !bRerunMD) {
        /* x and v have been collected in write_traj */
        fprintf(stderr,"\nWriting final coordinates.\n");
        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols && DOMAINDECOMP(cr)) {
          /* Make molecules whole only for confout writing */
          do_pbc_mtop(fplog, ir->ePBC, state->box, top_global, state_global->x);
        }
        write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                            *top_global->name, top_global,
                            state_global->x, state_global->v,
                            ir->ePBC, state->box);
      }
      wallcycle_stop(wcycle, ewcTRAJ);
    }
    GMX_MPE_LOG(ev_output_finish);
    clear_mat(shake_vir);

    /* Box is changed in update() when we do pressure coupling,
     * but we should still use the old box for energy corrections and when
     * writing it to the energy file, so it matches the trajectory files for
     * the same timestep above. Make a copy in a separate array.
     */
    copy_mat(state->box, lastbox);

    GMX_MPE_LOG(ev_update_start);
    /* This is also parallellized, but check code in update.c */
    /* bOK = update(nsb->natoms, START(nsb), HOMENR(nsb), step, state->lambda, &ener[F_DVDL], */
    bOK = TRUE;
    if (!bRerunMD || rerun_fr.bV || bForceUpdate) {
      wallcycle_start(wcycle, ewcUPDATE);
      dvdl = 0;
      /* We can only do Berendsen coupling after we have summed the kinetic
       * energy or virial. Since the happens in global_stat after update,
       * we should only do it at step % nstlist = 1 with bGStatEveryStep=FALSE.
       */
      update(fplog, step, &dvdl, ir, mdatoms, state, graph, f, buf, fcd,
             &top->idef, ekind, shake_vir, scale_tot,
             cr, nrnb, wcycle, sd, constr, bHaveConstr,
             bNEMD, bFirstStep && bStateFromTPX);
      if (fr->bSepDVDL && fplog && do_log) {
        fprintf(fplog, sepdvdlformat,"Constraint", 0.0, dvdl);
      }
      enerd->term[F_DGDL_CON] += dvdl;
      wallcycle_stop(wcycle, ewcUPDATE);
    } else if (graph) {
      /* Need to unshift here */
      unshift_self(graph, state->box, state->x);
    }

    GMX_BARRIER(cr->mpi_comm_mygroup);
    GMX_MPE_LOG(ev_update_finish);

    if (!bOK)
      gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");

    if (vsite) {
      wallcycle_start(wcycle, ewcVSITECONSTR);
      if (graph)
        shift_self(graph, state->box, state->x);

      construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
                       top->idef.iparams, top->idef.il,
                       fr->ePBC, fr->bMolPBC, graph, cr, state->box);

      if (graph)
        unshift_self(graph, state->box, state->x);
      wallcycle_stop(wcycle, ewcVSITECONSTR);
    }

    /* Non-equilibrium MD: this is parallellized, but only does communication
     * when there really is NEMD.
     */
    if (PAR(cr) && bNEMD)
      accumulate_u(cr, &(ir->opts), ekind);

    if (ekind->cosacc.cos_accel == 0) {
      calc_ke_part(state->v, &(ir->opts), mdatoms, ekind, nrnb, state->lambda);
    } else {
      calc_ke_part_visc(state->box, state->x, state->v, &(ir->opts),
                        mdatoms, ekind, nrnb, state->lambda);
    }

    /* since we use the new coordinates in calc_ke_part_visc, we should use
     * the new box too. Still, won't this be offset by one timestep in the
     * energy file? / EL 20040121
     */
    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM && !bRerunMD)
      calc_vcm_grp(fplog, mdatoms->start, mdatoms->homenr, mdatoms,
                   state->x, state->v, vcm);

    /* Determine the wallclock run time up till now */
    run_time = (double)time(NULL) - (double)start_t;

    /* Check whether everything is still allright */
    if (bGotTermSignal || bGotUsr1Signal) {
      if (bGotTermSignal || ir->nstlist == 0)
        terminate = 1;
      else
        terminate = -1;
      if (!PAR(cr))
        terminate_now = terminate;
      if (fplog) {
        fprintf(fplog,
                "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                bGotTermSignal ? "TERM" : "USR1", terminate == -1 ? "NS " : "");
        fflush(fplog);
      }
      fprintf(stderr,
              "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
              bGotTermSignal ? "TERM" : "USR1", terminate == -1 ? "NS " : "");
      fflush(stderr);
      bGotTermSignal = FALSE;
      bGotUsr1Signal = FALSE;
    } else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
               (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
               terminate == 0) {
      /* Signal to terminate the run */
      terminate = (ir->nstlist == 0 ? 1.0f : -1.0f);
      if (!PAR(cr))
        terminate_now = terminate;
     if (fplog)
        fprintf(fplog,"\nStep %d: Run time exceeded %.3f hours, will terminate the run\n", step, max_hours*0.99);
      fprintf(stderr, "\nStep %d: Run time exceeded %.3f hours, will terminate the run\n", step, max_hours*0.99);
    }

    bGStat = (bGStatEveryStep || bStopCM);

    if (ir->nstlist == -1 && !bRerunMD) {
      /* When bGStatEveryStep=FALSE, global_stat is only called
       * when we check the atom displacements, not at NS steps.
       * This means that also the bonded interaction count check is not
       * performed immediately after NS. Therefore a few MD steps could
       * be performed with missing interactions. But wrong energies are never
       * written to file, since energies are only written after global_stat
       * has been called.
       */
      if (step >= step_nscheck) {
        nabnsb = natoms_beyond_ns_buffer(ir, fr, &top->cgs, *scale_tot, state->x);
        bGStat = TRUE;
      } else {
        /* This is not necessarily true,
         * but step_nscheck is determined quite conservatively.
         */
        nabnsb = 0;
      }
    } else {
      if (bNS) {
        bGStat = TRUE;
      }
    }

    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication, otherwise the other nodes don't know.
     */
    if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                       cpt_period >= 0 &&
                       (cpt_period == 0 ||
                        run_time >= nchkpt*cpt_period*60.0))) {
      if (chkpt == 0) {
        nchkpt++;
      }
      /* Write checkpoint at the next energy output step (if there is one),
       * or after 0.2*cpt_period at any step.
       */
      if (!bGStatEveryStep || ir->nstenergy == 0 || cpt_period == 0 ||
          run_time >= (nchkpt + 0.2)*cpt_period*60.0) {
        chkpt = 1;
      } else {
        chkpt = -1;
      }
    }

    /* With exact energy averages (bGStatEveryStep=TRUE)
     * we should also write energy at first, last and continuation steps
     * such that we can get exact averages over a series of runs.
     * We therefore try to checkpoint at energy output frames.
     *
     * This is not necessary when we use the append-file-feature, so we avoid
     * the extra first frame in that case.
     */
    do_ene = (do_per_step(step, ir->nstenergy) ||
              (bGStatEveryStep && ((bFirstStep && !bAppend) ||
                                   bLastStep || bCPT)));

    if (do_ene || do_log) {
      bGStat = TRUE;
    }

    if (!bGStat) {
      /* We will not sum ekinh_old, so signal that we still have to do it */
      bSumEkinhOld = TRUE;
    } else {
      if (PAR(cr)) {
        wallcycle_start(wcycle, ewcMoveE);
        /* Globally (over all NODEs) sum energy, virial etc.
         * This includes communication
         */
        global_stat(fplog, cr, enerd, force_vir, shake_vir, mu_tot,
                    ir, ekind, bSumEkinhOld, constr, vcm,
                    ir->nstlist == -1 ? &nabnsb : NULL, &chkpt, &terminate);
        if (terminate != 0) {
          terminate_now = terminate;
          terminate = 0;
        }

        wallcycle_stop(wcycle, ewcMoveE);
        bSumEkinhOld = FALSE;
      }

      /* This is just for testing. Nothing is actually done to Ekin
       * since that would require extra communication.
       */
      if (!bNEMD && debug && (vcm->nr > 0)) {
        correct_ekin(debug, mdatoms->start, mdatoms->start+mdatoms->homenr,
                     state->v, vcm->group_p[0],
                     mdatoms->massT, mdatoms->tmass, ekin);
      }

      /* Do center of mass motion removal */
      if (bStopCM && !bRerunMD) {
        check_cm_grp(fplog, vcm, 1);
        do_stopcm_grp(fplog, mdatoms->start, mdatoms->homenr, mdatoms->cVCM,
                      state->x, state->v, vcm);
        inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);
        /*
          calc_vcm_grp(fplog,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);
          check_cm_grp(fplog,vcm);
          do_stopcm_grp(fplog,START(nsb),HOMENR(nsb),x,v,vcm);
          check_cm_grp(fplog,vcm);
        */
      }

      /* Add force and shake contribution to the virial */
      m_add(force_vir, shake_vir, total_vir);

      /* Calculate the amplitude of the cosine velocity profile */
      ekind->cosacc.vcos = ekind->cosacc.mvcos/mdatoms->tmass;

      /* Sum the kinetic energies of the groups & calc temp */
      enerd->term[F_TEMP] = sum_ekin((bRerunMD && !rerun_fr.bV),
                                     &(ir->opts), ekind, ekin,
                                     &(enerd->term[F_DKDL]));
      enerd->term[F_EKIN] = trace(ekin);

      /* Calculate pressure and apply LR correction if PPPM is used.
       * Use the box from last timestep since we already called update().
       */
      enerd->term[F_PRES] =
        calc_pres(fr->ePBC, ir->nwall, lastbox, ekin, total_vir, pres,
                  (fr->eeltype == eelPPPM) ? enerd->term[F_COUL_RECIP] : 0.0f);

      /* Calculate long range corrections to pressure and energy */
      /* based on average C6 and average C12, disabled by default */
      calc_dispcorr(fplog, ir, fr, step, top_global->natoms,
                    lastbox, state->lambda,
                    pres, total_vir, enerd->term);

      enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

      switch (ir->etc) {
      case etcNO:
      case etcBERENDSEN:
        break;
      case etcNOSEHOOVER:
        enerd->term[F_ECONSERVED] =
          enerd->term[F_ETOT] + nosehoover_energy(&(ir->opts), ekind,
                                                  state->nosehoover_xi,
                                                  state->therm_integral);
        break;
      case etcVRESCALE:
        enerd->term[F_ECONSERVED] =
          enerd->term[F_ETOT] + vrescale_energy(&(ir->opts),
                                                state->therm_integral);
        break;
      }

      /* Complicated conditional when bGStatEveryStep=FALSE.
       * We can not just use bGStat, since then the simulation results
       * would depend on nstenergy and nstlog or step_nscheck.
       */
      if ((state->flags & (1<<estPRES_PREV)) &&
          (bGStatEveryStep ||
           (ir->nstlist > 0 && step % ir->nstlist == 0) ||
           (ir->nstlist < 0 && nabnsb > 0) ||
           (ir->nstlist == 0 && bGStat))) {
        /* Store the pressure in t_state for pressure coupling
         * at the next MD step.
         */
        copy_mat(pres, state->pres_prev);
      }
    }

    /* update data, temperature and change at->scale */
    if (0 != atgmx_move(at, enerd,
         step, bFirstStep, bLastStep, bGStat,
         bXTC, bNS, cr) ) {
      exit(1);
    }

    /* Time for performance */
    if (((step % stepout) == 0) || bLastStep)
      update_time();

    /* Output stuff */
    if (MASTER(cr)) {
      bool do_dr, do_or;

      /* If we do reruns, the step numbers in the output energy frames
       * cannot be used for averages (since energies are only calculated
       * for trajectory frames).
       */
       upd_mdebin(mdebin, fp_dgdl, bGStatEveryStep && !bRerunMD,
                 mdatoms->tmass, step_ene, t, enerd, state, lastbox,
                 shake_vir, force_vir, total_vir, pres,
                 ekind, mu_tot, constr);

      do_dr  = do_per_step(step, ir->nstdisreout);
      do_or  = do_per_step(step, ir->nstorireout);

      print_ebin(fp_ene, do_ene, do_dr, do_or, do_log ? fplog : NULL, step, step_ene, t,
                 eprNORMAL, bCompact, mdebin, fcd, groups, &(ir->opts));

      if (ir->ePull != epullNO)
        {
          pull_print_output(ir->pull, step, t);
        }

      if (do_per_step(step, ir->nstlog))
        {
                if (fflush(fplog) != 0)
                {
                        gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
                }
        }
    }

    /* Remaining runtime */
    if (MULTIMASTER(cr) && do_verbose) {
      if (shellfc)
        fprintf(stderr,"\n");
      print_time(stderr, start_t, step, ir);
    }

    bFirstStep = FALSE;

    if (bRerunMD)
      /* read next frame from input trajectory */
      bNotLastFrame = read_next_frame(status, &rerun_fr);

    if (!bRerunMD || !rerun_fr.bStep) {
      /* increase the MD step number */
      step++;
      step_rel++;
    }

    cycles = wallcycle_stop(wcycle, ewcSTEP);
    if (DOMAINDECOMP(cr) && wcycle)
      dd_cycles_add(cr->dd, (float)cycles, ddCyclStep);
  }
  /* End of main MD loop */

  if (bRerunMD)
    close_trj(status);

  if (!(cr->duty & DUTY_PME)) {
    /* Tell the PME only node to finish */
    gmx_pme_finish(cr);
  }

  if (MASTER(cr)) {
    if (bGStatEveryStep && !bRerunMD) {
      print_ebin(fp_ene, FALSE, FALSE, FALSE, fplog, step, step_ene, t,
                 eprAVER, FALSE, mdebin, fcd, groups, &(ir->opts));
      print_ebin(fp_ene, FALSE, FALSE, FALSE, fplog, step, step_ene, t,
                 eprRMS, FALSE, mdebin, fcd, groups, &(ir->opts));
    }
    close_enx(fp_ene);
    if (ir->nstxtcout)
      close_xtc(fp_xtc);
    close_trn(fp_trn);
    if (fp_dgdl)
      gmx_fio_fclose(fp_dgdl);
    if (fp_field)
      gmx_fio_fclose(fp_field);
  }

  if (ir->nstlist == -1 && nns > 0 && fplog) {
    fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n", ns_s1/nns, sqrt(ns_s2/nns - dblsqr(ns_s1/nns)));
    fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n", ns_ab/nns);
  }

  if (shellfc && fplog) {
    fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
            (nconverged*100.0)/step_rel);
    fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
            tcount/step_rel);
  }

  *nsteps_done = step_rel;

  return start_t;
}


static
void runner(FILE *fplog, t_commrec *cr, int nfile, t_filenm fnm[],
    bool bVerbose, bool bCompact,
    ivec ddxyz, int dd_node_order, real rdd, real rconstr,
    const char *dddlb_opt, real dlb_scale,
    const char *ddcsx, const char *ddcsy, const char *ddcsz,
    int nstepout,
    real pforce, real cpt_period, real max_hours,
    unsigned long Flags, int atpremode)
{
  double     nodetime = 0, realtime;
  t_inputrec *inputrec;
  t_state    *state = NULL;
  matrix     box;
  rvec       *buf = NULL, *f = NULL;
  t_nrnb     *nrnb;
  gmx_mtop_t *mtop = NULL;
  t_mdatoms  *mdatoms = NULL;
  t_forcerec *fr = NULL;
  t_fcdata   *fcd = NULL;
  real       ewaldcoeff = 0;
  gmx_pme_t  *pmedata = NULL;
  time_t     start_t = 0;
  gmx_vsite_t *vsite = NULL;
  gmx_constr_t constr;
  int        nChargePerturbed = -1, status;
  gmx_wallcycle_t wcycle;
  bool       bReadRNG, bReadEkin;
  int        list;
  int        nsteps_done = 0;
  at_t *at;

  /* say hello to global unused variables */
  CALL_UNUSED_VEC_H();
  CALL_UNUSED_STATUTIL_H();
  CALL_UNUSED_PME_H();

  snew(inputrec, 1);
  snew(mtop, 1);

  if (Flags & MD_APPENDFILES)
  {
          fplog = NULL;
  }

  if (PAR(cr)) {
    /* The master thread on the master node reads from disk,
     * then distributes everything to the other processors.
     */

    list = (SIMMASTER(cr) && !(Flags & MD_APPENDFILES)) ?  (LIST_SCALARS | LIST_INPUTREC) : 0;

    snew(state, 1);
    init_parallel(fplog, ftp2fn(efTPX, nfile, fnm), cr,
                  inputrec, mtop, state, list);

  } else {
    /* Read a file for a single processor */
    snew(state, 1);
    init_single(fplog, inputrec, ftp2fn(efTPX, nfile, fnm), mtop, state);
  }
  if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC)) {
    cr->npmenodes = 0;
  }

  /* NMR restraints must be initialized before load_checkpoint,
   * since with time averaging the history is added to t_state.
   * For proper consistency check we therefore need to extend
   * t_state here.
   * So the PME-only nodes (if present) will also initialize
   * the distance restraints.
   */
  snew(fcd, 1);

  /* This needs to be called before read_checkpoint to extend the state */
  init_disres(fplog, mtop, inputrec, cr, Flags & MD_PARTDEC, fcd, state);

  if (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0) {
    if (PAR(cr) && !(Flags & MD_PARTDEC)) {
      gmx_fatal(FARGS,"Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
    }
    /* Orientation restraints */
    if (MASTER(cr)) {
      init_orires(fplog, mtop, state->x, inputrec, cr->ms, &(fcd->orires), state);
    }
  }

  if (DEFORM(*inputrec)) {
    /* Store the deform reference box before reading the checkpoint */
    if (SIMMASTER(cr)) {
      copy_mat(state->box, box);
    }
    if (PAR(cr)) {
      gmx_bcast(sizeof(box), box, cr);
    }
    set_deform_reference_box(inputrec->init_step, box);
  }

  if (opt2bSet("-cpi", nfile, fnm))
  {
      /* Check if checkpoint file exists before doing continuation.
       * This way we can use identical input options for the first and subsequent runs...
       */
      if (fexist(opt2fn("-cpi", nfile, fnm)) )
      {
          load_checkpoint(opt2fn("-cpi", nfile, fnm), fplog,
                          cr, Flags & MD_PARTDEC, ddxyz,
                          inputrec, state, &bReadRNG, &bReadEkin,
                          (Flags & MD_APPENDFILES));

          if (bReadRNG)
          {
              Flags |= MD_READ_RNG;
          }
          if (bReadEkin)
          {
              Flags |= MD_READ_EKIN;
          }
      }
  }

  if ((MASTER(cr) || (Flags & MD_SEPPOT)) && (Flags & MD_APPENDFILES))
  {
      fplog = gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr, !(Flags & MD_SEPPOT) , Flags);
  }

  if (SIMMASTER(cr))
  {
      copy_mat(state->box, box);
  }

  if (PAR(cr))
  {
      gmx_bcast(sizeof(box), box, cr);
  }

  if (PAR(cr) && !((Flags & MD_PARTDEC) || EI_TPI(inputrec->eI))) {
    cr->dd = init_domain_decomposition(fplog, cr, Flags, ddxyz, rdd, rconstr,
                                       dddlb_opt, dlb_scale,
                                       ddcsx, ddcsy, ddcsz,
                                       mtop, inputrec,
                                       box, state->x);

    make_dd_communicators(fplog, cr, dd_node_order);

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);
  } else {
    cr->duty = (DUTY_PP | DUTY_PME);

    if (inputrec->ePBC == epbcSCREW)
      gmx_fatal(FARGS,"pbc=%s is only implemented with domain decomposition",
                epbc_names[inputrec->ePBC]);
  }

  if (PAR(cr)) {
    /* After possible communicator splitting in make_dd_communicators.
     * we can set up the intra/inter node communication.
     */
    gmx_setup_nodecomm(fplog, cr);
  }

  wcycle = wallcycle_init(fplog, cr);

  /* initialize aT, cr->duty PP/PME has been assigned */
  at = atgmx_init(atgmx_opt2fn("-at", nfile, fnm),
      Flags & MD_STARTFROMCPT,
      mtop, inputrec, cr, atpremode);
  if ((cr->duty & DUTY_PP) && at == NULL) {
    fprintf(stderr, "Failed to initialize aT\n");
    exit(1);
  }

  snew(nrnb, 1);
  if (cr->duty & DUTY_PP) {
    /* For domain decomposition we allocate dynamically
     * in dd_partition_system.
     */
    if (DOMAINDECOMP(cr)) {
      bcast_state_setup(cr, state);
    } else {
      if (PAR(cr)) {
        if (!MASTER(cr)) {
          snew(state, 1);
        }
        bcast_state(cr, state, TRUE);
      }

      snew(buf, mtop->natoms);
      snew(f, mtop->natoms);
    }

    /* Dihedral Restraints */
    if (gmx_mtop_ftype_count(mtop, F_DIHRES) > 0) {
      init_dihres(fplog, mtop, inputrec, fcd);
    }

    /* Initiate forcerecord */
    fr = mk_forcerec();
    init_forcerec(fplog, fr, fcd, inputrec, mtop, cr, box, FALSE,
                  opt2fn("-table", nfile, fnm), opt2fn("-tablep", nfile, fnm),
                  opt2fn("-tableb", nfile, fnm), FALSE, pforce);
    fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);


    /* Initialize the mdatoms structure.
     * mdatoms is not filled with atom data,
     * as this can not be done now with domain decomposition.
     */
    mdatoms = init_mdatoms(fplog, mtop, inputrec->efep != efepNO);

    /* Initialize the virtual site communication */
    vsite = init_vsite(mtop, cr);

    calc_shifts(box, fr->shift_vec);

    /* With periodic molecules the charge groups should be whole at start up
     * and the virtual sites should not be far from their proper positions.
     */
    if (!inputrec->bContinuation && MASTER(cr) &&
        !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols)) {
      /* Make molecules whole at start of run */
      if (fr->ePBC != epbcNONE)  {
       	do_pbc_first_mtop(fplog, inputrec->ePBC, box, mtop, state->x);
      }
      if (vsite) {
        /* Correct initial vsite positions are required
         * for the initial distribution in the domain decomposition
         * and for the initial shell prediction.
         */
        construct_vsites_mtop(fplog, vsite, mtop, state->x);
      }
    }
    /* Initiate PPPM if necessary */
    if (fr->eeltype == eelPPPM) {
      if (mdatoms->nChargePerturbed)
        gmx_fatal(FARGS,"Free energy with %s is not implemented",
                  eel_names[fr->eeltype]);
      status = gmx_pppm_init(fplog, cr, FALSE, TRUE, box,
                             getenv("GMXGHAT"), inputrec, (Flags & MD_REPRODUCIBLE));
      if (status != 0)
        gmx_fatal(FARGS,"Error %d initializing PPPM", status);
    }

    if (EEL_PME(fr->eeltype)) {
      ewaldcoeff = fr->ewaldcoeff;
      pmedata = &fr->pmedata;
    } else {
      pmedata = NULL;
    }
  } else {
    /* This is a PME only node */

    /* We don't need the state */
    done_state(state);

    ewaldcoeff = calc_ewaldcoeff(inputrec->rcoulomb, inputrec->ewald_rtol);
    snew(pmedata, 1);
  }

  /* Initiate PME if necessary,
   * either on all nodes or on dedicated PME nodes only. */
  if (EEL_PME(inputrec->coulombtype)) {
    if (mdatoms) {
      nChargePerturbed = mdatoms->nChargePerturbed;
    }
    if (cr->npmenodes > 0) {
      /* The PME only nodes need to know nChargePerturbed */
      gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, cr);
    }
    if (cr->duty & DUTY_PME) {
      status = gmx_pme_init(pmedata, cr, inputrec,
                            mtop ? mtop->natoms : 0, nChargePerturbed,
                            (Flags & MD_REPRODUCIBLE));
      if (status != 0)
        gmx_fatal(FARGS,"Error %d initializing PME", status);
    }
  }

  {
    /* Turn on signal handling on all nodes */
    /*
     * (A user signal from the PME nodes (if any)
     * is communicated to the PP nodes.
     */
    if (getenv("GMX_NO_TERM") == NULL) {
      if (debug)
        fprintf(debug,"Installing signal handler for SIGTERM\n");
      signal(SIGTERM, signal_handler);
    }
    if (getenv("GMX_NO_USR1") == NULL) {
      if (debug)
        fprintf(debug,"Installing signal handler for SIGUSR1\n");
      signal(SIGUSR1, signal_handler);
    }
  }

  if (cr->duty & DUTY_PP) {
    if (inputrec->ePull != epullNO) {
      /* Initialize pull code */
      init_pull(fplog, inputrec, nfile, fnm, mtop, cr,
                EI_DYNAMICS(inputrec->eI) && MASTER(cr), Flags);
    }

    constr = init_constraints(fplog, mtop, inputrec, NULL, state, cr);

    if (DOMAINDECOMP(cr)) {
      dd_init_bondeds(fplog, cr->dd, mtop, vsite, constr, inputrec,
                      Flags & MD_DDBONDCHECK, fr->cginfo_global);

      set_dd_parameters(fplog, cr->dd, dlb_scale, inputrec, fr, box);

      setup_dd_grid(fplog, cr->dd);
    }

    /* Now do whatever the user wants us to do (how flexible...) */
    start_t = md(fplog, cr, nfile, fnm,
        bVerbose, bCompact,
        vsite, constr,
        nstepout, inputrec, mtop,
        fcd, state, f, buf,
        mdatoms, nrnb, wcycle, fr,
        cpt_period, max_hours,
        Flags,
        &nsteps_done, at);
    if (inputrec->ePull != epullNO)
      finish_pull(fplog, inputrec->pull);
  } else {
    /* do PME only */
    gmx_pmeonly(*pmedata, cr, nrnb, wcycle, ewaldcoeff, FALSE);
  }

  /* Some timing stats */
  if (MASTER(cr)) {
    realtime = difftime(time(NULL), start_t);
    if ((nodetime = node_time()) == 0)
      nodetime = realtime;
  }
  else
    realtime = 0;

  wallcycle_stop(wcycle, ewcRUN);

  /* Finish up, write some stuff
   * if rerunMD, don't write last frame again
   */
  finish_run(fplog, cr, ftp2fn(efSTO, nfile, fnm),
             inputrec, nrnb, wcycle, nodetime, realtime, nsteps_done,
             EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

  /* Does what it says */
  print_date_and_time(fplog, cr->nodeid,"Finished mdrun");

  /* Close logfile already here if we were appending to it */
  if (MASTER(cr) && (Flags & MD_APPENDFILES))
    gmx_log_close(fplog);

  if (at != NULL)
    at_close(at); /* release memory */
}
#include "macros.h"
#include "main.h"
#include "futil.h"

int main(int argc, char *argv[])
{
  static const char *desc[] = {"self-contained GROMACS"};
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efXTC, "-x",      NULL,       ffOPTWR },
    { efCPT, "-cpi",    NULL,       ffOPTRD },
    { efCPT, "-cpo",    NULL,       ffOPTWR },
    { efSTO, "-c",      "confout",  ffWRITE },
    { efENX, "-e",      "ener",     ffWRITE },
    { efLOG, "-g",      "md",       ffWRITE },
    { efXVG, "-dgdl",   "dgdl",     ffOPTWR },
    { efXVG, "-field",  "field",    ffOPTWR },
    { efXVG, "-table",  "table",    ffOPTRD },
    { efXVG, "-tablep", "tablep",   ffOPTRD },
    { efXVG, "-tableb", "table",    ffOPTRD },
    { efTRX, "-rerun",  "rerun",    ffOPTRD },
    { efXVG, "-tpi",    "tpi",      ffOPTWR },
    { efXVG, "-tpid",   "tpidist",  ffOPTWR },
    { efEDI, "-ei",     "sam",      ffOPTRD },
    { efEDO, "-eo",     "sam",      ffOPTWR },
    { efGCT, "-j",      "wham",     ffOPTRD },
    { efGCT, "-jo",     "bam",      ffOPTWR },
    { efXVG, "-ffout",  "gct",      ffOPTWR },
    { efXVG, "-devout", "deviatie", ffOPTWR },
    { efXVG, "-runav",  "runaver",  ffOPTWR },
    { efMTX, "-mtx",    "nm",       ffOPTWR },
    { efNDX, "-dn",     "dipole",   ffOPTWR },
    { efMDP, "-at",     NULL,       ffOPTRD }  /* at.cfg */
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bPartDec     = FALSE;
  static bool bDDBondCheck = TRUE;
  static bool bDDBondComm  = TRUE;
  static bool bSumEner     = TRUE;
  static bool bVerbose     = FALSE;
  static bool bCompact     = TRUE;
  static bool bSepPot      = FALSE;

  static bool bRerunVSite  = FALSE;
  static bool bConfout     = TRUE;
  static bool bReproducible = FALSE;

  static int  npme = -1;
  static int  nstepout = 100;

  static rvec realddxyz = { 0, 0, 0 };
  static const char *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
  static const char *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
  static real rdd = 0.0, rconstr = 0.0, dlb_scale = 0.8, pforce = -1;
  static const char *ddcsx = NULL, *ddcsy = NULL, *ddcsz = NULL;
  static real cpt_period = 15.0, max_hours = -1;
  static bool bAppendFiles = FALSE, bAddPart = TRUE;
  static int  atpremode = 0;
  static t_pargs pa[] = {
    { "-pd",      FALSE, etBOOL, {&bPartDec },
      "Use particle decompostion" },
    { "-dd",      FALSE, etRVEC, {&realddxyz },
      "Domain decomposition grid, 0 is optimize" },
    { "-npme",    FALSE, etINT, {&npme },
      "Number of separate nodes to be used for PME, -1 is guess" },
    { "-ddorder", FALSE, etENUM, { ddno_opt },
      "DD node order" },
    { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck },
      "Check for all bonded interactions with DD" },
    { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm },
      "HIDDENUse special bonded atom communication when -rdd > cut-off" },
    { "-rdd",     FALSE, etREAL, {&rdd },
      "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
    { "-rcon",    FALSE, etREAL, {&rconstr },
      "Maximum distance for P-LINCS (nm), 0 is estimate" },
    { "-dlb",     FALSE, etENUM, { dddlb_opt },
      "Dynamic load balancing (with DD)" },
    { "-dds",     FALSE, etREAL, {&dlb_scale },
      "Minimum allowed dlb scaling of the DD cell size" },
    { "-ddcsx",   FALSE, etSTR, {&ddcsx },
      "HIDDENThe DD cell sizes in x" },
    { "-ddcsy",   FALSE, etSTR, {&ddcsy },
      "HIDDENThe DD cell sizes in y" },
    { "-ddcsz",   FALSE, etSTR, {&ddcsz },
      "HIDDENThe DD cell sizes in z" },
    { "-sum",     FALSE, etBOOL, {&bSumEner },
      "Sum the energies at every step" },
    { "-v",       FALSE, etBOOL, {&bVerbose },
      "Be loud and noisy" },
    { "-compact", FALSE, etBOOL, {&bCompact },
      "Write a compact log file" },
    { "-seppot",  FALSE, etBOOL, {&bSepPot },
      "Write separate V and dVdl terms for each interaction type and node to the log file(s)" },
    { "-pforce",  FALSE, etREAL, {&pforce },
      "Print all forces larger than this (kJ/mol nm)" },
    { "-reprod",  FALSE, etBOOL, {&bReproducible },
      "Try to avoid optimizations that affect binary reproducibility" },
    { "-cpt",     FALSE, etREAL, {&cpt_period },
      "Checkpoint interval (minutes)" },
    { "-append",  FALSE, etBOOL, {&bAppendFiles },
      "Append to previous output files when continuing from checkpoint" },
    { "-addpart",  FALSE, etBOOL, {&bAddPart },
      "Add the simulation part number to all output files when continuing from checkpoint" },
    { "-maxh",   FALSE, etREAL, {&max_hours },
      "Terminate after 0.99 times this time (hours)" },
    { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite },
      "HIDDENRecalculate virtual site coordinates with -rerun" },
    { "-confout", FALSE, etBOOL, {&bConfout },
      "HIDDENWrite the last configuration with -c" },
    { "-stepout", FALSE, etINT, {&nstepout },
      "HIDDENFrequency of writing the remaining runtime" },
    { "-atpre", FALSE, etINT, {&atpremode},
      "a preparation run" }
  };
  unsigned long Flags, PCA_Flags;
  ivec     ddxyz;
  int      dd_node_order;
  FILE     *fplog;
  int      sim_part;
  char     suffix[STRLEN];

  /* say hello to global variables */
  CALL_UNUSED_MDRUN_H();
  CALL_UNUSED_STATUTIL_H();

  cr = init_par(&argc, &argv);

  PCA_Flags = (PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_CAN_SET_DEFFNM
               | (MASTER(cr) ? 0 : PCA_QUIET));

  parse_common_args(&argc, argv, PCA_Flags,
                    NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL);

  dd_node_order = nenum(ddno_opt);
  cr->npmenodes = npme;

  /* Check if there is ANY checkpoint file available */
  sim_part = 1;
  if (opt2bSet("-cpi", NFILE, fnm))
  {
          sim_part = read_checkpoint_simulation_part(opt2fn("-cpi", NFILE, fnm)) + 1;
          /* sim_part will now be 1 if no checkpoint file was found */
          if (sim_part == 1 && MASTER(cr))
          {
                  fprintf(stdout,"No previous checkpoint file present, assuming this is a new run.\n");
          }
  }

  if (sim_part <= 1)
  {
          bAppendFiles = FALSE;
  }

  if (!bAppendFiles && bAddPart && sim_part > 1)
  {
          /* This is a continuation run, rename trajectory output files (except checkpoint files) */
          /* create new part name first (zero-filled) */
          if (sim_part < 10)
                  sprintf(suffix,"part000%d", sim_part);
          else if (sim_part < 100)
                  sprintf(suffix,"part00%d", sim_part);
          else if (sim_part < 1000)
                  sprintf(suffix,"part0%d", sim_part);
          else
                  sprintf(suffix,"part%d", sim_part);

          add_suffix_to_output_names(fnm, NFILE, suffix);
          fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed %s.\n", sim_part-1, suffix);
  }

  Flags = opt2bSet("-rerun", NFILE, fnm) ? MD_RERUN : 0;
  Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
  Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
  Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
  Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
  Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
  Flags = Flags | (!bSumEner     ? MD_NOGSTAT      : 0);
  Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
  Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
  Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0);
  Flags = Flags | (sim_part > 1    ? MD_STARTFROMCPT : 0);


  /* We postpone opening the log file if we are appending, so we can first truncate
   * the old log file and append to the correct position there instead.
   */
  if ((MASTER(cr) || bSepPot) && !bAppendFiles)
    fplog = gmx_log_open(ftp2fn(efLOG, NFILE, fnm), cr, !bSepPot, Flags);
  else
    fplog = NULL;

  ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
  ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
  ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);

  runner(fplog, cr, NFILE, fnm, bVerbose, bCompact,
           ddxyz, dd_node_order, rdd, rconstr,
           dddlb_opt[0], dlb_scale, ddcsx, ddcsy, ddcsz,
           nstepout, pforce,
           cpt_period, max_hours, Flags, atpremode);

  if (gmx_parallel_env)
    gmx_finalize(cr);

  /* Log file has to be closed in runner if we are appending to it (fplog not set here) */
  if (MASTER(cr) && !bAppendFiles)
    gmx_log_close(fplog);

  return 0;
}


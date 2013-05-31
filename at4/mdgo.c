#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
#endif

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
#include "sighandler.h"
#include "string2.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#define GMXVERSION 40505
#include "mdgoutil.h"


double md(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose,gmx_bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite,gmx_constr_t constr,
             int stepout,t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb,gmx_wallcycle_t wcycle,
             gmx_edsam_t ed,t_forcerec *fr,
             int repl_ex_nst,int repl_ex_seed,
             real cpt_period,real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime, ago_t *ago)
{
  gmx_mdoutf_t *outf;
  gmx_large_int_t step,step_rel;
  double     run_time;
  double     t,t0,lam0;
  gmx_bool       bGStatEveryStep,bGStat,bNstEner,bCalcEnerPres;
  gmx_bool       bNS,bNStList,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bStateFromTPX,bInitStep,bLastStep,
             bBornRadii,bStartingFromCpt;
  gmx_bool       bDoDHDL=FALSE;
  gmx_bool       do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
             bForceUpdate=FALSE,bCPT;
  int        mdof_flags;
  gmx_bool       bMasterState;
  int        force_flags,cglo_flags;
  tensor     force_vir,shake_vir,total_vir,tmp_vir,pres;
  int        i,m;
  t_trxstatus *status;
  rvec       mu_tot;
  t_vcm      *vcm;
  t_state    *bufstate=NULL;
  matrix     *scale_tot,pcoupl_mu,M;
  gmx_nlheur_t nlh;
  t_trxframe rerun_fr;
  int        nchkpt=1;

  gmx_localtop_t *top;
  t_mdebin *mdebin=NULL;
  t_state    *state=NULL;
  rvec       *f_global=NULL;
  int        n_xtc=-1;
  rvec       *x_xtc=NULL;
  gmx_enerdata_t *enerd;
  rvec       *f=NULL;
  gmx_global_stat_t gstat;
  gmx_update_t upd=NULL;
  t_graph    *graph=NULL;
  globsig_t   gs;

  gmx_groups_t *groups;
  gmx_ekindata_t *ekind, *ekind_save;
  gmx_shellfc_t shellfc;
  int         count,nconverged=0;
  real        timestep=0;
  double      tcount=0;
  gmx_bool        bConverged=TRUE,bOK,bSumEkinhOld;
  gmx_bool        bResetCountersHalfMaxH=FALSE;
  gmx_bool        bVV,bIterations,bFirstIterate,bTemp,bPres,bTrotter;
  real        dvdl;
  int         a0,a1;
  rvec        *cbuf=NULL;
  matrix      lastbox;
  real        veta_save,pcurr,scalevir,tracevir;
  real        vetanew = 0;
  double      cycles;
  real        saved_conserved_quantity = 0;
  real        last_ekin = 0;
  t_extmass   MassQ;
  int         **trotter_seq;
  char        sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
  int         handled_stop_condition=gmx_stop_cond_none; /* compare to get_stop_condition*/
  gmx_iterate_t iterate;
  gmx_large_int_t multisim_nsteps=-1; /* number of steps to do  before first multisim
                                        simulation stops. If equal to zero, don't
                                        communicate any more between multisims.*/
#ifdef GMX_FAHCORE
  /* Temporary addition for FAHCORE checkpointing */
  int chkpt_ret;
#endif
  (void) deviceOptions; (void) repl_ex_seed; (void) ed;
  /* Check for special mdrun options */
  bRerunMD = (Flags & MD_RERUN);
  if (Flags & MD_RESETCOUNTERSHALFWAY)
  {
    if (ir->nsteps > 0)
    {
      /* Signal to reset the counters half the simulation steps. */
      wcycle_set_reset_counters(wcycle,ir->nsteps/2);
    }
    /* Signal to reset the counters halfway the simulation time. */
    bResetCountersHalfMaxH = (max_hours > 0);
  }

  /* md-vv uses averaged full step velocities for T-control
     md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
     md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
  bVV = EI_VV(ir->eI);
  if (bVV) /* to store the initial velocities while computing virial */
  {
    snew(cbuf,top_global->natoms);
  }
  /* all the iteratative cases - only if there are constraints */
  bIterations = ((IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
  bTrotter = (bVV && (IR_NPT_TROTTER(ir) || (IR_NVT_TROTTER(ir))));

  if (bRerunMD)
  {
    /* Since we don't know if the frames read are related in any way,
     * rebuild the neighborlist at every step.
     */
    ir->nstlist       = 1;
    ir->nstcalcenergy = 1;
    nstglobalcomm     = 1;
  }

  check_ir_old_tpx_versions(cr,fplog,ir,top_global);

  nstglobalcomm = check_nstglobalcomm(fplog,cr,nstglobalcomm,ir);
  bGStatEveryStep = (nstglobalcomm == 1);

  if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
  {
    fprintf(fplog,
            "To reduce the energy communication with nstlist = -1\n"
            "the neighbor list validity should not be checked at every step,\n"
            "this means that exact integration is not guaranteed.\n"
            "The neighbor list validity is checked after:\n"
            "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
            "In most cases this will result in exact integration.\n"
            "This reduces the energy communication by a factor of 2 to 3.\n"
            "If you want less energy communication, set nstlist > 3.\n\n");
  }

  if (bRerunMD)
  {
    ir->nstxtcout = 0;
  }
  groups = &top_global->groups;

  /* Initial values */
  init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
          nrnb,top_global,&upd,
          nfile,fnm,&outf,&mdebin,
          force_vir,shake_vir,mu_tot,&bSimAnn,&vcm,state_global,Flags);

  clear_mat(total_vir);
  clear_mat(pres);
  /* Energy terms and groups */
  snew(enerd,1);
  init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
  if (DOMAINDECOMP(cr))
  {
    f = NULL;
  }
  else
  {
    snew(f,top_global->natoms);
  }

  /* Kinetic energy data */
  snew(ekind,1);
  init_ekindata(fplog,top_global,&(ir->opts),ekind);
  /* needed for iteration of constraints */
  snew(ekind_save,1);
  init_ekindata(fplog,top_global,&(ir->opts),ekind_save);
  /* Copy the cos acceleration to the groups struct */
  ekind->cosacc.cos_accel = ir->cos_accel;

  gstat = global_stat_init(ir);

  /* Check for polarizable models and flexible constraints */
  shellfc = init_shell_flexcon(fplog,
                               top_global,n_flexible_constraints(constr),
                               (ir->bContinuation ||
                                (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                               NULL : state_global->x);

  if (DEFORM(*ir))
  {
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
    set_deform_reference_box(upd,
                             deform_init_init_step_tpx,
                             deform_init_box_tpx);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
  }

  {
    double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
    if ((io > 2000) && MASTER(cr))
        fprintf(stderr,
                "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                io);
  }

  if (DOMAINDECOMP(cr)) {
    top = dd_init_local_top(top_global);

    snew(state,1);
    dd_init_local_state(cr->dd,state_global,state);

    if (DDMASTER(cr->dd) && ir->nstfout) {
      snew(f_global,state_global->natoms);
    }
  } else {
    if (PAR(cr)) {
      /* Initialize the particle decomposition and split the topology */
      top = split_system(fplog,top_global,ir,cr);

      pd_cg_range(cr,&fr->cg0,&fr->hcg);
      pd_at_range(cr,&a0,&a1);
    } else {
      top = gmx_mtop_generate_local_top(top_global,ir);

      a0 = 0;
      a1 = top_global->natoms;
    }

    state = partdec_init_local_state(cr,state_global);
    f_global = f;

    atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

    if (vsite) {
      set_vsite_top(vsite,top,mdatoms,cr);
    }

    if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
      graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
    }

    if (shellfc) {
      make_local_shells(cr,mdatoms,shellfc);
    }

    if (ir->pull && PAR(cr)) {
      dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
    }
  }

  if (DOMAINDECOMP(cr))
  {
    /* Distribute the charge groups over the nodes from the master node */
    dd_partition_system(fplog,ir->init_step,cr,TRUE,1,
                        state_global,top_global,ir,
                        state,&f,mdatoms,top,fr,
                        vsite,shellfc,constr,
                        nrnb,wcycle,FALSE);
  }

  update_mdatoms(mdatoms,state->lambda);

  if (MASTER(cr))
  {
    if (opt2bSet("-cpi",nfile,fnm))
    {
      /* Update mdebin with energy history if appending to output files */
      if ( Flags & MD_APPENDFILES )
      {
        restore_energyhistory_from_state(mdebin,&state_global->enerhist);
      }
      else
      {
        /* We might have read an energy history from checkpoint,
         * free the allocated memory and reset the counts.
         */
        done_energyhistory(&state_global->enerhist);
        init_energyhistory(&state_global->enerhist);
      }
    }
    /* Set the initial energy history in state by updating once */
    update_energyhistory(&state_global->enerhist,mdebin);
  }

  if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG)) {
    /* Set the random state if we read a checkpoint file */
    set_stochd_state(upd,state);
  }

  /* Initialize constraints */
  if (constr) {
    if (!DOMAINDECOMP(cr))
        set_constraints(constr,top,ir,mdatoms,cr);
  }

  if (repl_ex_nst > 0)
  {
    /* We need to be sure replica exchange can only occur
     * when the energies are current */
    check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                    "repl_ex_nst",&repl_ex_nst);
    /* This check needs to happen before inter-simulation
     * signals are initialized, too */
  }
  if (!ir->bContinuation && !bRerunMD)
  {
    if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
    {
      /* Set the velocities of frozen particles to zero */
      for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
      {
        for(m=0; m<DIM; m++)
        {
          if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
          {
            state->v[i][m] = 0;
          }
        }
      }
    }

    if (constr)
    {
      /* Constrain the initial coordinates and velocities */
      do_constrain_first(fplog,constr,ir,mdatoms,state,f,
                         graph,cr,nrnb,fr,top,shake_vir);
    }
    if (vsite)
    {
      /* Construct the virtual sites for the initial configuration */
      construct_vsites(fplog,vsite,state->x,nrnb,(real) ir->delta_t,NULL,
                       top->idef.iparams,top->idef.il,
                       fr->ePBC,fr->bMolPBC,graph,cr,state->box);
    }
  }


  /* I'm assuming we need global communication the first time! MRS */
  cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                | (bVV ? CGLO_PRESSURE:0)
                | (bVV ? CGLO_CONSTRAINT:0)
                | (bRerunMD ? CGLO_RERUNMD:0)
                | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN:0));

  bSumEkinhOld = FALSE;
  compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                  wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                  constr,NULL,FALSE,state->box,
                  top_global,&pcurr,top_global->natoms,&bSumEkinhOld,cglo_flags);
  if (ir->eI == eiVVAK) {
    /* a second call to get the half step temperature initialized as well */
    /* we do the same call as above, but turn the pressure off -- internally to
       compute_globals, this is recognized as a velocity verlet half-step
       kinetic energy calculation.  This minimized excess variables, but
       perhaps loses some logic?*/

    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                    constr,NULL,FALSE,state->box,
                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                    cglo_flags &~ CGLO_PRESSURE);
  }

  /* Calculate the initial half step temperature, and save the ekinh_old */
  if (!(Flags & MD_STARTFROMCPT))
  {
    for(i=0; (i<ir->opts.ngtc); i++)
    {
      copy_mat(ekind->tcstat[i].ekinh,ekind->tcstat[i].ekinh_old);
    }
  }
  if (ir->eI != eiVV)
  {
    enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                 and there is no previous step */
  }

  /* if using an iterative algorithm, we need to create a working directory for the state. */
  if (bIterations)
  {
    bufstate = init_bufstate(state);
  }
  /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
     temperature control */
  trotter_seq = init_npt_vars(ir,state,&MassQ,bTrotter);

  if (MASTER(cr))
  {
    if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
    {
      fprintf(fplog,
              "RMS relative constraint deviation after constraining: %.2e\n",
              constr_rmsd(constr,FALSE));
    }
    fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
    if (bRerunMD)
    {
      fprintf(stderr,"starting md rerun '%s', reading coordinates from"
              " input trajectory '%s'\n\n",
              *(top_global->name),opt2fn("-rerun",nfile,fnm));
      if (bVerbose)
      {
        fprintf(stderr,"Calculated time to finish depends on nsteps from "
                "run input file,\nwhich may not correspond to the time "
                "needed to process input trajectory.\n\n");
      }
    }
    else
    {
      char tbuf[20];
      fprintf(stderr,"starting mdrun '%s'\n",
              *(top_global->name));
      if (ir->nsteps >= 0)
      {
        sprintf(tbuf,"%8.1f",(ir->init_step+ir->nsteps)*ir->delta_t);
      }
      else
      {
        sprintf(tbuf,"%s","infinite");
      }
      if (ir->init_step > 0)
      {
        fprintf(stderr,"%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                gmx_step_str(ir->init_step+ir->nsteps,sbuf),tbuf,
                gmx_step_str(ir->init_step,sbuf2),
                ir->init_step*ir->delta_t);
      }
      else
      {
        fprintf(stderr,"%s steps, %s ps.\n",
                gmx_step_str(ir->nsteps,sbuf),tbuf);
      }
    }
    fprintf(fplog,"\n");
  }

  /* Set and write start time */
  runtime_start(runtime);
  print_date_and_time(fplog,cr->nodeid,"Started mdrun",runtime);
  wallcycle_start(wcycle,ewcRUN);
  if (fplog)
      fprintf(fplog,"\n");

  /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
  chkpt_ret=fcCheckPointParallel( cr->nodeid,
                                  NULL,0);
  if ( chkpt_ret == 0 )
      gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", 0 );
#endif

  /* must be called after atoms2md() assigns mdatoms->start and mdatoms->homenr */
  agox_assign_mkls(ago, mdatoms->start, mdatoms->start+mdatoms->homenr, cr);

  /***********************************************************
   *
   *             Loop over MD steps
   *
   ************************************************************/

  /* if rerunMD then read coordinates and velocities from input trajectory */
  if (bRerunMD)
  {
    if (getenv("GMX_FORCE_UPDATE"))
    {
      bForceUpdate = TRUE;
    }

    rerun_fr.natoms = 0;
    if (MASTER(cr))
    {
      bNotLastFrame = read_first_frame(oenv,&status,
                                       opt2fn("-rerun",nfile,fnm),
                                       &rerun_fr,TRX_NEED_X | TRX_READ_V);
      if (rerun_fr.natoms != top_global->natoms)
      {
        gmx_fatal(FARGS,
                  "Number of atoms in trajectory (%d) does not match the "
                  "run input file (%d)\n",
                  rerun_fr.natoms,top_global->natoms);
      }
      if (ir->ePBC != epbcNONE)
      {
        if (!rerun_fr.bBox)
        {
          gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
        }
        if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
        {
          gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
        }
      }
    }

    if (PAR(cr))
    {
      rerun_parallel_comm(cr,&rerun_fr,&bNotLastFrame);
    }

    if (ir->ePBC != epbcNONE)
    {
      /* Set the shift vectors.
       * Necessary here when have a static box different from the tpr box.
       */
      calc_shifts(rerun_fr.box,fr->shift_vec);
    }
  }

  /* loop over MD steps or if rerunMD to end of input trajectory */
  bFirstStep = TRUE;
  /* Skip the first Nose-Hoover integration when we get the state from tpx */
  bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
  bInitStep = bFirstStep && (bStateFromTPX || bVV);
  bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
  bLastStep    = FALSE;
  bSumEkinhOld = FALSE;

  init_global_signals(&gs,cr,ir,repl_ex_nst);

  step = ir->init_step;
  step_rel = 0;

  if (ir->nstlist == -1)
  {
    init_nlistheuristics(&nlh,bGStatEveryStep,step);
  }

  if (MULTISIM(cr) && (repl_ex_nst <=0 ))
  {
    /* check how many steps are left in other sims */
    multisim_nsteps=get_multisim_nsteps(cr, ir->nsteps);
  }


  /* and stop now if we should */
  bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps) ||
               ((multisim_nsteps >= 0) && (step_rel >= multisim_nsteps )));
  while (!bLastStep || (bRerunMD && bNotLastFrame)) {

    wallcycle_start(wcycle,ewcSTEP);

    GMX_MPE_LOG(ev_timestep1);

    if (bRerunMD) {
      if (rerun_fr.bStep) {
        step = rerun_fr.step;
        step_rel = step - ir->init_step;
      }
      if (rerun_fr.bTime) {
        t = rerun_fr.time;
      }
      else
      {
        t = (double) step;
      }
    }
    else
    {
      bLastStep = (step_rel == ir->nsteps);
      t = t0 + step*ir->delta_t;
    }

    if (ir->efep != efepNO)
    {
      if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
      {
        state_global->lambda = rerun_fr.lambda;
      }
      else
      {
        state_global->lambda = (real)( lam0 + step*ir->delta_lambda );
      }
      state->lambda = state_global->lambda;
      bDoDHDL = do_per_step(step,ir->nstdhdl);
    }

    if (bSimAnn)
    {
      update_annealing_target_temp(&(ir->opts), (real) t);
    }

    if (bRerunMD)
    {
      if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
      {
        for(i=0; i<state_global->natoms; i++)
        {
          copy_rvec(rerun_fr.x[i],state_global->x[i]);
        }
        if (rerun_fr.bV)
        {
          for(i=0; i<state_global->natoms; i++)
          {
            copy_rvec(rerun_fr.v[i],state_global->v[i]);
          }
        }
        else
        {
          for(i=0; i<state_global->natoms; i++)
          {
            clear_rvec(state_global->v[i]);
          }
          if (bRerunWarnNoV)
          {
            fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                    "         Ekin, temperature and pressure are incorrect,\n"
                    "         the virial will be incorrect when constraints are present.\n"
                    "\n");
            bRerunWarnNoV = FALSE;
          }
        }
      }
      copy_mat(rerun_fr.box,state_global->box);
      copy_mat(state_global->box,state->box);

      if (vsite && (Flags & MD_RERUN_VSITE))
      {
        if (DOMAINDECOMP(cr))
        {
          gmx_fatal(FARGS,"Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
        }
        if (graph)
        {
          /* Following is necessary because the graph may get out of sync
           * with the coordinates if we only have every N'th coordinate set
           */
          mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
          shift_self(graph,state->box,state->x);
        }
        construct_vsites(fplog,vsite,state->x,nrnb,(real) ir->delta_t,state->v,
                         top->idef.iparams,top->idef.il,
                         fr->ePBC,fr->bMolPBC,graph,cr,state->box);
        if (graph)
        {
          unshift_self(graph,state->box,state->x);
        }
      }
    }

    /* Stop Center of Mass motion */
    bStopCM = (ir->comm_mode != ecmNO && do_per_step(step,ir->nstcomm));

    if (bRerunMD)
    {
      /* for rerun MD always do Neighbour Searching */
      bNS = (bFirstStep || ir->nstlist != 0);
      bNStList = bNS;
    }
    else
    {
      /* Determine whether or not to do Neighbour Searching and LR */
      bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

      bNS = (bFirstStep || bNStList ||
             (ir->nstlist == -1 && nlh.nabnsb > 0));

      if (bNS && ir->nstlist == -1)
      {
        set_nlistheuristics(&nlh,bFirstStep,step);
      }
    }

    /* check whether we should stop because another simulation has
       stopped. */
    if (MULTISIM(cr))
    {
      if ( (multisim_nsteps >= 0) &&  (step_rel >= multisim_nsteps)  &&
           (multisim_nsteps != ir->nsteps) )
      {
        if (bNS)
        {
          if (MASTER(cr))
          {
            fprintf(stderr,
                    "Stopping simulation %d because another one has finished\n",
                    cr->ms->sim);
          }
          bLastStep=TRUE;
          gs.sig[eglsCHKPT] = 1;
        }
      }
    }

    /* < 0 means stop at next step, > 0 means stop at next NS step */
    if ( (gs.set[eglsSTOPCOND] < 0 ) ||
         ( (gs.set[eglsSTOPCOND] > 0 ) && ( bNS || ir->nstlist==0)) )
    {
      bLastStep = TRUE;
    }

    /* Determine whether or not to update the Born radii if doing GB */
    bBornRadii=bFirstStep;
    if (ir->implicit_solvent && (step % ir->nstgbradii==0))
    {
      bBornRadii=TRUE;
    }

    do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
    do_verbose = bVerbose &&
              (step % stepout == 0 || bFirstStep || bLastStep);

    if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
    {
      if (bRerunMD)
      {
        bMasterState = TRUE;
      }
      else
      {
        bMasterState = FALSE;
        /* Correct the new box if it is too skewed */
        if (DYNAMIC_BOX(*ir))
        {
          if (correct_box(fplog,step,state->box,graph))
          {
            bMasterState = TRUE;
          }
        }
        if (DOMAINDECOMP(cr) && bMasterState)
        {
          dd_collect_state(cr->dd,state,state_global);
        }
      }

      if (DOMAINDECOMP(cr))
      {
        /* Repartition the domain decomposition */
        wallcycle_start(wcycle,ewcDOMDEC);
        dd_partition_system(fplog,step,cr,
                            bMasterState,nstglobalcomm,
                            state_global,top_global,ir,
                            state,&f,mdatoms,top,fr,
                            vsite,shellfc,constr,
                            nrnb,wcycle,do_verbose);
        wallcycle_stop(wcycle,ewcDOMDEC);
        /* If using an iterative integrator, reallocate space to match the decomposition */
      }
    }

    if (MASTER(cr) && do_log)
    {
      print_ebin_header(fplog,step,t,state->lambda);
    }

    if (ir->efep != efepNO)
    {
      update_mdatoms(mdatoms,state->lambda);
    }

    if (bRerunMD && rerun_fr.bV)
    {

      /* We need the kinetic energy at minus the half step for determining
       * the full step kinetic energy and possibly for T-coupling.*/
      /* This may not be quite working correctly yet . . . . */
      compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                      wcycle,enerd,NULL,NULL,NULL,NULL,mu_tot,
                      constr,NULL,FALSE,state->box,
                      top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                      CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
    }
    clear_mat(force_vir);

    GMX_MPE_LOG(ev_timestep2);

    /* We write a checkpoint at this MD step when:
     * either at an NS step when we signalled through gs,
     * or at the last step (but not when we do not want confout),
     * but never at the first step or with rerun.
     */
    bCPT = (((gs.set[eglsCHKPT] && (bNS || ir->nstlist == 0)) ||
             (bLastStep && (Flags & MD_CONFOUT))) &&
            step > ir->init_step && !bRerunMD);
    if (bCPT)
    {
      gs.set[eglsCHKPT] = 0;
    }

    /* Determine the energy and pressure:
     * at nstcalcenergy steps and at energy output steps (set below).
     */
    bNstEner = do_per_step(step,ir->nstcalcenergy);
    bCalcEnerPres =
        (bNstEner ||
         (ir->epc != epcNO && do_per_step(step,ir->nstpcouple)));

    /* Do we need global communication ? */
    bGStat = (bCalcEnerPres || bStopCM ||
              do_per_step(step,nstglobalcomm) ||
              (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

    do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

    if (do_ene || do_log)
    {
      bCalcEnerPres = TRUE;
      bGStat        = TRUE;
    }

    /* these CGLO_ options remain the same throughout the iteration */
    cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                  (bStopCM ? CGLO_STOPCM : 0) |
                  (bGStat ? CGLO_GSTAT : 0)
        );

    force_flags = (GMX_FORCE_STATECHANGED |
                   ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                   GMX_FORCE_ALLFORCES |
                   (bNStList ? GMX_FORCE_DOLR : 0) |
                   GMX_FORCE_SEPLRF |
                   (bCalcEnerPres ? GMX_FORCE_VIRIAL : 0) |
                   (bDoDHDL ? GMX_FORCE_DHDL : 0)
        );

    if (shellfc)
    {
      /* Now is the time to relax the shells */
      count=relax_shell_flexcon(fplog,cr,bVerbose, step,
                                ir,bNS,force_flags,
                                bStopCM,top,top_global,
                                constr,enerd,fcd,
                                state,f,force_vir,mdatoms,
                                nrnb,wcycle,graph,groups,
                                shellfc,fr,bBornRadii,t,mu_tot,
                                state->natoms,&bConverged,vsite,
                                outf->fp_field);
      tcount+=count;

      if (bConverged)
      {
        nconverged++;
      }
    }
    else
    {
      /* The coordinates (x) are shifted (to get whole molecules)
       * in do_force.
       * This is parallellized as well, and does communication too.
       * Check comments in sim_util.c
       */

      agox_doforce(fplog,cr,ir,step,nrnb,wcycle,top,top_global,groups,
               state->box,state->x,&state->hist,
               f,force_vir,mdatoms,enerd,fcd,
               state->lambda,graph,
               fr,vsite,mu_tot,t,outf->fp_field,NULL,bBornRadii,
               (bNS ? GMX_FORCE_NS : 0) | force_flags, ago);
    }

    GMX_BARRIER(cr->mpi_comm_mygroup);

    if (bVV && !bStartingFromCpt && !bRerunMD)
    /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
    {
      if (ir->eI==eiVV && bInitStep)
      {
        /* if using velocity verlet with full time step Ekin,
         * take the first half step only to compute the
         * virial for the first step. From there,
         * revert back to the initial coordinates
         * so that the input is actually the initial step.
         */
        copy_rvecn(state->v,cbuf,0,state->natoms); /* should make this better for parallelizing? */
      } else {
        /* this is for NHC in the Ekin(t+dt/2) version of vv */
        trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ1);
      }

      update_coords(fplog,step,ir,mdatoms,state,
                    f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                    ekind,M,wcycle,upd,bInitStep,etrtVELOCITY1,
                    cr,nrnb,constr,&top->idef);

      if (bIterations)
      {
        gmx_iterate_init(&iterate,bIterations && !bInitStep);
      }
      /* for iterations, we save these vectors, as we will be self-consistently iterating
         the calculations */

      /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */

      /* save the state */
      if (bIterations && iterate.bIterate) {
        copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
      }

      bFirstIterate = TRUE;
      while (bFirstIterate || (bIterations && iterate.bIterate))
      {
        if (bIterations && iterate.bIterate)
        {
          copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
          if (bFirstIterate && bTrotter)
          {
            /* The first time through, we need a decent first estimate
               of veta(t+dt) to compute the constraints.  Do
               this by computing the box volume part of the
               trotter integration at this time. Nothing else
               should be changed by this routine here.  If
               !(first time), we start with the previous value
               of veta.  */

            veta_save = state->veta;
            trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ0);
            vetanew = state->veta;
            state->veta = veta_save;
          }
        }

        bOK = TRUE;
        if ( !bRerunMD || rerun_fr.bV || bForceUpdate) {  /* Why is rerun_fr.bV here?  Unclear. */
          dvdl = 0;

          update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                             &top->idef,shake_vir,NULL,
                             cr,nrnb,wcycle,upd,constr,
                             bInitStep,TRUE,bCalcEnerPres,vetanew);

          if (!bOK)
          {
            gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
          }

        }
        else if (graph)
        { /* Need to unshift here if a do_force has been
          called in the previous step */
         unshift_self(graph,state->box,state->x);
        }


        /* if VV, compute the pressure and constraints */
        /* For VV2, we strictly only need this if using pressure
         * control, but we really would like to have accurate pressures
         * printed out.
         * Think about ways around this in the future?
         * For now, keep this choice in comments.
         */
        /*bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir)); */
            /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));*/
        bPres = TRUE;
        bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK));
        compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                        wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                        constr,NULL,FALSE,state->box,
                        top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                        cglo_flags
                        | CGLO_ENERGY
                        | (bTemp ? CGLO_TEMPERATURE:0)
                        | (bPres ? CGLO_PRESSURE : 0)
                        | (bPres ? CGLO_CONSTRAINT : 0)
                        | ((bIterations && iterate.bIterate) ? CGLO_ITERATE : 0)
                        | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                        | CGLO_SCALEEKIN
            );
        /* explanation of above:
           a) We compute Ekin at the full time step
           if 1) we are using the AveVel Ekin, and it's not the
           initial step, or 2) if we are using AveEkin, but need the full
           time step kinetic energy for the pressure (always true now, since we want accurate statistics).
           b) If we are using EkinAveEkin for the kinetic energy for the temperture control, we still feed in
           EkinAveVel because it's needed for the pressure */

        /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
        if (!bInitStep)
        {
          if (bTrotter)
          {
            trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ2);
          }
          else
          {
            update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
          }
        }

        if (bIterations &&
            done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                           state->veta,&vetanew))
        {
          break;
        }
        bFirstIterate = FALSE;
      }

      if (bTrotter && !bInitStep) {
        copy_mat(shake_vir,state->svir_prev);
        copy_mat(force_vir,state->fvir_prev);
        if (IR_NVT_TROTTER(ir) && ir->eI==eiVV) {
          /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
          enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,NULL,(ir->eI==eiVV),FALSE,FALSE);
          enerd->term[F_EKIN] = trace(ekind->ekin);
        }
      }
      /* if it's the initial step, we performed this first step just to get the constraint virial */
      if (bInitStep && ir->eI==eiVV) {
        copy_rvecn(cbuf,state->v,0,state->natoms);
      }

      if (fr->bSepDVDL && fplog && do_log)
      {
        fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
      }
      enerd->term[F_DHDL_CON] += dvdl;

      GMX_MPE_LOG(ev_timestep1);
    }

    /* MRS -- now done iterating -- compute the conserved quantity */
    if (bVV) {
      saved_conserved_quantity = compute_conserved_from_auxiliary(ir,state,&MassQ);
      if (ir->eI==eiVV)
      {
        last_ekin = enerd->term[F_EKIN]; /* does this get preserved through checkpointing? */
      }
      if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
      {
        saved_conserved_quantity -= enerd->term[F_DISPCORR];
      }
    }

    /* ########  END FIRST UPDATE STEP  ############## */
    /* ########  If doing VV, we now have v(dt) ###### */

    /* ################## START TRAJECTORY OUTPUT ################# */

    /* Now we have the energies and forces corresponding to the
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */
    GMX_MPE_LOG(ev_output_start);

    mdof_flags = 0;
    if (do_per_step(step,ir->nstxout)) { mdof_flags |= MDOF_X; }
    if (do_per_step(step,ir->nstvout)) { mdof_flags |= MDOF_V; }
    if (do_per_step(step,ir->nstfout)) { mdof_flags |= MDOF_F; }
    if (do_per_step(step,ir->nstxtcout)) { mdof_flags |= MDOF_XTC; }
    if (bCPT) { mdof_flags |= MDOF_CPT; };

#if defined(GMX_FAHCORE) || defined(GMX_WRITELASTSTEP)
    if (bLastStep)
    {
      /* Enforce writing positions and velocities at end of run */
      mdof_flags |= (MDOF_X | MDOF_V);
    }
#endif
#ifdef GMX_FAHCORE
    if (MASTER(cr))
        fcReportProgress( ir->nsteps, step );

    /* sync bCPT and fc record-keeping */
    if (bCPT && MASTER(cr))
        fcRequestCheckPoint();
#endif

    if (mdof_flags != 0)
    {
      wallcycle_start(wcycle,ewcTRAJ);
      if (bCPT)
      {
        if (state->flags & (1<<estLD_RNG))
        {
          get_stochd_state(upd,state);
        }
        if (MASTER(cr))
        {
          if (bSumEkinhOld)
          {
            state_global->ekinstate.bUpToDate = FALSE;
          }
          else
          {
            update_ekinstate(&state_global->ekinstate,ekind);
            state_global->ekinstate.bUpToDate = TRUE;
          }
          update_energyhistory(&state_global->enerhist,mdebin);
        }
      }
      write_traj(fplog,cr,outf,mdof_flags,top_global,
                 step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
      if (bCPT)
      {
        nchkpt++;
        bCPT = FALSE;
      }
      if (bLastStep && step_rel == ir->nsteps &&
          (Flags & MD_CONFOUT) && MASTER(cr) &&
          !bRerunMD)
      {
        /* x and v have been collected in write_traj,
         * because a checkpoint file will always be written
         * at the last step.
         */
        fprintf(stderr,"\nWriting final coordinates.\n");
        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols &&
            DOMAINDECOMP(cr))
        {
          /* Make molecules whole only for confout writing */
          do_pbc_mtop(fplog,ir->ePBC,state->box,top_global,state_global->x);
        }
        write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                            *top_global->name,top_global,
                            state_global->x,state_global->v,
                            ir->ePBC,state->box);
      }
      wallcycle_stop(wcycle,ewcTRAJ);
    }
    GMX_MPE_LOG(ev_output_finish);

    /* kludge -- virial is lost with restart for NPT control. Must restart */
    if (bStartingFromCpt && bVV)
    {
      copy_mat(state->svir_prev,shake_vir);
      copy_mat(state->fvir_prev,force_vir);
    }
    /*  ################## END TRAJECTORY OUTPUT ################ */

    /* Determine the wallclock run time up till now */
    run_time = gmx_gettime() - (double)runtime->real;

    /* Check whether everything is still allright */
    if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREADS
        && MASTER(cr)
#endif
        )
    {
      /* this is just make gs.sig compatible with the hack
         of sending signals around by MPI_Reduce with together with
         other floats */
      if ( gmx_get_stop_condition() == gmx_stop_cond_next_ns )
          gs.sig[eglsSTOPCOND]=1;
      if ( gmx_get_stop_condition() == gmx_stop_cond_next )
          gs.sig[eglsSTOPCOND]=-1;
      /* < 0 means stop at next step, > 0 means stop at next NS step */
      if (fplog)
      {
        fprintf(fplog,
                "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                gmx_get_signal_name(),
                gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
        fflush(fplog);
      }
      fprintf(stderr,
              "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
              gmx_get_signal_name(),
              gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
      fflush(stderr);
      handled_stop_condition=(int)gmx_get_stop_condition();
    }
    else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
             (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
             gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
    {
      /* Signal to terminate the run */
      gs.sig[eglsSTOPCOND] = 1;
      if (fplog)
      {
        fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
      }
      fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
    }

    if (bResetCountersHalfMaxH && MASTER(cr) &&
        run_time > max_hours*60.0*60.0*0.495)
    {
      gs.sig[eglsRESETCOUNTERS] = 1;
    }

    if (ir->nstlist == -1 && !bRerunMD)
    {
      /* When bGStatEveryStep=FALSE, global_stat is only called
       * when we check the atom displacements, not at NS steps.
       * This means that also the bonded interaction count check is not
       * performed immediately after NS. Therefore a few MD steps could
       * be performed with missing interactions.
       * But wrong energies are never written to file,
       * since energies are only written after global_stat
       * has been called.
       */
      if (step >= nlh.step_nscheck)
      {
        nlh.nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,
                                             nlh.scale_tot,state->x);
      }
      else
      {
        /* This is not necessarily true,
         * but step_nscheck is determined quite conservatively.
         */
        nlh.nabnsb = 0;
      }
    }

    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication,
     *  otherwise the other nodes don't know.
     */
    if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                       cpt_period >= 0 &&
                       (cpt_period == 0 ||
                        run_time >= nchkpt*cpt_period*60.0)) &&
        gs.set[eglsCHKPT] == 0)
    {
      gs.sig[eglsCHKPT] = 1;
    }

    if (bIterations)
    {
      gmx_iterate_init(&iterate,bIterations);
    }

    /* for iterations, we save these vectors, as we will be redoing the calculations */
    if (bIterations && iterate.bIterate)
    {
      copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
    }
    bFirstIterate = TRUE;
    while (bFirstIterate || (bIterations && iterate.bIterate))
    {
      /* We now restore these vectors to redo the calculation with improved extended variables */
      if (bIterations)
      {
        copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
      }

      /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
         so scroll down for that logic */

      /* #########   START SECOND UPDATE STEP ################# */
      GMX_MPE_LOG(ev_update_start);
      /* Box is changed in update() when we do pressure coupling,
       * but we should still use the old box for energy corrections and when
       * writing it to the energy file, so it matches the trajectory files for
       * the same timestep above. Make a copy in a separate array.
       */
      copy_mat(state->box,lastbox);

      bOK = TRUE;
      if (!(bRerunMD && !rerun_fr.bV && !bForceUpdate))
      {
        wallcycle_start(wcycle,ewcUPDATE);
        dvdl = 0;
        /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
        if (bTrotter)
        {
          if (bIterations && iterate.bIterate)
          {
            if (bFirstIterate)
            {
              scalevir = 1;
            }
            else
            {
              /* we use a new value of scalevir to converge the iterations faster */
              scalevir = tracevir/trace(shake_vir);
            }
            msmul(shake_vir,scalevir,shake_vir);
            m_add(force_vir,shake_vir,total_vir);
            clear_mat(shake_vir);
          }
          trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ3);
      /* We can only do Berendsen coupling after we have summed
       * the kinetic energy or virial. Since the happens
       * in global_state after update, we should only do it at
       * step % nstlist = 1 with bGStatEveryStep=FALSE.
       */
        }
        else
        {
          update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
          update_pcouple(fplog,step,ir,state,pcoupl_mu,M,wcycle,
                         upd,bInitStep);
        }

        if (bVV)
        {
          /* velocity half-step update */
          update_coords(fplog,step,ir,mdatoms,state,f,
                        fr->bTwinRange && bNStList,fr->f_twin,fcd,
                        ekind,M,wcycle,upd,FALSE,etrtVELOCITY2,
                        cr,nrnb,constr,&top->idef);
        }

        /* Above, initialize just copies ekinh into ekin,
         * it doesn't copy position (for VV),
         * and entire integrator for MD.
         */

        if (ir->eI==eiVVAK)
        {
          copy_rvecn(state->x,cbuf,0,state->natoms);
        }

        update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                      ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
        wallcycle_stop(wcycle,ewcUPDATE);

        update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                           &top->idef,shake_vir,force_vir,
                           cr,nrnb,wcycle,upd,constr,
                           bInitStep,FALSE,bCalcEnerPres,state->veta);

        if (ir->eI==eiVVAK)
        {
          /* erase F_EKIN and F_TEMP here? */
          /* just compute the kinetic energy at the half step to perform a trotter step */
          compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                          wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                          constr,NULL,FALSE,lastbox,
                          top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                          cglo_flags | CGLO_TEMPERATURE
              );
          wallcycle_start(wcycle,ewcUPDATE);
          trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ4);
          /* now we know the scaling, we can compute the positions again again */
          copy_rvecn(cbuf,state->x,0,state->natoms);

          update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                        ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
          wallcycle_stop(wcycle,ewcUPDATE);

          /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
          /* are the small terms in the shake_vir here due
           * to numerical errors, or are they important
           * physically? I'm thinking they are just errors, but not completely sure.
           * For now, will call without actually constraining, constr=NULL*/
          update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                             &top->idef,tmp_vir,force_vir,
                             cr,nrnb,wcycle,upd,NULL,
                             bInitStep,FALSE,bCalcEnerPres,
                             state->veta);
        }
        if (!bOK)
        {
          gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
        }

        if (fr->bSepDVDL && fplog && do_log)
        {
          fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
        }
        enerd->term[F_DHDL_CON] += dvdl;
      }
      else if (graph)
      {
        /* Need to unshift here */
        unshift_self(graph,state->box,state->x);
      }

      GMX_BARRIER(cr->mpi_comm_mygroup);
      GMX_MPE_LOG(ev_update_finish);

      if (vsite != NULL)
      {
        wallcycle_start(wcycle,ewcVSITECONSTR);
        if (graph != NULL)
        {
          shift_self(graph,state->box,state->x);
        }
        construct_vsites(fplog,vsite,state->x,nrnb,(real) ir->delta_t,state->v,
                         top->idef.iparams,top->idef.il,
                         fr->ePBC,fr->bMolPBC,graph,cr,state->box);

        if (graph != NULL)
        {
          unshift_self(graph,state->box,state->x);
        }
        wallcycle_stop(wcycle,ewcVSITECONSTR);
      }

      /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints ############ */
      if (ir->nstlist == -1 && bFirstIterate)
      {
        gs.sig[eglsNABNSB] = nlh.nabnsb;
      }
      compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                      wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                      constr,
                      bFirstIterate ? &gs : NULL,
                      (step_rel % gs.nstms == 0) &&
                          (multisim_nsteps<0 || (step_rel<multisim_nsteps)),
                      lastbox,
                      top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                      cglo_flags
                      | (!EI_VV(ir->eI) ? CGLO_ENERGY : 0)
                      | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                      | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
                      | (bIterations && iterate.bIterate ? CGLO_ITERATE : 0)
                      | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                      | CGLO_CONSTRAINT
          );
      if (ir->nstlist == -1 && bFirstIterate)
      {
        nlh.nabnsb = gs.set[eglsNABNSB];
        gs.set[eglsNABNSB] = 0;
      }
      /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
      /* #############  END CALC EKIN AND PRESSURE ################# */

      /* Note: this is OK, but there are some numerical precision issues with using the convergence of
         the virial that should probably be addressed eventually. state->veta has better properies,
         but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
         generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

      if (bIterations &&
          done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                         trace(shake_vir),&tracevir))
      {
        break;
      }
      bFirstIterate = FALSE;
    }

    update_box(fplog,step,ir,mdatoms,state,graph,f,
               ir->nstlist==-1 ? &nlh.scale_tot : NULL,pcoupl_mu,nrnb,wcycle,upd,bInitStep,FALSE);

    /* ################# END UPDATE STEP 2 ################# */
    /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

    if (!bGStat)
    {
      /* We will not sum ekinh_old,
       * so signal that we still have to do it.
       */
      bSumEkinhOld = TRUE;
    }

    /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

    /* sum up the foreign energy and dhdl terms */
    sum_dhdl(enerd,state->lambda,ir);

    /* use the directly determined last velocity, not actually the averaged half steps */
    if (bTrotter && ir->eI==eiVV)
    {
      enerd->term[F_EKIN] = last_ekin;
    }
    enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

    if (bVV)
    {
      enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
    }
    else
    {
      enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + compute_conserved_from_auxiliary(ir,state,&MassQ);
    }
    /* #########  END PREPARING EDR OUTPUT  ###########  */

    /* update data, temperature and change at->scale */
    die_if (0 != agox_move(ago, f, enerd, step, bFirstStep, bLastStep, cr),
        "agox_move: failed " llfmt "\n", step);

    /* Time for performance */
    if (((step % stepout) == 0) || bLastStep)
    {
      runtime_upd_proc(runtime);
    }

    /* Output stuff */
    if (MASTER(cr))
    {
      gmx_bool do_dr,do_or;

      if (!(bStartingFromCpt && (EI_VV(ir->eI))))
      {
        if (bNstEner)
        {
          upd_mdebin(mdebin,bDoDHDL, TRUE,
                     t,mdatoms->tmass,enerd,state,lastbox,
                     shake_vir,force_vir,total_vir,pres,
                     ekind,mu_tot,constr);
        }
        else
        {
          upd_mdebin_step(mdebin);
        }

        do_dr  = do_per_step(step,ir->nstdisreout);
        do_or  = do_per_step(step,ir->nstorireout);

        print_ebin(outf->fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,
                   step,t,
                   eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
      }
      if (ir->ePull != epullNO)
      {
        pull_print_output(ir->pull,step,t);
      }

      if (do_per_step(step,ir->nstlog))
      {
        if(fflush(fplog) != 0)
        {
          gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
        }
      }
    }


    /* Remaining runtime */
    if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal() ))
    {
      if (shellfc)
      {
        fprintf(stderr,"\n");
      }
      print_time(stderr,runtime,step,ir,cr);
    }

    bFirstStep = FALSE;
    bInitStep = FALSE;
    bStartingFromCpt = FALSE;

    /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
    /* With all integrators, except VV, we need to retain the pressure
     * at the current step for coupling at the next step.
     */
    if ((state->flags & (1<<estPRES_PREV)) &&
        (bGStatEveryStep ||
         (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
    {
      /* Store the pressure in t_state for pressure coupling
       * at the next MD step.
       */
      copy_mat(pres,state->pres_prev);
    }

    /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

    if (bRerunMD)
    {
      if (MASTER(cr))
      {
        /* read next frame from input trajectory */
        bNotLastFrame = read_next_frame(oenv,status,&rerun_fr);
      }

      if (PAR(cr))
      {
        rerun_parallel_comm(cr,&rerun_fr,&bNotLastFrame);
      }
    }

    if (!bRerunMD || !rerun_fr.bStep)
    {
      /* increase the MD step number */
      step++;
      step_rel++;
    }

    cycles = wallcycle_stop(wcycle,ewcSTEP);
    if (DOMAINDECOMP(cr) && wcycle)
    {
      dd_cycles_add(cr->dd,(real) cycles,ddCyclStep);
    }

    if (step_rel == wcycle_get_reset_counters(wcycle) ||
        gs.set[eglsRESETCOUNTERS] != 0)
    {
      /* Reset all the counters related to performance over the run */
      reset_all_counters(fplog,cr,step,&step_rel,ir,wcycle,nrnb,runtime);
      wcycle_set_reset_counters(wcycle,-1);
      /* Correct max_hours for the elapsed time */
      max_hours -= (real) (run_time/(60.0*60.0));
      bResetCountersHalfMaxH = FALSE;
      gs.set[eglsRESETCOUNTERS] = 0;
    }

  }
  /* End of main MD loop */

  /* Stop the time */
  runtime_end(runtime);

  if (bRerunMD && MASTER(cr))
  {
    close_trj(status);
  }

  if (!(cr->duty & DUTY_PME))
  {
    /* Tell the PME only node to finish */
    gmx_pme_finish(cr);
  }

  if (MASTER(cr))
  {
    if (ir->nstcalcenergy > 0 && !bRerunMD)
    {
      print_ebin(outf->fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                 eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
    }
  }

  done_mdoutf(outf);


  if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
  {
    fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",nlh.s1/nlh.nns,sqrt(nlh.s2/nlh.nns - dblsqr(nlh.s1/nlh.nns)));
    fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",nlh.ab/nlh.nns);
  }

  if (shellfc && fplog)
  {
    fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
            (nconverged*100.0)/step_rel);
    fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
            tcount/step_rel);
  }

  runtime->nsteps_done = step_rel;

  return 0;
}

#include <signal.h>
#include <stdlib.h>


#include "tpxio.h"
#include "txtdump.h"

/* declare runner() before mdrunner_start_fn() uses it */
static
int runner(int nthreads_requested, FILE *fplog,t_commrec *cr,int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
             int repl_ex_seed, real pforce,real cpt_period,real max_hours,
             const char *deviceOptions, unsigned long Flags, int agomode);







gmx_large_int_t     deform_init_init_step_tpx;
matrix              deform_init_box_tpx;
#ifdef GMX_THREADS
tMPI_Thread_mutex_t deform_init_box_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


#ifdef GMX_THREADS
struct mdrunner_arglist
{
  FILE *fplog;
  t_commrec *cr;
  int nfile;
  const t_filenm *fnm;
  output_env_t oenv;
  gmx_bool bVerbose;
  gmx_bool bCompact;
  int nstglobalcomm;
  ivec ddxyz;
  int dd_node_order;
  real rdd;
  real rconstr;
  const char *dddlb_opt;
  real dlb_scale;
  const char *ddcsx;
  const char *ddcsy;
  const char *ddcsz;
  int nstepout;
  int resetstep;
  int nmultisim;
  int repl_ex_nst;
  int repl_ex_seed;
  real pforce;
  real cpt_period;
  real max_hours;
  const char *deviceOptions;
  unsigned long Flags;
  int agomode;
  int ret; /* return value */
};


/* The function used for spawning threads. Extracts the runner()
   arguments from its one argument and calls runner(), after making
   a commrec. */
static void mdrunner_start_fn(void *arg)
{
  struct mdrunner_arglist *mda=(struct mdrunner_arglist*)arg;
  struct mdrunner_arglist mc=*mda; /* copy the arg list to make sure
                                      that it's thread-local. This doesn't
                                      copy pointed-to items, of course,
                                      but those are all const. */
  t_commrec *cr;  /* we need a local version of this */
  FILE *fplog=NULL;
  t_filenm *fnm;

  fnm = dup_tfn(mc.nfile, mc.fnm);

  cr = init_par_threads(mc.cr);

  if (MASTER(cr))
  {
    fplog=mc.fplog;
  }

  mda->ret=runner(cr->nnodes, fplog, cr, mc.nfile, fnm, mc.oenv,
                    mc.bVerbose, mc.bCompact, mc.nstglobalcomm,
                    mc.ddxyz, mc.dd_node_order, mc.rdd,
                    mc.rconstr, mc.dddlb_opt, mc.dlb_scale,
                    mc.ddcsx, mc.ddcsy, mc.ddcsz, mc.nstepout, mc.resetstep,
                    mc.nmultisim, mc.repl_ex_nst, mc.repl_ex_seed, mc.pforce,
                    mc.cpt_period, mc.max_hours, mc.deviceOptions, mc.Flags, mc.agomode);
}

/* called by runner() to start a specific number of threads (including
   the main thread) for thread-parallel runs. This in turn calls runner()
   for each thread.
   All options besides nthreads are the same as for runner(). */
static t_commrec *mdrunner_start_threads(int nthreads,
              FILE *fplog,t_commrec *cr,int nfile,
              const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
              gmx_bool bCompact, int nstglobalcomm,
              ivec ddxyz,int dd_node_order,real rdd,real rconstr,
              const char *dddlb_opt,real dlb_scale,
              const char *ddcsx,const char *ddcsy,const char *ddcsz,
              int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
              int repl_ex_seed, real pforce,real cpt_period, real max_hours,
              const char *deviceOptions, unsigned long Flags, int agomode)
{
  int ret;
  struct mdrunner_arglist *mda;
  t_commrec *crn; /* the new commrec */
  t_filenm *fnmn;

  /* first check whether we even need to start tMPI */
  if (nthreads<2)
      return cr;

  /* a few small, one-time, almost unavoidable memory leaks: */
  snew(mda,1);
  fnmn=dup_tfn(nfile, fnm);

  /* fill the data structure to pass as void pointer to thread start fn */
  mda->fplog=fplog;
  mda->cr=cr;
  mda->nfile=nfile;
  mda->fnm=fnmn;
  mda->oenv=oenv;
  mda->bVerbose=bVerbose;
  mda->bCompact=bCompact;
  mda->nstglobalcomm=nstglobalcomm;
  mda->ddxyz[XX]=ddxyz[XX];
  mda->ddxyz[YY]=ddxyz[YY];
  mda->ddxyz[ZZ]=ddxyz[ZZ];
  mda->dd_node_order=dd_node_order;
  mda->rdd=rdd;
  mda->rconstr=rconstr;
  mda->dddlb_opt=dddlb_opt;
  mda->dlb_scale=dlb_scale;
  mda->ddcsx=ddcsx;
  mda->ddcsy=ddcsy;
  mda->ddcsz=ddcsz;
  mda->nstepout=nstepout;
  mda->resetstep=resetstep;
  mda->nmultisim=nmultisim;
  mda->repl_ex_nst=repl_ex_nst;
  mda->repl_ex_seed=repl_ex_seed;
  mda->pforce=pforce;
  mda->cpt_period=cpt_period;
  mda->max_hours=max_hours;
  mda->deviceOptions=deviceOptions;
  mda->Flags=Flags;
  mda->agomode = agomode;

  fprintf(stderr, "Starting %d threads\n",nthreads);
  fflush(stderr);
  /* now spawn new threads that start mdrunner_start_fn(), while
     the main thread returns */
  ret=tMPI_Init_fn(TRUE, nthreads, mdrunner_start_fn, (void*)(mda) );
  if (ret!=TMPI_SUCCESS)
      return NULL;

  /* make a new comm_rec to reflect the new situation */
  crn=init_par_threads(cr);
  return crn;
}


/* get the number of threads based on how many there were requested,
   which algorithms we're using, and how many particles there are. */
static int get_nthreads(int nthreads_requested, t_inputrec *inputrec,
                        gmx_mtop_t *mtop)
{
  int nthreads,nthreads_new;
  int min_atoms_per_thread;
  char *env;

  nthreads = nthreads_requested;

  /* determine # of hardware threads. */
  if (nthreads_requested < 1)
  {
    if ((env = getenv("GMX_MAX_THREADS")) != NULL)
    {
      nthreads = 0;
      sscanf(env,"%d",&nthreads);
      if (nthreads < 1)
      {
        gmx_fatal(FARGS,"GMX_MAX_THREADS (%d) should be larger than 0",
                  nthreads);
      }
    }
    else
    {
      nthreads = tMPI_Thread_get_hw_number();
    }
  }

  if (inputrec->eI == eiNM || EI_TPI(inputrec->eI))
  {
    /* Steps are divided over the nodes iso splitting the atoms */
    min_atoms_per_thread = 0;
  }
  else
  {
    min_atoms_per_thread = MIN_ATOMS_PER_THREAD;
  }

  /* Check if an algorithm does not support parallel simulation.  */
  if (nthreads != 1 &&
      ( inputrec->eI == eiLBFGS ||
        inputrec->coulombtype == eelEWALD ) )
  {
    fprintf(stderr,"\nThe integration or electrostatics algorithm doesn't support parallel runs. Not starting any threads.\n");
    nthreads = 1;
  }
  else if (nthreads_requested < 1 &&
           mtop->natoms/nthreads < min_atoms_per_thread)
  {
    /* the thread number was chosen automatically, but there are too many
       threads (too few atoms per thread) */
    nthreads_new = max(1,mtop->natoms/min_atoms_per_thread);

    if (nthreads_new > 8 || (nthreads == 8 && nthreads_new > 4))
    {
      /* Use only multiples of 4 above 8 threads
       * or with an 8-core processor
       * (to avoid 6 threads on 8 core processors with 4 real cores).
       */
      nthreads_new = (nthreads_new/4)*4;
    }
    else if (nthreads_new > 4)
    {
      /* Avoid 5 or 7 threads */
      nthreads_new = (nthreads_new/2)*2;
    }

    nthreads = nthreads_new;

    fprintf(stderr,"\n");
    fprintf(stderr,"NOTE: Parallelization is limited by the small number of atoms,\n");
    fprintf(stderr,"      only starting %d threads.\n",nthreads);
    fprintf(stderr,"      You can use the -nt option to optimize the number of threads.\n\n");
  }
  return nthreads;
}
#endif


int runner(int nthreads_requested, FILE *fplog,t_commrec *cr,int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
             int repl_ex_seed, real pforce,real cpt_period,real max_hours,
             const char *deviceOptions, unsigned long Flags, int agomode)
{
  t_inputrec *inputrec;
  t_state    *state=NULL;
  matrix     box;
  gmx_ddbox_t ddbox={0};
  int        npme_major,npme_minor;
  t_nrnb     *nrnb;
  gmx_mtop_t *mtop=NULL;
  t_mdatoms  *mdatoms=NULL;
  t_forcerec *fr=NULL;
  t_fcdata   *fcd=NULL;
  real       ewaldcoeff=0;
  gmx_pme_t  *pmedata=NULL;
  gmx_vsite_t *vsite=NULL;
  gmx_constr_t constr;
  int        nChargePerturbed=-1,status;
  gmx_wallcycle_t wcycle;
  gmx_bool       bReadRNG,bReadEkin;
  gmx_runtime_t runtime;
  int        rc;
  gmx_large_int_t reset_counters;
  gmx_edsam_t ed=NULL;
  t_commrec   *cr_old=cr;
  int         nthreads=1;
  ago_t *ago;

  /* CAUTION: threads may be started later on in this function, so
     cr doesn't reflect the final parallel state right now */
  snew(inputrec,1);
  snew(mtop,1);

  if (Flags & MD_APPENDFILES)
  {
    fplog = NULL;
  }

  snew(state,1);
  if (MASTER(cr))
  {
    /* Read (nearly) all data required for the simulation */
    read_tpx_state(ftp2fn(efTPX,nfile,fnm),inputrec,state,NULL,mtop);

    /* NOW the threads will be started: */
#ifdef GMX_THREADS
    nthreads = get_nthreads(nthreads_requested, inputrec, mtop);

    if (nthreads > 1)
    {
      /* now start the threads. */
      cr=mdrunner_start_threads(nthreads, fplog, cr_old, nfile, fnm,
                                oenv, bVerbose, bCompact, nstglobalcomm,
                                ddxyz, dd_node_order, rdd, rconstr,
                                dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                nstepout, resetstep, nmultisim,
                                repl_ex_nst, repl_ex_seed, pforce,
                                cpt_period, max_hours, deviceOptions,
                                Flags, agomode);
      /* the main thread continues here with a new cr. We don't deallocate
         the old cr because other threads may still be reading it. */
      if (cr == NULL)
      {
        gmx_comm("Failed to spawn threads");
      }
    }
#endif
  }
  /* END OF CAUTION: cr is now reliable */

  if (PAR(cr))
  {
    /* now broadcast everything to the non-master nodes/threads: */
    init_parallel(fplog, cr, inputrec, mtop);
  }
  if (fplog != NULL)
  {
    pr_inputrec(fplog,0,"Input Parameters",inputrec,FALSE);
  }

  /* now make sure the state is initialized and propagated */
  set_state_entries(state,inputrec,cr->nnodes);

  /* A parallel command line option consistency check that we can
     only do after any threads have started. */
  if (!PAR(cr) &&
      (ddxyz[XX] > 1 || ddxyz[YY] > 1 || ddxyz[ZZ] > 1 || cr->npmenodes > 0))
  {
    gmx_fatal(FARGS,
              "The -dd or -npme option request a parallel simulation, "
#ifndef GMX_MPI
              "but mdrun was compiled without threads or MPI enabled"
#else
#ifdef GMX_THREADS
              "but the number of threads (option -nt) is 1"
#else
              "but mdrun was not started through mpirun/mpiexec or only one process was requested through mpirun/mpiexec"
#endif
#endif
        );
  }

  if ((Flags & MD_RERUN) &&
      (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
  {
    gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
  }

  if (can_use_allvsall(inputrec,mtop,TRUE,cr,fplog))
  {
    /* All-vs-all loops do not work with domain decomposition */
    Flags |= MD_PARTDEC;
  }

  if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC))
  {
    cr->npmenodes = 0;
  }

#ifdef GMX_FAHCORE
  fcRegisterSteps(inputrec->nsteps,inputrec->init_step);
#endif

  /* NMR restraints must be initialized before load_checkpoint,
   * since with time averaging the history is added to t_state.
   * For proper consistency check we therefore need to extend
   * t_state here.
   * So the PME-only nodes (if present) will also initialize
   * the distance restraints.
   */
  snew(fcd,1);

  /* This needs to be called before read_checkpoint to extend the state */
  init_disres(fplog,mtop,inputrec,cr,Flags & MD_PARTDEC,fcd,state);

  if (gmx_mtop_ftype_count(mtop,F_ORIRES) > 0)
  {
    if (PAR(cr) && !(Flags & MD_PARTDEC))
    {
      gmx_fatal(FARGS,"Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
    }
    /* Orientation restraints */
    if (MASTER(cr))
    {
      init_orires(fplog,mtop,state->x,inputrec,cr->ms,&(fcd->orires),
                  state);
    }
  }

  if (DEFORM(*inputrec))
  {
    /* Store the deform reference box before reading the checkpoint */
    if (SIMMASTER(cr))
    {
      copy_mat(state->box,box);
    }
    if (PAR(cr))
    {
      gmx_bcast(sizeof(box),box,cr);
    }
    /* Because we do not have the update struct available yet
     * in which the reference values should be stored,
     * we store them temporarily in static variables.
     * This should be thread safe, since they are only written once
     * and with identical values.
     */
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
    deform_init_init_step_tpx = inputrec->init_step;
    copy_mat(box,deform_init_box_tpx);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
  }

  if (opt2bSet("-cpi",nfile,fnm))
  {
    /* Check if checkpoint file exists before doing continuation.
     * This way we can use identical input options for the first and subsequent runs...
     */
    if( gmx_fexist_master(opt2fn_master("-cpi",nfile,fnm,cr),cr) )
    {
      load_checkpoint(opt2fn_master("-cpi",nfile,fnm,cr),&fplog,
                      cr,Flags & MD_PARTDEC,ddxyz,
                      inputrec,state,&bReadRNG,&bReadEkin,
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

  if (((MASTER(cr) || (Flags & MD_SEPPOT)) && (Flags & MD_APPENDFILES))
#ifdef GMX_THREADS
      /* With thread MPI only the master node/thread exists in mdrun.c,
       * therefore non-master nodes need to open the "seppot" log file here.
       */
      || (!MASTER(cr) && (Flags & MD_SEPPOT))
#endif
      )
  {
    gmx_log_open(ftp2fn(efLOG,nfile,fnm),cr,!(Flags & MD_SEPPOT),
                         Flags,&fplog);
  }

  if (SIMMASTER(cr))
  {
    copy_mat(state->box,box);
  }

  if (PAR(cr))
  {
    gmx_bcast(sizeof(box),box,cr);
  }

  /* Essential dynamics */
  if (opt2bSet("-ei",nfile,fnm))
  {
    /* Open input and output files, allocate space for ED data structure */
    ed = ed_open(nfile,fnm,Flags,cr);
  }

  if (PAR(cr) && !((Flags & MD_PARTDEC) ||
                   EI_TPI(inputrec->eI) ||
                   inputrec->eI == eiNM))
  {
    cr->dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
                                       dddlb_opt,dlb_scale,
                                       ddcsx,ddcsy,ddcsz,
                                       mtop,inputrec,
                                       box,state->x,
                                       &ddbox,&npme_major,&npme_minor);

    make_dd_communicators(fplog,cr,dd_node_order);

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);
  }
  else
  {
    /* PME, if used, is done on all nodes with 1D decomposition */
    cr->npmenodes = 0;
    cr->duty = (DUTY_PP | DUTY_PME);
    npme_major = 1;
    npme_minor = 1;
    if (!EI_TPI(inputrec->eI))
    {
      npme_major = cr->nnodes;
    }

    if (inputrec->ePBC == epbcSCREW)
    {
      gmx_fatal(FARGS,
                "pbc=%s is only implemented with domain decomposition",
                epbc_names[inputrec->ePBC]);
    }
  }

  if (PAR(cr))
  {
    /* After possible communicator splitting in make_dd_communicators.
     * we can set up the intra/inter node communication.
     */
    gmx_setup_nodecomm(fplog,cr);
  }

  wcycle = wallcycle_init(fplog,resetstep,cr);
  if (PAR(cr))
  {
    /* Master synchronizes its value of reset_counters with all nodes
     * including PME only nodes */
    reset_counters = wcycle_get_reset_counters(wcycle);
    gmx_bcast_sim(sizeof(reset_counters),&reset_counters,cr);
    wcycle_set_reset_counters(wcycle, reset_counters);
  }


  snew(nrnb,1);
  if (cr->duty & DUTY_PP)
  {
    /* For domain decomposition we allocate dynamically
     * in dd_partition_system.
     */
    if (DOMAINDECOMP(cr))
    {
      bcast_state_setup(cr,state);
    }
    else
    {
      if (PAR(cr))
      {
        bcast_state(cr,state,TRUE);
      }
    }

    /* Dihedral Restraints */
    if (gmx_mtop_ftype_count(mtop,F_DIHRES) > 0)
    {
      init_dihres(fplog,mtop,inputrec,fcd);
    }

    /* Initiate forcerecord */
    fr = mk_forcerec();
    init_forcerec(fplog,oenv,fr,fcd,inputrec,mtop,cr,box,FALSE,
                  opt2fn("-table",nfile,fnm),
                  opt2fn("-tablep",nfile,fnm),
                  opt2fn("-tableb",nfile,fnm),FALSE,pforce);

    /* version for PCA_NOT_READ_NODE (see md.c) */
    /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
      "nofile","nofile","nofile",FALSE,pforce);
      */
    fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);

    /* Initialize QM-MM */
    if(fr->bQMMM)
    {
      init_QMMMrec(cr,box,mtop,inputrec,fr);
    }

    /* Initialize the mdatoms structure.
     * mdatoms is not filled with atom data,
     * as this can not be done now with domain decomposition.
     */
    mdatoms = init_mdatoms(fplog,mtop,inputrec->efep!=efepNO);

    /* Initialize the virtual site communication */
    vsite = init_vsite(mtop,cr);

    calc_shifts(box,fr->shift_vec);

    /* With periodic molecules the charge groups should be whole at start up
     * and the virtual sites should not be far from their proper positions.
     */
    if (!inputrec->bContinuation && MASTER(cr) &&
        !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
    {
      /* Make molecules whole at start of run */
      if (fr->ePBC != epbcNONE)
      {
        do_pbc_first_mtop(fplog,inputrec->ePBC,box,mtop,state->x);
      }
      if (vsite)
      {
        /* Correct initial vsite positions are required
         * for the initial distribution in the domain decomposition
         * and for the initial shell prediction.
         */
        construct_vsites_mtop(fplog,vsite,mtop,state->x);
      }
    }

    /* Initiate PPPM if necessary */
    if (fr->eeltype == eelPPPM)
    {
      if (mdatoms->nChargePerturbed)
      {
        gmx_fatal(FARGS,"Free energy with %s is not implemented",
                  eel_names[fr->eeltype]);
      }
      status = gmx_pppm_init(fplog,cr,oenv,FALSE,TRUE,box,
                             getenv("GMXGHAT"),inputrec, (Flags & MD_REPRODUCIBLE));
      if (status != 0)
      {
        gmx_fatal(FARGS,"Error %d initializing PPPM",status);
      }
    }

    if (EEL_PME(fr->eeltype))
    {
      ewaldcoeff = fr->ewaldcoeff;
      pmedata = &fr->pmedata;
    }
    else
    {
      pmedata = NULL;
    }
  }
  else
  {
    /* This is a PME only node */

    /* We don't need the state */
    done_state(state);

    ewaldcoeff = calc_ewaldcoeff(inputrec->rcoulomb, inputrec->ewald_rtol);
    snew(pmedata,1);
  }

  /* Initiate PME if necessary,
   * either on all nodes or on dedicated PME nodes only. */
  if (EEL_PME(inputrec->coulombtype))
  {
    if (mdatoms)
    {
      nChargePerturbed = mdatoms->nChargePerturbed;
    }
    if (cr->npmenodes > 0)
    {
      /* The PME only nodes need to know nChargePerturbed */
      gmx_bcast_sim(sizeof(nChargePerturbed),&nChargePerturbed,cr);
    }
    if (cr->duty & DUTY_PME)
    {
      status = gmx_pme_init(pmedata,cr,npme_major,npme_minor,inputrec,
                            mtop ? mtop->natoms : 0,nChargePerturbed,
                            (Flags & MD_REPRODUCIBLE));
      if (status != 0)
      {
        gmx_fatal(FARGS,"Error %d initializing PME",status);
      }
    }
  }


  {
    /* Turn on signal handling on all nodes */
    /*
     * (A user signal from the PME nodes (if any)
     * is communicated to the PP nodes.
     */
    signal_handler_install();
  }

  /* initialize ago, cr->duty PP/PME has been assigned */
  ago = agox_init(agox_opt2fn("-ago", nfile, fnm),
      Flags & MD_STARTFROMCPT, state,
      mtop, inputrec, cr, agomode);
  if ((cr->duty & DUTY_PP) && ago == NULL) {
    fprintf(stderr, "Failed to initialize aT\n");
    exit(1);
  }

  if (cr->duty & DUTY_PP)
  {
    if (inputrec->ePull != epullNO)
    {
      /* Initialize pull code */
      init_pull(fplog,inputrec,nfile,fnm,mtop,cr,oenv,
                EI_DYNAMICS(inputrec->eI) && MASTER(cr),Flags);
    }

    constr = init_constraints(fplog,mtop,inputrec,ed,state,cr);

    if (DOMAINDECOMP(cr))
    {
      dd_init_bondeds(fplog,cr->dd,mtop,vsite,constr,inputrec,
                      Flags & MD_DDBONDCHECK,fr->cginfo_mb);

      set_dd_parameters(fplog,cr->dd,dlb_scale,inputrec,fr,&ddbox);

      setup_dd_grid(fplog,cr->dd);
    }

    /* Now do whatever the user wants us to do (how flexible...) */
    md(fplog,cr,nfile,fnm,
       oenv,bVerbose,bCompact,
       nstglobalcomm,
       vsite,constr,
       nstepout,inputrec,mtop,
       fcd,state,
       mdatoms,nrnb,wcycle,ed,fr,
       repl_ex_nst,repl_ex_seed,
       cpt_period,max_hours,
       deviceOptions,
       Flags,
       &runtime, ago);

    if (inputrec->ePull != epullNO)
    {
      finish_pull(fplog,inputrec->pull);
    }
  }
  else
  {
    /* do PME only */
    gmx_pmeonly(*pmedata,cr,nrnb,wcycle,ewaldcoeff,FALSE,inputrec);
  }

  if (EI_DYNAMICS(inputrec->eI) || EI_TPI(inputrec->eI))
  {
    /* Some timing stats */
    if (SIMMASTER(cr))
    {
      if (runtime.proc == 0)
      {
        runtime.proc = (clock_t) runtime.real;
      }
    }
    else
    {
      runtime.real = 0;
    }
  }

  wallcycle_stop(wcycle,ewcRUN);

  /* Finish up, write some stuff
   * if rerunMD, don't write last frame again
   */
  finish_run(fplog,cr,ftp2fn(efSTO,nfile,fnm),
             inputrec,nrnb,wcycle,&runtime,
             EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

  /* Does what it says */
  print_date_and_time(fplog,cr->nodeid,"Finished mdrun",&runtime);

  /* Close logfile already here if we were appending to it */
  if (MASTER(cr) && (Flags & MD_APPENDFILES))
  {
    gmx_log_close(fplog);
  }

  rc=(int)gmx_get_stop_condition();

#ifdef GMX_THREADS
  /* we need to join all threads. The sub-threads join when they
     exit this function, but the master thread needs to be told to
     wait for that. */
  if (PAR(cr) && MASTER(cr))
  {
    tMPI_Finalize();
  }
#endif

  if(ago != NULL)
    agox_close(ago); /* release memory */
  return rc;
}

#include "macros.h"
#include "main.h"
#include "futil.h"
#include "edsam.h"
#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif


int main(int argc,char *argv[])
{
  const char *desc[] = {"self-contained GROMACS"};
  t_commrec    *cr;
  t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efXTC, "-x",      NULL,       ffOPTWR },
    { efCPT, "-cpi",    NULL,       ffOPTRD },
    { efCPT, "-cpo",    NULL,       ffOPTWR },
    { efSTO, "-c",      "confout",  ffWRITE },
    { efEDR, "-e",      "ener",     ffWRITE },
    { efLOG, "-g",      "md",       ffWRITE },
    { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
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
    { efRND, "-multidir",NULL,      ffOPTRDMULT},
    { efMDP, "-ago",     NULL,      ffOPTRD } /* ago.cfg */
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  gmx_bool bPartDec     = FALSE;
  gmx_bool bDDBondCheck = TRUE;
  gmx_bool bDDBondComm  = TRUE;
  gmx_bool bVerbose     = FALSE;
  gmx_bool bCompact     = TRUE;
  gmx_bool bSepPot      = FALSE;
  gmx_bool bRerunVSite  = FALSE;
  gmx_bool bConfout     = TRUE;
  gmx_bool bReproducible = FALSE;

  int  npme=-1;
  int  nmultisim=0;
  int  nstglobalcomm=-1;
  int  repl_ex_nst=0;
  int  repl_ex_seed=-1;
  int  nstepout=100;
  int  nthreads=0; /* set to determine # of threads automatically */
  int  resetstep=-1;

  rvec realddxyz={0,0,0};
  const char *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
    const char *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
  real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
  char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
  real cpt_period=15.0,max_hours=-1;
  gmx_bool bAppendFiles=TRUE;
  gmx_bool bKeepAndNumCPT=FALSE;
  gmx_bool bResetCountersHalfWay=FALSE;
  output_env_t oenv=NULL;
  const char *deviceOptions = "";

  int agomode = 0;
  t_pargs pa[] = {

    { "-pd",      FALSE, etBOOL,{&bPartDec},
      "Use particle decompostion" },
    { "-dd",      FALSE, etRVEC,{&realddxyz},
      "Domain decomposition grid, 0 is optimize" },
#ifdef GMX_THREADS
    { "-nt",      FALSE, etINT, {&nthreads},
      "Number of threads to start (0 is guess)" },
#endif
    { "-npme",    FALSE, etINT, {&npme},
      "Number of separate nodes to be used for PME, -1 is guess" },
    { "-ddorder", FALSE, etENUM, {ddno_opt},
      "DD node order" },
    { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
      "Check for all bonded interactions with DD" },
    { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
      "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
    { "-rdd",     FALSE, etREAL, {&rdd},
      "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
    { "-rcon",    FALSE, etREAL, {&rconstr},
      "Maximum distance for P-LINCS (nm), 0 is estimate" },
    { "-dlb",     FALSE, etENUM, {dddlb_opt},
      "Dynamic load balancing (with DD)" },
    { "-dds",     FALSE, etREAL, {&dlb_scale},
      "Minimum allowed dlb scaling of the DD cell size" },
    { "-ddcsx",   FALSE, etSTR, {&ddcsx},
      "HIDDENThe DD cell sizes in x" },
    { "-ddcsy",   FALSE, etSTR, {&ddcsy},
      "HIDDENThe DD cell sizes in y" },
    { "-ddcsz",   FALSE, etSTR, {&ddcsz},
      "HIDDENThe DD cell sizes in z" },
    { "-gcom",    FALSE, etINT,{&nstglobalcomm},
      "Global communication frequency" },
    { "-v",       FALSE, etBOOL,{&bVerbose},
      "Be loud and noisy" },
    { "-compact", FALSE, etBOOL,{&bCompact},
      "Write a compact log file" },
    { "-seppot",  FALSE, etBOOL, {&bSepPot},
      "Write separate V and dVdl terms for each interaction type and node to the log file(s)" },
    { "-pforce",  FALSE, etREAL, {&pforce},
      "Print all forces larger than this (kJ/mol nm)" },
    { "-reprod",  FALSE, etBOOL,{&bReproducible},
      "Try to avoid optimizations that affect binary reproducibility" },
    { "-cpt",     FALSE, etREAL, {&cpt_period},
      "Checkpoint interval (minutes)" },
    { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
      "Keep and number checkpoint files" },
    { "-append",  FALSE, etBOOL, {&bAppendFiles},
      "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
    { "-maxh",   FALSE, etREAL, {&max_hours},
      "Terminate after 0.99 times this time (hours)" },
    { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
      "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
    { "-confout", FALSE, etBOOL, {&bConfout},
      "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
    { "-stepout", FALSE, etINT, {&nstepout},
      "HIDDENFrequency of writing the remaining runtime" },
    { "-resetstep", FALSE, etINT, {&resetstep},
      "HIDDENReset cycle counters after these many time steps" },
    { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
      "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" },
    { "-mode", FALSE, etINT, {&agomode},
      "a preparation run" }
  };
  unsigned long Flags, PCA_Flags;
  ivec     ddxyz;
  int      dd_node_order;
  gmx_bool     bAddPart;
  FILE     *fplog;
  int      sim_part,sim_part_fn;
  const char *part_suffix=".part";
  char     suffix[STRLEN];
  int      rc;
  char **multidir=NULL;


  cr = init_par(&argc,&argv);

  PCA_Flags = (PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_CAN_SET_DEFFNM
               | (MASTER(cr) ? 0 : PCA_QUIET));


  /* Comment this in to do fexist calls only on master
   * works not with rerun or tables at the moment
   * also comment out the version of init_forcerec in md.c
   * with NULL instead of opt2fn
   */
  /*
     if (!MASTER(cr))
     {
     PCA_Flags |= PCA_NOT_READ_NODE;
     }
     */

  parse_common_args(&argc,argv,PCA_Flags, NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL, &oenv);



  /* we set these early because they might be used in init_multisystem()
     Note that there is the potential for npme>nnodes until the number of
     threads is set later on, if there's thread parallelization. That shouldn't
     lead to problems. */
  dd_node_order = nenum(ddno_opt);
  cr->npmenodes = npme;

#ifndef GMX_THREADS
  nthreads=1;
#endif

  /* now check the -multi and -multidir option */
  if (opt2bSet("-multidir", NFILE, fnm))
  {
    if (nmultisim > 0)
    {
      gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually exclusive.");
    }
    nmultisim = opt2fns(&multidir, "-multidir", NFILE, fnm);
  }


  bAddPart = !bAppendFiles;

  /* Check if there is ANY checkpoint file available */
  sim_part    = 1;
  sim_part_fn = sim_part;
  if (opt2bSet("-cpi",NFILE,fnm))
  {
    if (bSepPot && bAppendFiles)
    {
      gmx_fatal(FARGS,"Output file appending is not supported with -seppot");
    }

    bAppendFiles =
              read_checkpoint_simulation_part(opt2fn_master("-cpi", NFILE,
                                                            fnm,cr),
                                              &sim_part_fn,NULL,cr,
                                              bAppendFiles,NFILE,fnm,
                                              part_suffix,&bAddPart);
    if (sim_part_fn==0 && MASTER(cr))
    {
      fprintf(stdout,"No previous checkpoint file present, assuming this is a new run.\n");
    }
    else
    {
      sim_part = sim_part_fn + 1;
    }

    if (MULTISIM(cr) && MASTER(cr))
    {
      check_multi_int(stdout,cr->ms,sim_part,"simulation part");
    }
  }
  else
  {
    bAppendFiles = FALSE;
  }

  if (!bAppendFiles)
  {
    sim_part_fn = sim_part;
  }

  if (bAddPart)
  {
    /* Rename all output files (except checkpoint files) */
    /* create new part name first (zero-filled) */
    sprintf(suffix,"%s%04d",part_suffix,sim_part_fn);

    add_suffix_to_output_names(fnm,NFILE,suffix);
    if (MASTER(cr))
    {
      fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed '%s'.\n",sim_part-1,suffix);
    }
  }

  Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
  Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
  Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
  Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
  Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
  Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
  Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
  Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
  Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0);
  Flags = Flags | (bKeepAndNumCPT ? MD_KEEPANDNUMCPT : 0);
  Flags = Flags | (sim_part>1    ? MD_STARTFROMCPT : 0);
  Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);


  /* We postpone opening the log file if we are appending, so we can
     first truncate the old log file and append to the correct position
     there instead.  */
  if ((MASTER(cr) || bSepPot) && !bAppendFiles)
  {
    gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags,&fplog);
  }
  else if (!MASTER(cr) && bSepPot)
  {
    gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags,&fplog);
  }
  else
  {
    fplog = NULL;
  }

  ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
  ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
  ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);

  rc = runner(nthreads, fplog,cr,NFILE,fnm,oenv,bVerbose,bCompact,
              nstglobalcomm, ddxyz,dd_node_order,rdd,rconstr,
              dddlb_opt[0],dlb_scale,ddcsx,ddcsy,ddcsz,
              nstepout,resetstep,nmultisim,repl_ex_nst,repl_ex_seed,
              pforce, cpt_period,max_hours,deviceOptions,Flags, agomode);

  if (gmx_parallel_env_initialized())
    gmx_finalize();

  /* Log file has to be closed in runner if we are appending to it
     (fplog not set here) */
  if (MASTER(cr) && !bAppendFiles)
  {
    gmx_log_close(fplog);
  }

  return rc;
}

/* =========================================================================================

         GROMACS integration

 * ========================================================================================= */

#include "mdgocore.h"

INLINE const char * yesno(int x) { return (x) ? "yes" : "no"; }

INLINE int every(int on, llong_t step, int freq)
  { return on && (freq > 0) && (step % freq == 0); }

INLINE void agox_close(ago_t *ago)
{
  if (ago->mpi_rank == 0) {
    if (ago->log) { log_close(ago->log); ago->log = NULL; }
    if (ago->tmh) { tmh_close(ago->tmh); ago->tmh = NULL; }
    mtsave(ago->fnrng);
  }
  //if (ago->atnb) free(ago->atnb);
  ago_close(ago);
}

/* check reference structure against GROMACS mtop
 * assume a single topology in mtop */
static void agox_checkxreftop(ago_t *ago, const gmx_mtop_t *mtop)
{
  t_atoms *atoms = &mtop->moltype[0].atoms;
  pdbmodel_t *pm = ago->pm;
  const char *resnm, *atnm;
  int i;

  die_if (atoms->nr != ago->nat,
      "# of atoms mismatch, top: %d, ref: %d\n", atoms->nr, ago->nat);

  die_if (atoms->nres != ago->nres,
      "# of residues mismatch, top: %d, ref: %d\n", atoms->nres, ago->nres);

  for (i = 0; i < atoms->nr; i++) {
    atnm = atoms->atomname[i][0];
    die_if (strcmp(atnm, pm->atm[i].atnm) != 0,
      "atom %d name mismatch, top: %s, ref: %s\n",
          i, atnm, pm->atm[i].atnm);

    resnm = atoms->resinfo[ atoms->atom[i].resind ].name[0];
    die_if (strcmp(resnm, pm->atm[i].resnm) != 0,
        "atom %d residue mismatch, top %s, ref: %s\n",
        i, resnm, pm->atm[i].atnm);
  }
}

/* make a neighbor list of contacting residues
 * n0: mdatoms->start
 * n1: mdatoms->start + mdatoms->homenr
 * use after atoms2md() which assigns mdatoms->start and mdatoms->homenr */
static void agox_assign_mkls(ago_t *ago, int n0, int n1, t_commrec *cr)
{
  int i, j, nbtot = 0, nat = ago->nat;
  atnbls_t *nb;

  ago->atcnt = 0;
  i = n1 - n0; if (i < 1) i = 1;
  xnew(ago->atnb, i);
  for (i = n0; i < n1; i++) {
    if (ago->wt[i] <= 0.0) continue;

    /* generate neighborlist for particle i */
    nb = ago->atnb + ago->atcnt;
    nb->id = i;
    nb->nbcnt = 0;
    for (j = 0; j < nat; j++) {
      if (!ago->isct[i*nat + j] || ago->wt[j] <= 0.0) continue;
      /* if j is an home atom, require j > i */
      if (j >= n0 && j < n1) {
        if (j < i) continue;
      } else { /* external atom */
        /* if j > i, include it only if j + i is even */
        if (j > i && (j + i) % 2 == 1) continue;
        /* if j < i, include it only if j + i is odd */
        if (j < i && (j + i) % 2 == 0) continue;
      }

      nb->nb[ nb->nbcnt++ ] = j;

      die_if (nb->nbcnt > ATNBMAX, "%d: too many neighbors for atom %d\n", cr->nodeid, i);
    }
    ago->atcnt++;
    nbtot += nb->nbcnt;
  }
  printf("node %d: has %d go atoms, %d contacts\n", cr->nodeid, ago->atcnt, nbtot);
}

/* initialize aT, every node calls this */
static ago_t *agox_init(const char *fncfg, unsigned fromcpt,
    const t_state *state,
    const gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr, int mode)
{
  ago_t *ago;

  if (!(cr->duty & DUTY_PP)) return NULL; /* for a PME-only node */

  die_if (DOMAINDECOMP(cr),
      "only support particle decomposition\n");

  if (SIMMASTER(cr)) {
    /* init. the master, passing temperature, time step, and mode */
    ago = ago_open(fncfg, fromcpt, ir->opts.ref_t[0], ir->delta_t, mode);
    agox_checkxreftop(ago, mtop);
    ago->rmsd = ago_carmsd(ago, state->x);
    ago->goep = ago_goepot(ago, state->x);
    ago_dump(ago, state->x, "mdgo.manifest", 3);
  } else { /* for non-master, simply allocate space */
    xnew(ago, 1);
  }

#ifdef GMX_MPI
  /* tell everyone settings on the master, valid only for PP only node
   * we pass MPI_COMM_NULL to avoid the case of one-node-mpi */
  ago_initmpi(ago, PAR(cr) ? cr->mpi_comm_mygroup : MPI_COMM_NULL);
#endif

  if (PAR(cr)) gmx_barrier(cr);
  return ago;
}

/* compute the force from the reference structure
 * fscal is scaling factor for the force */
double agox_goforce(ago_t *ago, rvec x[], rvec f[], real fscal,
    t_commrec *cr, llong_t step)
{
  int ii, jj, i, j, prid, nat = ago->nat;
  int check = every(1, step, ago->nstcheckf);
  real dx[3], dr2, nbe = ago->nbe;
  real ene = 0, amp;
  rv3_t *fgo = ago->fgo;

  if (every(MASTER(cr), step, ago->nstgoadj))
    ago->rmsd = ago_carmsd(ago, x);

  if (check) { /* clear go force */
    for (i = 0; i < nat; i++) rv3_zero(fgo[i]);
  }

  for (ii = 0; ii < ago->atcnt; ii++) { /* loop over home atoms */
    atnbls_t *nb = ago->atnb + ii;
    i = nb->id;
    for (jj = 0; jj < nb->nbcnt; jj++) {
      j = nb->nb[jj];
      dr2 = rv3_sqr( rv3_diff(dx, x[i], x[j]) );
      prid = i * nat + j;
      ene += nbe * ago_pairene(dr2, ago->drref[prid], &amp);
      amp *= nbe * fscal;
      rv3_sinc(f[i], dx, amp);
      rv3_sinc(f[j], dx, -amp);
      if (check) {
        rv3_sinc(fgo[i], dx, amp);
        rv3_sinc(fgo[j], dx, -amp);
      }
    }
  }
  ago->goepvdw = ene;

  /* CA dihedral energy */
  ene = 0;
  if (ago->edih > 0) {
    int nn = ago->mpi_size, ni = ago->mpi_rank, id, id0, id1, idsz, dihcnt = ago->dihcnt, i0, i1, i2, i3;
    real dih, edih = ago->edih, ampdih;
    dihcalc_t dc[1];

    idsz = (dihcnt + nn - 1) / nn;
    id0 = ni * idsz;
    id1 = (ni + 1) * idsz;
    if (id1 > dihcnt) id1 = dihcnt;
    memset(dc, '\0', sizeof(*dc));
    dc->szreal = sizeof(real);
    for (id = id0; id < id1; id++) {
      dih = rv3_calcdihv(dc, x, ago->dihid + 4 * id, DIH_GRAD|DIH_FOUR|DIH_POLYMER);
      i0 = ago->dihid[4*id];
      i1 = ago->dihid[4*id + 1];
      i2 = ago->dihid[4*id + 2];
      i3 = ago->dihid[4*id + 3];
      dr2 = rv3_sqr( rv3_diff(dx, x[i0], x[i3]) );
      ene += edih * ago_dihene(dih, ago->dihref[id], dr2, ago->dihdis[id], &ampdih, &amp);
      ampdih *= edih * fscal;
      amp *= edih * fscal;
      rv3_sinc(f[i0], dc->g[0], ampdih);
      rv3_sinc(f[i1], dc->g[1], ampdih);
      rv3_sinc(f[i2], dc->g[2], ampdih);
      rv3_sinc(f[i3], dc->g[3], ampdih);
      rv3_sinc(f[i0], dx, amp);
      rv3_sinc(f[i3], dx, -amp);
      if (check) {
        rv3_sinc(fgo[i0], dc->g[0], ampdih);
        rv3_sinc(fgo[i1], dc->g[1], ampdih);
        rv3_sinc(fgo[i2], dc->g[2], ampdih);
        rv3_sinc(fgo[i3], dc->g[3], ampdih);
        rv3_sinc(fgo[i0], dx, amp);
        rv3_sinc(fgo[i3], dx, -amp);
      }
    }
  }
  ago->goepdih = ene;

  /* dihedral bias potential (towards a single sign)
   * the energy is always fully applied no matter what fscal is */
  ene = 0;
  if (ago->dihbias > 0) {
    int nn = ago->mpi_size, ni = ago->mpi_rank, id, id0, id1, idsz, dihcnt = ago->dihcnt, i0, i1, i2, i3;
    real dih, dihbias = ago->dihbias, ampdih;
    dihcalc_t dc[1];

    idsz = (dihcnt + nn - 1) / nn;
    id0 = ni * idsz;
    id1 = (ni + 1) * idsz;
    if (id1 > dihcnt) id1 = dihcnt;
    memset(dc, '\0', sizeof(*dc));
    dc->szreal = sizeof(real);
    for (id = id0; id < id1; id++) {
      dih = rv3_calcdihv(dc, x, ago->dihid + 4 * id, DIH_GRAD|DIH_FOUR|DIH_POLYMER);
      i0 = ago->dihid[4*id];
      i1 = ago->dihid[4*id + 1];
      i2 = ago->dihid[4*id + 2];
      i3 = ago->dihid[4*id + 3];
      dr2 = rv3_sqr( rv3_diff(dx, x[i0], x[i3]) );
      ene += dihbias * ago_dihene(dih, ago->dihref[id], dr2, ago->dihdis[id], &ampdih, &amp);
      ampdih *= dihbias;
      amp *= dihbias;
      rv3_sinc(f[i0], dc->g[0], ampdih);
      rv3_sinc(f[i1], dc->g[1], ampdih);
      rv3_sinc(f[i2], dc->g[2], ampdih);
      rv3_sinc(f[i3], dc->g[3], ampdih);
      rv3_sinc(f[i0], dx, amp);
      rv3_sinc(f[i3], dx, -amp);
    }
  }
  ago->goepdihbias = ene;

  ago->goep = ago->goepvdw + ago->goepdih;
  //printf("node %d: step "llfmt": energy check = %f vs. %g\n", cr->nodeid, step, ago_goepot(ago, x), ago->goep); getchar();
  return ago->goep;
}

/* check if the force is correct */
static void agox_forcecheck(ago_t *ago, rvec x[], t_commrec *cr, llong_t step)
{
  int i, nat = ago->nat;
  double ene1, ene2;
  real f2, invf2, del = 1.0f;
  rv3_t *f;

  if (PAR(cr)) gmx_barrier(cr);
  xnew(f, nat);
  for (i = 0; i < nat; i++) rv3_zero(f[i]);
  ene1 = agox_goforce(ago, x, f, 1.0f, cr, step);
  for (f2 = 0, i = 0; i < nat; i++) f2 += rv3_sqr(f[i]);
  invf2 = del/f2;
  if (PAR(cr)) gmx_barrier(cr);
  /* move along the force */
  for (i = 0; i < nat; i++) rv3_sinc(x[i], f[i], invf2);
  ene2 = agox_goforce(ago, x, f, 1.0f, cr, step);
  for (i = 0; i < nat; i++) rv3_sinc(x[i], f[i], -invf2);
  if (PAR(cr)) gmx_barrier(cr);
  printf("%d: force check: ene1 %g ene2 %g, dene %g vs. %g\n", cr->nodeid, ene1, ene2, ene1 - ene2, del);
  free(f);
}

/* low level force checking */
static void agox_checkf0(ago_t *ago, rv3_t *f,
    real big, const char *tag)
{
  int i, imax = 0, nat = ago->nat;
  real f2, fmax = 0, big2 = big * big, fc[3];

  rv3_zero(fc);
  for (i = 0; i < nat; i++) {
    rv3_inc(fc, f[i]);
    f2 = rv3_sqr(f[i]);
    if (f2 > fmax) { fmax = f2; imax = i; }
    if (f2 > big2)
      fprintf(stderr, "large force on atom %d: %g, %g %g\n",
          i + 1, f[i][0], f[i][1], f[i][2]);
  }
  rv3_smul(fc, 1.0f/nat);
  fmax = (real) sqrt(fmax);
  if (tag)
    printf("%4s: max. %4d: %12.4e, fc %+8.3f, %+8.3f, %+8.3f\n",
      tag, imax + 1, fmax, fc[0], fc[1], fc[2]);
}

/* check force */
static void agox_checkf(ago_t *ago, rv3_t *f)
{
  agox_checkf0(ago, f, 2e4f, "f");
  agox_checkf0(ago, ago->fgo, 2e3f, "fgo");
}

/* only a PP-processor calls this function
 * also assume global_stat() has been called
 * SIMMASTER(cr) should be equivalent to MASTER(cr) */
int agox_move(ago_t *ago, rv3_t *f, gmx_enerdata_t *enerd,
    llong_t step, int bFirstStep, int bLastStep, t_commrec *cr)
{
  int ismaster = MASTER(cr);

  die_if(!(cr->duty & DUTY_PP), "pme node %d cannot call agox_sum\n", cr->nodeid);

  if (!every(1, step, ago->nstgoadj)) return 0;

  ago->H0 = enerd->term[F_EPOT];
  if (PAR(cr)) {
    gmx_sumd(1, &ago->goepdih, cr);
    gmx_sumd(1, &ago->goepdihbias, cr);
    gmx_sumd(1, &ago->goep, cr);
    ago->goepvdw = ago->goep - ago->goepdih;
  }

  if (PAR(cr) && every(1, step, ago->nstcheckf)) { /* sum ago->fgo */
    /*
    if (sizeof(real) == sizeof(double)) {
      gmx_sumd(ago->nat * 3, (double *) ago->fgo, cr);
    } else {
      gmx_sumf(ago->nat * 3, (float *) ago->fgo, cr);
    }*/
    }

    if (ismaster) {
      tmh_t *tmh = ago->tmh;
      double realtp = 1.0, hdif = 0.0;

      if (ago->mode == 1) { /* adjust lambda */
        ago->lambda += (ago->goep - ago->goeptarget)*ago->lambdadt;
      } else if (ago->mode == 2) {
        hdif = tmh_hdif(tmh, ago->goep, tmh->ec);
        if (ago->tmh_entropic) { /* entropic move */
          tmh_ezmoves(tmh, ago->goep, 1.0);
          ago->tmh_dhde = tmh_getdhde(tmh, ago->goep);
          tmh->tp = tmh->tp0; /* use tp0 */
        } else { /* tempering move */
          tmh_ezmove(tmh, ago->goep, 1.0, ago->tmh_lgvdt);
          ago->tmh_dhde = tmh_getdhde(tmh, ago->goep);
        }

        if (tmh->dtp > 0) {
          ago->lambda = -1.0 / ago->tmh_tpc + ago->tmh_dhde / tmh->tp;
        } else {
          ago->lambda = -ago->tmh_tpc + ago->tmh_dhde * tmh->tp;
        }
      }
      ago->lambda = dblconfine(ago->lambda, ago->lambda_min, ago->lambda_max);

      if (ago->mode < 2) {
        /* collect simple statistical data */
        av_add(ago->avgoep, ago->goep);
        av_add(ago->avrmsd, ago->rmsd);
        av_add(ago->avlamb, ago->lambda);
      }
      if (every(1, step, ago->nstgorep) || bLastStep || bFirstStep) {
        if (ago->mode < 2) {
          double rmsdav = av_getave(ago->avrmsd);
          double goepav = av_getave(ago->avgoep);
          double lambav = av_getave(ago->avlamb);
          printf(llfmt ": goepav %g/%g, RMSDav %g/%g A, lambda %g/%g\n",
            step, ago->goep, goepav, ago->rmsd*10.0, rmsdav*10.0, ago->lambda, lambav);
        } else if (ago->mode == 2) {

          if (tmh->dtp > 0) realtp = tmh->tp; else realtp = 1.0/tmh->tp;
          if (ago->tmh_entropic) {
            printf(llfmt ": goep %9.3f dh %9.3f (iec %2d/%2d), RMSD %7.4f, amp %8.6f, %6.2f%%, %2d, "
                "lambda %8.3f, dhde %7.4f %c tmh->tp %8.5f = %8.5f\n",
              step, ago->goep, hdif,
              tmh->iec, tmh->ergn, ago->rmsd*10.0,
              tmh->wl->lnf, tmh->wl->perc*100.0, tmh->wl->stage,
              ago->lambda, ago->tmh_dhde, (tmh->dtp > 0 ? '/' : '*'), tmh->tp, ago->tmh_dhde/realtp);
          } else {
            printf(llfmt ": goep %9.3f%+9.3f = %9.3f/%9.3f (iec %2d/%2d), RMSD %7.4f, amp %8.6f, %6.2f%%, %2d, "
                "lambda %8.3f, dhde %7.4f %c tmh->tp %8.5f = %8.5f (itp %d/%d)\n",
              step, ago->goep, -tmh->ec, ago->goep - tmh->ec, hdif,
              tmh->iec, tmh->ergn, ago->rmsd*10.0,
              tmh->wl->lnf, tmh->wl->perc*100.0, tmh->wl->stage,
              ago->lambda, ago->tmh_dhde, (tmh->dtp > 0 ? '/' : '*'), tmh->tp, ago->tmh_dhde/realtp, tmh->itp, tmh->tpn);
          }
          if (ago->log) {
            log_printf(ago->log, llfmt " %10.4f %10.4f %10.4f %10.7f %12.10f "
              "%10.8f %3d %10.4f %3d %10.4f %9.6f %6.2f%% %3d\n",
              step, ago->goep, ago->goepdihbias, ago->lambda, ago->tmh_dhde, tmh->wl->lnf,
              tmh->tp, tmh->itp, tmh->ec, tmh->iec, hdif, ago->rmsd*10.0,
              tmh->wl->perc*100, tmh->wl->stage);
          }
        }
      }

      if (every(1, step, ago->nstcheckf)) /* check force */
        agox_checkf(ago, f);

      if (ago->mode == 2 && (every(!bFirstStep, step, ago->tmh_nstsave) || bLastStep) ) {
        tmh_save(tmh, ago->tmh_fntp, ago->tmh_fnehis, ago->tmh_fndhde,
            ago->tmh->wl->lnf, (double) step);
        mtsave(ago->fnrng);
      }

      if (bLastStep) {
        if (ago->log) { printf("closing log\n"); log_close(ago->log); ago->log = NULL; }
      }
    }

    if (PAR(cr)) {
      if (ago->mode > 0) {
        gmx_bcast(sizeof(ago->lambda), &ago->lambda, cr);
      }
    }

    return 0;
}

/* to be used as a replacement of opt2fn(),
 * it will replace the file extension from .mdp to .cfg */
char *agox_opt2fn(const char *opt, int nfile, const t_filenm fnm[])
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
}

#include "txtdump.h"
#include "bondf.h"
#include "qmmm.h"
#include "nonbonded.h"
#include "chargegroup.h"

/* compute local force for a PP node
 * correpsonds to mdlib/force.c/do_force_lowlevel()
 * 1. Wall force
 * 2. Bonded
 * 3. nonbonded
 * 4. electrostatic correction (PP-PME node)
 * 5. PME (PP-PME node)
 *
 * + In 4. and 5. forces are saved to fr->f_novirsum
 * + 4. and 5. can be replaced by their RF counterparts
 * + We assume flags=GMX_FORCE_ALLFORCES,
 *
 * scale[] is equal to ago->scale (no offset)
 */
void agox_forcelow(FILE       *fplog,   gmx_large_int_t step,
                       t_forcerec *fr,      t_inputrec *ir,
                       t_idef     *idef,    t_commrec  *cr,
                       t_nrnb     *nrnb,    gmx_wallcycle_t wcycle,
                       t_mdatoms  *md,
                       t_grpopts  *opts,
                       rvec       x[],      history_t  *hist,
                       rvec       f[],
                       gmx_enerdata_t *enerd,
                       t_fcdata   *fcd,
                       gmx_mtop_t     *mtop,
                       gmx_localtop_t *top,
                       gmx_genborn_t *born,
                       t_atomtypes *atype,
                       gmx_bool       bBornRadii,
                       matrix     box,
                       real       lambda,
                       t_graph    *graph,
                       t_blocka   *excl,
                       rvec       mu_tot[],
                       int        flags,
                       float      *cycles_pme, ago_t *ago)
{
  int     i,status;
  int     donb_flags;
  gmx_bool    bSepDVDL,bSB;
  int     pme_flags;
  matrix  boxs;
  rvec    box_size;
  real    dvdlambda,Vsr,Vlr,Vcorr=0,vdip,vcharge;
  t_pbc   pbc;
  char    buf[22];
  gmx_enerdata_t ed_lam;
  double  lam_i;
  real    dvdl_dum;

#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,sepdvdlformat,s,v,dvdl);
  (void) opts;
  GMX_MPE_LOG(ev_force_start);
  set_pbc(&pbc,fr->ePBC,box);

  /* Reset box */
  for(i=0; (i<DIM); i++)
  {
    box_size[i]=box[i][i];
  }

  bSepDVDL=(fr->bSepDVDL && do_per_step(step,ir->nstlog));

  /* do QMMM first if requested */
  if(fr->bQMMM)
  {
    enerd->term[F_EQM] = calculate_QMMM(cr,x,f,fr,md);
  }

  if (bSepDVDL)
  {
    fprintf(fplog,"Step %s: non-bonded V and dVdl for node %d:\n",
            gmx_step_str(step,buf),cr->nodeid);
  }

  //agox_forcecheck(ago, x, cr, step); getchar();
  agox_goforce(ago, x, f, (real) ago->lambda, cr, step);

  /* Call the short range functions all in one go. */
  GMX_MPE_LOG(ev_do_fnbf_start);

  dvdlambda = 0;


  if (ir->nwall)
  {
    dvdlambda = do_walls(ir,fr,box,md,x,f,lambda,
                         enerd->grpp.ener[egLJSR],nrnb);
    PRINT_SEPDVDL("Walls",0.0,dvdlambda);
    enerd->dvdl_lin += dvdlambda;
  }
  /* If doing GB, reset dvda and calculate the Born radii */
  if (ir->implicit_solvent)
  {
    /* wallcycle_start(wcycle,ewcGB); */

    for(i=0;i<born->nr;i++)
    {
      fr->dvda[i]=0;
    }

    if(bBornRadii)
    {
      calc_gb_rad(cr,fr,ir,top,atype,x,&(fr->gblist),born,md,nrnb);
    }

    /* wallcycle_stop(wcycle, ewcGB); */
  }

  donb_flags = 0;
  if (flags & GMX_FORCE_FORCES)
  {
    donb_flags |= GMX_DONB_FORCES;
  }
  do_nonbonded(cr,fr,x,f,md,excl,
               fr->bBHAM ?
               enerd->grpp.ener[egBHAMSR] :
               enerd->grpp.ener[egLJSR],
               enerd->grpp.ener[egCOULSR],
               enerd->grpp.ener[egGB],box_size,nrnb,
               lambda,&dvdlambda,-1,-1,donb_flags);

  /* If we do foreign lambda and we have soft-core interactions
   * we have to recalculate the (non-linear) energies contributions.
   */
  if (ir->n_flambda > 0 && (flags & GMX_FORCE_DHDL) && ir->sc_alpha != 0)
  {
    init_enerdata(mtop->groups.grps[egcENER].nr,ir->n_flambda,&ed_lam);

    for(i=0; i<enerd->n_lambda; i++)
    {
      lam_i = (i==0 ? lambda : ir->flambda[i-1]);
      dvdl_dum = 0;
      reset_enerdata(&ir->opts,fr,TRUE,&ed_lam,FALSE);
      do_nonbonded(cr,fr,x,f,md,excl,
                   fr->bBHAM ?
                   ed_lam.grpp.ener[egBHAMSR] :
                   ed_lam.grpp.ener[egLJSR],
                   ed_lam.grpp.ener[egCOULSR],
                   enerd->grpp.ener[egGB], box_size,nrnb,
                   (real) lam_i,&dvdl_dum,-1,-1,
                   GMX_DONB_FOREIGNLAMBDA);
      sum_epot(&ir->opts,&ed_lam);
      enerd->enerpart_lambda[i] += ed_lam.term[F_EPOT];
    }
    destroy_enerdata(&ed_lam);
  }

  /* If we are doing GB, calculate bonded forces and apply corrections
   * to the solvation forces */
  if (ir->implicit_solvent)  {
    calc_gb_forces(cr,md,born,top,atype,x,f,fr,idef,
                   ir->gb_algorithm,ir->sa_algorithm,nrnb,bBornRadii,&pbc,graph,enerd);
  }

  if (ir->sc_alpha != 0)
  {
    enerd->dvdl_nonlin += dvdlambda;
  }
  else
  {
    enerd->dvdl_lin    += dvdlambda;
  }
  Vsr = 0;
  if (bSepDVDL)
  {
    for(i=0; i<enerd->grpp.nener; i++)
    {
      Vsr +=
          (fr->bBHAM ?
           enerd->grpp.ener[egBHAMSR][i] :
           enerd->grpp.ener[egLJSR][i])
          + enerd->grpp.ener[egCOULSR][i] + enerd->grpp.ener[egGB][i];
    }
  }
  PRINT_SEPDVDL("VdW and Coulomb SR particle-p.",Vsr,dvdlambda);

  GMX_MPE_LOG(ev_do_fnbf_finish);

  /* Shift the coordinates. Must be done before bonded forces and PPPM,
   * but is also necessary for SHAKE and update, therefore it can NOT
   * go when no bonded forces have to be evaluated.
   */

  /* Here sometimes we would not need to shift with NBFonly,
   * but we do so anyhow for consistency of the returned coordinates.
   */
  if (graph)
  {
    shift_self(graph,box,x);
    if (TRICLINIC(box))
    {
      inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
    }
    else
    {
      inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
    }
  }
  /* Check whether we need to do bondeds or correct for exclusions */
  if (fr->bMolPBC &&
      ((flags & GMX_FORCE_BONDED)
       || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)))
  {
    /* Since all atoms are in the rectangular or triclinic unit-cell,
     * only single box vector shifts (2 in x) are required.
     */
    set_pbc_dd(&pbc,fr->ePBC,cr->dd,TRUE,box);
  }

  if (flags & GMX_FORCE_BONDED)
  {
    GMX_MPE_LOG(ev_calc_bonds_start);
    calc_bonds(fplog,cr->ms,
               idef,x,hist,f,fr,&pbc,graph,enerd,nrnb,lambda,md,fcd,
               DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL, atype, born,
               fr->bSepDVDL && do_per_step(step,ir->nstlog),step);

    /* Check if we have to determine energy differences
     * at foreign lambda's.
     */
    if (ir->n_flambda > 0 && (flags & GMX_FORCE_DHDL) &&
        idef->ilsort != ilsortNO_FE)
    {
      if (idef->ilsort != ilsortFE_SORTED)
      {
        gmx_incons("The bonded interactions are not sorted for free energy");
      }
      init_enerdata(mtop->groups.grps[egcENER].nr,ir->n_flambda,&ed_lam);

      for(i=0; i<enerd->n_lambda; i++)
      {
        lam_i = (i==0 ? lambda : ir->flambda[i-1]);
        dvdl_dum = 0;
        reset_enerdata(&ir->opts,fr,TRUE,&ed_lam,FALSE);
        calc_bonds_lambda(fplog,
                          idef,x,fr,&pbc,graph,&ed_lam,nrnb,(real) lam_i,md,
                          fcd,
                          DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL);
        sum_epot(&ir->opts,&ed_lam);
        enerd->enerpart_lambda[i] += ed_lam.term[F_EPOT];
      }
      destroy_enerdata(&ed_lam);
    }
    GMX_MPE_LOG(ev_calc_bonds_finish);
  }


  *cycles_pme = 0;
  if (EEL_FULL(fr->eeltype))
  {
    bSB = (ir->nwall == 2);
    if (bSB)
    {
      copy_mat(box,boxs);
      svmul(ir->wall_ewald_zfac,boxs[ZZ],boxs[ZZ]);
      box_size[ZZ] *= ir->wall_ewald_zfac;
    }

    clear_mat(fr->vir_el_recip);

    if (fr->bEwald)
    {
      if (fr->n_tpi == 0)
      {
        dvdlambda = 0;
        Vcorr = ewald_LRcorrection(fplog,md->start,md->start+md->homenr,
                                   cr,fr,
                                   md->chargeA,
                                   md->nChargePerturbed ? md->chargeB : NULL,
                                   excl,x,bSB ? boxs : box,mu_tot,
                                   ir->ewald_geometry,
                                   ir->epsilon_surface,
                                   lambda,&dvdlambda,&vdip,&vcharge);
        PRINT_SEPDVDL("Ewald excl./charge/dip. corr.",Vcorr,dvdlambda);
        enerd->dvdl_lin += dvdlambda;
      }
      else
      {
        if (ir->ewald_geometry != eewg3D || ir->epsilon_surface != 0)
        {
          gmx_fatal(FARGS,"TPI with PME currently only works in a 3D geometry with tin-foil boundary conditions");
        }
        /* The TPI molecule does not have exclusions with the rest
         * of the system and no intra-molecular PME grid contributions
         * will be calculated in gmx_pme_calc_energy.
         */
        Vcorr = 0;
      }
    }
    else
    {
      Vcorr = shift_LRcorrection(fplog,md->start,md->homenr,cr,fr,
                                 md->chargeA,excl,x,TRUE,box,
                                 fr->vir_el_recip);
    }

    dvdlambda = 0;
    status = 0;
    switch (fr->eeltype)
    {
    case eelPPPM:
        status = gmx_pppm_do(fplog,fr->pmedata,FALSE,x,fr->f_novirsum,
                             md->chargeA,
                             box_size,fr->phi,cr,md->start,md->homenr,
                             nrnb,ir->pme_order,&Vlr);
        break;
    case eelPME:
    case eelPMESWITCH:
    case eelPMEUSER:
    case eelPMEUSERSWITCH:
        if (cr->duty & DUTY_PME)
        {
          if (fr->n_tpi == 0 || (flags & GMX_FORCE_STATECHANGED))
          {
            pme_flags = GMX_PME_SPREAD_Q | GMX_PME_SOLVE;
            if (flags & GMX_FORCE_FORCES)
            {
              pme_flags |= GMX_PME_CALC_F;
            }
            if (flags & GMX_FORCE_VIRIAL)
            {
              pme_flags |= GMX_PME_CALC_ENER_VIR;
            }
            if (fr->n_tpi > 0)
            {
              /* We don't calculate f, but we do want the potential */
              pme_flags |= GMX_PME_CALC_POT;
            }
            wallcycle_start(wcycle,ewcPMEMESH);
            status = gmx_pme_do(fr->pmedata,
                                md->start,md->homenr - fr->n_tpi,
                                x,fr->f_novirsum,
                                md->chargeA,md->chargeB,
                                bSB ? boxs : box,cr,
                                DOMAINDECOMP(cr) ? dd_pme_maxshift_x(cr->dd) : 0,
                                DOMAINDECOMP(cr) ? dd_pme_maxshift_y(cr->dd) : 0,
                                nrnb,wcycle,
                                fr->vir_el_recip,fr->ewaldcoeff,
                                &Vlr,lambda,&dvdlambda,
                                pme_flags);
            *cycles_pme = (float) wallcycle_stop(wcycle,ewcPMEMESH);

            /* We should try to do as little computation after
             * this as possible, because parallel PME synchronizes
             * the nodes, so we want all load imbalance of the rest
             * of the force calculation to be before the PME call.
             * DD load balancing is done on the whole time of
             * the force call (without PME).
             */
          }
          if (fr->n_tpi > 0)
          {
            /* Determine the PME grid energy of the test molecule
             * with the PME grid potential of the other charges.
             */
            gmx_pme_calc_energy(fr->pmedata,fr->n_tpi,
                                x + md->homenr - fr->n_tpi,
                                md->chargeA + md->homenr - fr->n_tpi,
                                &Vlr);
          }
          PRINT_SEPDVDL("PME mesh",Vlr,dvdlambda);
        }
        else
        {
          /* Energies and virial are obtained later from the PME nodes */
          /* but values have to be zeroed out here */
          Vlr=0.0;
        }
        break;
    case eelEWALD:
        Vlr = do_ewald(fplog,FALSE,ir,x,fr->f_novirsum,
                       md->chargeA,md->chargeB,
                       box_size,cr,md->homenr,
                       fr->vir_el_recip,fr->ewaldcoeff,
                       lambda,&dvdlambda,fr->ewald_table);
        PRINT_SEPDVDL("Ewald long-range",Vlr,dvdlambda);
        break;
    default:
        Vlr = 0;
        gmx_fatal(FARGS,"No such electrostatics method implemented %s",
                  eel_names[fr->eeltype]);
    }
    if (status != 0)
    {
      gmx_fatal(FARGS,"Error %d in long range electrostatics routine %s",
                status,EELTYPE(fr->eeltype));
    }
    enerd->dvdl_lin += dvdlambda;
    enerd->term[F_COUL_RECIP] = Vlr + Vcorr;
  }
  else
  {
    if (EEL_RF(fr->eeltype))
    {
      dvdlambda = 0;

      if (fr->eeltype != eelRF_NEC)
      {
        enerd->term[F_RF_EXCL] =
            RF_excl_correction(fplog,fr,graph,md,excl,x,f,
                               fr->fshift,&pbc,lambda,&dvdlambda);
      }

      enerd->dvdl_lin += dvdlambda;
      PRINT_SEPDVDL("RF exclusion correction",
                    enerd->term[F_RF_EXCL],dvdlambda);
    }
  }

#ifdef GMX_MPI
#endif

  GMX_MPE_LOG(ev_force_finish);

}

#ifdef PRINT_SEPDVDL
#undef PRINT_SEPDVDL
#endif
#define PRINT_SEPDVDL(yn,msg,V,dVdl) if(fplog && yn) fprintf(fplog,"%s: dVdl=%g(this) %g(total)\n",msg,dVdl,V);

/* calculate global force
 * the energies, however, are still local
 * always assumes bDoForces = TRUE,bStateChanged=TRUE,
 * as it is called in md(), bFillGrid==bNS
 * */
static void sum_forces(int start,int end,rvec f[],rvec flr[])
{
  int i;

  if (gmx_debug_at) {
    pr_rvecs(debug,0,"fsr",f+start,end-start);
    pr_rvecs(debug,0,"flr",flr+start,end-start);
  }
  for(i=start; (i<end); i++)
    rvec_inc(f[i],flr[i]);
}

static void calc_f_el(FILE *fp,int  start,int homenr,
                      real charge[],rvec x[],rvec f[],
                      t_cosines Ex[],t_cosines Et[],double t)
{
  rvec Ext;
  real t0;
  int  i,m;

  (void) x;
  for(m=0; (m<DIM); m++)
  {
    if (Et[m].n > 0)
    {
      if (Et[m].n == 3)
      {
        t0 = Et[m].a[1];
        Ext[m] = (real)( cos(Et[m].a[0]*(t-t0))*exp(-dblsqr(t-t0)/(2.0*sqr(Et[m].a[2]))) );
      }
      else
      {
        Ext[m] = (real) cos(Et[m].a[0]*t);
      }
    }
    else
    {
      Ext[m] = 1.0;
    }
    if (Ex[m].n > 0)
    {
      /* Convert the field strength from V/nm to MD-units */
      Ext[m] *= (real)( Ex[m].a[0]*FIELDFAC );
      for(i=start; (i<start+homenr); i++)
          f[i][m] += charge[i]*Ext[m];
    }
    else
    {
      Ext[m] = 0;
    }
  }
  if (fp != NULL)
  {
    fprintf(fp,"%10g  %10g  %10g  %10g #FIELD\n",t,
            Ext[XX]/FIELDFAC,Ext[YY]/FIELDFAC,Ext[ZZ]/FIELDFAC);
  }
}

static void calc_virial(FILE *fplog,int start,int homenr,rvec x[],rvec f[],
            tensor vir_part,t_graph *graph,matrix box,
            t_nrnb *nrnb,const t_forcerec *fr,int ePBC)
{
  int i;

  /* The short-range virial from surrounding boxes */
  clear_mat(vir_part);
  calc_vir(fplog,SHIFTS,fr->shift_vec,fr->fshift,vir_part,ePBC==epbcSCREW,box);
  inc_nrnb(nrnb,eNR_VIRIAL,SHIFTS);

  /* Calculate partial virial, for local atoms only, based on short range.
   * Total virial is computed in global_stat, called from do_md
   */
  f_calc_vir(fplog,start,start+homenr,x,f,vir_part,graph,box);
  inc_nrnb(nrnb,eNR_VIRIAL,homenr);

  /* Add position restraint contribution */
  for(i=0; i<DIM; i++) {
    vir_part[i][i] += fr->vir_diag_posres[i];
  }

  /* Add wall contribution */
  for(i=0; i<DIM; i++) {
    vir_part[i][ZZ] += (real) fr->vir_wall_z[i];
  }

}

static void print_large_forces(FILE *fp,t_mdatoms *md,t_commrec *cr,
                   gmx_large_int_t step,real pforce,rvec *x,rvec *f)
{
  int  i;
  real pf2,fn2;
  char buf[STEPSTRSIZE];

  pf2 = sqr(pforce);
  for(i=md->start; i<md->start+md->homenr; i++) {
    fn2 = norm2(f[i]);
    /* We also catch NAN, if the compiler does not optimize this away. */
    if (fn2 >= pf2 || fn2 != fn2) {
      fprintf(fp,"step %s  atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
          gmx_step_str(step,buf),
          ddglatnr(cr->dd,i),x[i][XX],x[i][YY],x[i][ZZ],sqrt(fn2));
    }
  }
}

void agox_doforce(FILE *fplog,t_commrec *cr,
              t_inputrec *inputrec,
              gmx_large_int_t step,t_nrnb *nrnb,gmx_wallcycle_t wcycle,
              gmx_localtop_t *top,
              gmx_mtop_t *mtop,
              gmx_groups_t *groups,
              matrix box,rvec x[],history_t *hist,
              rvec f[],
              tensor vir_force,
              t_mdatoms *mdatoms,
              gmx_enerdata_t *enerd,t_fcdata *fcd,
              real lambda,t_graph *graph,
              t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
              double t,FILE *field,gmx_edsam_t ed,
              gmx_bool bBornRadii,
              int flags, ago_t *ago)
{
  int    cg0,cg1,i,j;
  int    start,homenr;
  double mu[2*DIM];
  gmx_bool   bSepDVDL,bStateChanged,bNS,bFillGrid,bCalcCGCM,bBS;
  gmx_bool   bDoLongRange,bDoForces,bSepLRF;
  matrix boxs;
  real   e,v,dvdl;
  t_pbc  pbc;
  float  cycles_ppdpme,cycles_pme,cycles_seppme,cycles_force;

  start  = mdatoms->start;
  homenr = mdatoms->homenr;

  bSepDVDL = (fr->bSepDVDL && do_per_step(step,inputrec->nstlog));

  clear_mat(vir_force);

  if (PARTDECOMP(cr))
  {
    pd_cg_range(cr,&cg0,&cg1);
  }
  else
  {
    cg0 = 0;
    if (DOMAINDECOMP(cr))
    {
      cg1 = cr->dd->ncg_tot;
    }
    else
    {
      cg1 = top->cgs.nr;
    }
    if (fr->n_tpi > 0)
    {
      cg1--;
    }
  }

  bStateChanged = (flags & GMX_FORCE_STATECHANGED);
  bNS           = (flags & GMX_FORCE_NS) && (fr->bAllvsAll==FALSE);
  bFillGrid     = (bNS && bStateChanged);
  bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
  bDoLongRange  = (fr->bTwinRange && bNS && (flags & GMX_FORCE_DOLR));
  bDoForces     = (flags & GMX_FORCE_FORCES);
  bSepLRF       = (bDoLongRange && bDoForces && (flags & GMX_FORCE_SEPLRF));

  if (bStateChanged)
  {
    update_forcerec(fplog,fr,box);

    /* Calculate total (local) dipole moment in a temporary common array.
     * This makes it possible to sum them over nodes faster.
     */
    calc_mu(start,homenr,
            x,mdatoms->chargeA,mdatoms->chargeB,mdatoms->nChargePerturbed,
            mu,mu+DIM);
  }

if (fr->ePBC != epbcNONE) {
  /* Compute shift vectors every step,
   * because of pressure coupling or box deformation!
   */
  if ((flags & GMX_FORCE_DYNAMICBOX) && bStateChanged)
    calc_shifts(box,fr->shift_vec);

  if (bCalcCGCM) {
    put_charge_groups_in_box(fplog,cg0,cg1,fr->ePBC,box,
                 &(top->cgs),x,fr->cg_cm);
    inc_nrnb(nrnb,eNR_CGCM,homenr);
    inc_nrnb(nrnb,eNR_RESETX,cg1-cg0);
  }
  else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph) {
    unshift_self(graph,box,x);
  }
}
else if (bCalcCGCM) {
  calc_cgcm(fplog,cg0,cg1,&(top->cgs),x,fr->cg_cm);
  inc_nrnb(nrnb,eNR_CGCM,homenr);
}

if (bCalcCGCM) {
  if (PAR(cr)) {
    move_cgcm(fplog,cr,fr->cg_cm);
  }
  if (gmx_debug_at)
    pr_rvecs(debug,0,"cgcm",fr->cg_cm,top->cgs.nr);
}

#ifdef GMX_MPI
if (!(cr->duty & DUTY_PME)) {
  /* Send particle coordinates to the pme nodes.
   * Since this is only implemented for domain decomposition
   * and domain decomposition does not use the graph,
   * we do not need to worry about shifting.
   */

  wallcycle_start(wcycle,ewcPP_PMESENDX);
  GMX_MPE_LOG(ev_send_coordinates_start);

  bBS = (inputrec->nwall == 2);
  if (bBS) {
    copy_mat(box,boxs);
    svmul(inputrec->wall_ewald_zfac,boxs[ZZ],boxs[ZZ]);
  }

  gmx_pme_send_x(cr,bBS ? boxs : box,x,
                 mdatoms->nChargePerturbed,lambda,
                 ( flags & GMX_FORCE_VIRIAL),step);

  GMX_MPE_LOG(ev_send_coordinates_finish);
  wallcycle_stop(wcycle,ewcPP_PMESENDX);
}
#endif /* GMX_MPI */

  /* Communicate coordinates and sum dipole if necessary */
  if (PAR(cr))
  {
    wallcycle_start(wcycle,ewcMOVEX);
    if (DOMAINDECOMP(cr))
    {
      dd_move_x(cr->dd,box,x);
    }
    else
    {
      move_x(fplog,cr,GMX_LEFT,GMX_RIGHT,x,nrnb);
    }
    /* When we don't need the total dipole we sum it in global_stat */
    if (bStateChanged && NEED_MUTOT(*inputrec))
    {
      gmx_sumd(2*DIM,mu,cr);
    }
    wallcycle_stop(wcycle,ewcMOVEX);
  }
  if (bStateChanged)
  {
    for(i=0; i<2; i++)
    {
      for(j=0;j<DIM;j++)
      {
        fr->mu_tot[i][j] = (real) mu[i*DIM + j];
      }
    }
  }
  if (fr->efep == efepNO)
  {
    copy_rvec(fr->mu_tot[0],mu_tot);
  }
  else
  {
    for(j=0; j<DIM; j++)
    {
      mu_tot[j] = (real)(
          (1.0 - lambda)*fr->mu_tot[0][j] + lambda*fr->mu_tot[1][j] );
    }
  }

  /* Reset energies */
  reset_enerdata(&(inputrec->opts),fr,bNS,enerd,MASTER(cr));
  clear_rvecs(SHIFTS,fr->fshift);

  if (bNS)
  {
    wallcycle_start(wcycle,ewcNS);

    if (graph && bStateChanged)
    {
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(fplog,graph,fr->ePBC,box,x);
    }

    /* Reset long range forces if necessary */
    if (fr->bTwinRange)
    {
      /* Reset the (long-range) forces if necessary */
      clear_rvecs(fr->natoms_force_constr,bSepLRF ? fr->f_twin : f);
    }

    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl = 0;
    ns(fplog,fr,x,box,
       groups,&(inputrec->opts),top,mdatoms,
       cr,nrnb,lambda,&dvdl,&enerd->grpp,bFillGrid,
       bDoLongRange,bDoForces,bSepLRF ? fr->f_twin : f);
    if (bSepDVDL)
    {
      fprintf(fplog,sepdvdlformat,"LR non-bonded",0.0,dvdl);
    }
    enerd->dvdl_lin += dvdl;

    wallcycle_stop(wcycle,ewcNS);
  }

  if (inputrec->implicit_solvent && bNS)
  {
    make_gb_nblist(cr,inputrec->gb_algorithm,inputrec->rlist,
                   x,box,fr,&top->idef,graph,fr->born);
  }

  if (DOMAINDECOMP(cr))
  {
    if (!(cr->duty & DUTY_PME))
    {
      wallcycle_start(wcycle,ewcPPDURINGPME);
      dd_force_flop_start(cr->dd,nrnb);
    }
  }

  /* Start the force cycle counter.
   * This counter is stopped in do_forcelow_level.
   * No parallel communication should occur while this counter is running,
   * since that will interfere with the dynamic load balancing.
   */
  wallcycle_start(wcycle,ewcFORCE);

  if (bDoForces)
  {
    /* Reset forces for which the virial is calculated separately:
     * PME/Ewald forces if necessary */
    if (fr->bF_NoVirSum)
    {
      if (flags & GMX_FORCE_VIRIAL)
      {
        fr->f_novirsum = fr->f_novirsum_alloc;
        GMX_BARRIER(cr->mpi_comm_mygroup);
        if (fr->bDomDec)
        {
          clear_rvecs(fr->f_novirsum_n,fr->f_novirsum);
        }
        else
        {
          clear_rvecs(homenr,fr->f_novirsum+start);
        }
        GMX_BARRIER(cr->mpi_comm_mygroup);
      }
      else
      {
        /* We are not calculating the pressure so we do not need
         * a separate array for forces that do not contribute
         * to the pressure.
         */
        fr->f_novirsum = f;
      }
    }

    if (bSepLRF)
    {
      /* Add the long range forces to the short range forces */
      for(i=0; i<fr->natoms_force_constr; i++)
      {
        copy_rvec(fr->f_twin[i],f[i]);
      }
    }
    else if (!(fr->bTwinRange && bNS))
    {
      /* Clear the short-range forces */
      clear_rvecs(fr->natoms_force_constr,f);
    }

    clear_rvec(fr->vir_diag_posres);

    GMX_BARRIER(cr->mpi_comm_mygroup);
  }
  if (inputrec->ePull == epullCONSTRAINT)
  {
    clear_pull_forces(inputrec->pull);
  }

  /* update QMMMrec, if necessary */
  if(fr->bQMMM)
  {
    update_QMMMrec(cr,fr,x,mdatoms,box,top);
  }

  if ((flags & GMX_FORCE_BONDED) && top->idef.il[F_POSRES].nr > 0)
  {
    /* Position restraints always require full pbc */
    set_pbc(&pbc,inputrec->ePBC,box);
    v = posres(top->idef.il[F_POSRES].nr,top->idef.il[F_POSRES].iatoms,
               top->idef.iparams_posres,
               (const rvec*)x,fr->f_novirsum,fr->vir_diag_posres,
               inputrec->ePBC==epbcNONE ? NULL : &pbc,lambda,&dvdl,
               fr->rc_scaling,fr->ePBC,fr->posres_com,fr->posres_comB);
    if (bSepDVDL)
    {
      fprintf(fplog,sepdvdlformat,
              interaction_function[F_POSRES].longname,v,dvdl);
    }
    enerd->term[F_POSRES] += v;
    /* This linear lambda dependence assumption is only correct
     * when only k depends on lambda,
     * not when the reference position depends on lambda.
     * grompp checks for this.
     */
    enerd->dvdl_lin += dvdl;
    inc_nrnb(nrnb,eNR_POSRES,top->idef.il[F_POSRES].nr/2);
  }

  /* Compute the bonded and non-bonded energies and optionally forces */
  agox_forcelow(fplog,step,fr,inputrec,&(top->idef),
                    cr,nrnb,wcycle,mdatoms,&(inputrec->opts),
                    x,hist,f,enerd,fcd,mtop,top,fr->born,
                    &(top->atomtypes),bBornRadii,box,
                    lambda,graph,&(top->excls),fr->mu_tot,
                    flags,&cycles_pme, ago);

  cycles_force = (float) wallcycle_stop(wcycle,ewcFORCE);
  GMX_BARRIER(cr->mpi_comm_mygroup);

  if (ed)
  {
    do_flood(fplog,cr,x,f,ed,box,step);
  }

  if (DOMAINDECOMP(cr))
  {
    dd_force_flop_stop(cr->dd,nrnb);
    if (wcycle)
    {
      dd_cycles_add(cr->dd,cycles_force-cycles_pme,ddCyclF);
    }
  }

  if (bDoForces)
  {
    if (IR_ELEC_FIELD(*inputrec))
    {
      /* Compute forces due to electric field */
      calc_f_el(MASTER(cr) ? field : NULL,
                start,homenr,mdatoms->chargeA,x,fr->f_novirsum,
                inputrec->ex,inputrec->et,t);
    }

    /* Communicate the forces */
    if (PAR(cr))
    {
      wallcycle_start(wcycle,ewcMOVEF);
      if (DOMAINDECOMP(cr))
      {
        dd_move_f(cr->dd,f,fr->fshift);
        /* Do we need to communicate the separate force array
         * for terms that do not contribute to the single sum virial?
         * Position restraints and electric fields do not introduce
         * inter-cg forces, only full electrostatics methods do.
         * When we do not calculate the virial, fr->f_novirsum = f,
         * so we have already communicated these forces.
         */
        if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl &&
            (flags & GMX_FORCE_VIRIAL))
        {
          dd_move_f(cr->dd,fr->f_novirsum,NULL);
        }
        if (bSepLRF)
        {
          /* We should not update the shift forces here,
           * since f_twin is already included in f.
           */
          dd_move_f(cr->dd,fr->f_twin,NULL);
        }
      }
      else
      {
        pd_move_f(cr,f,nrnb);
        if (bSepLRF)
        {
          pd_move_f(cr,fr->f_twin,nrnb);
        }
      }
      wallcycle_stop(wcycle,ewcMOVEF);
    }

    /* If we have NoVirSum forces, but we do not calculate the virial,
     * we sum fr->f_novirum=f later.
     */
    if (vsite && !(fr->bF_NoVirSum && !(flags & GMX_FORCE_VIRIAL)))
    {
      wallcycle_start(wcycle,ewcVSITESPREAD);
      spread_vsite_f(fplog,vsite,x,f,fr->fshift,nrnb,
                     &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
      wallcycle_stop(wcycle,ewcVSITESPREAD);

      if (bSepLRF)
      {
        wallcycle_start(wcycle,ewcVSITESPREAD);
        spread_vsite_f(fplog,vsite,x,fr->f_twin,NULL,
                       nrnb,
                       &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
        wallcycle_stop(wcycle,ewcVSITESPREAD);
      }
    }

    if (flags & GMX_FORCE_VIRIAL)
    {
      /* Calculation of the virial must be done after vsites! */
      calc_virial(fplog,mdatoms->start,mdatoms->homenr,x,f,
                  vir_force,graph,box,nrnb,fr,inputrec->ePBC);
    }
  }

  if (inputrec->ePull == epullUMBRELLA || inputrec->ePull == epullCONST_F)
  {
    /* Calculate the center of mass forces, this requires communication,
     * which is why pull_potential is called close to other communication.
     * The virial contribution is calculated directly,
     * which is why we call pull_potential after calc_virial.
     */
    set_pbc(&pbc,inputrec->ePBC,box);
    dvdl = 0;
    enerd->term[F_COM_PULL] =
        pull_potential(inputrec->ePull,inputrec->pull,mdatoms,&pbc,
                       cr,t,lambda,x,f,vir_force,&dvdl);
    if (bSepDVDL)
    {
      fprintf(fplog,sepdvdlformat,"Com pull",enerd->term[F_COM_PULL],dvdl);
    }
    enerd->dvdl_lin += dvdl;
  }

  if (PAR(cr) && !(cr->duty & DUTY_PME))
  {
    cycles_ppdpme = (float) wallcycle_stop(wcycle,ewcPPDURINGPME);
    dd_cycles_add(cr->dd,cycles_ppdpme,ddCyclPPduringPME);

    /* In case of node-splitting, the PP nodes receive the long-range
     * forces, virial and energy from the PME nodes here.
     */
    wallcycle_start(wcycle,ewcPP_PMEWAITRECVF);
    dvdl = 0;
    gmx_pme_receive_f(cr,fr->f_novirsum,fr->vir_el_recip,&e,&dvdl,
                      &cycles_seppme);
    if (bSepDVDL)
    {
      fprintf(fplog,sepdvdlformat,"PME mesh",e,dvdl);
    }
    enerd->term[F_COUL_RECIP] += e;
    enerd->dvdl_lin += dvdl;
    if (wcycle)
    {
      dd_cycles_add(cr->dd,cycles_seppme,ddCyclPME);
    }
    wallcycle_stop(wcycle,ewcPP_PMEWAITRECVF);
  }

  if (bDoForces && fr->bF_NoVirSum)
  {
    if (vsite)
    {
      /* Spread the mesh force on virtual sites to the other particles...
       * This is parallellized. MPI communication is performed
       * if the constructing atoms aren't local.
       */
      wallcycle_start(wcycle,ewcVSITESPREAD);
      spread_vsite_f(fplog,vsite,x,fr->f_novirsum,NULL,nrnb,
                     &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
      wallcycle_stop(wcycle,ewcVSITESPREAD);
    }
    if (flags & GMX_FORCE_VIRIAL)
    {
      /* Now add the forces, this is local */
      if (fr->bDomDec)
      {
        sum_forces(0,fr->f_novirsum_n,f,fr->f_novirsum);
      }
      else
      {
        sum_forces(start,start+homenr,f,fr->f_novirsum);
      }
      if (EEL_FULL(fr->eeltype))
      {
        /* Add the mesh contribution to the virial */
        m_add(vir_force,fr->vir_el_recip,vir_force);
      }
    }
  }

  /* scaling code for ago here */
  if (ago != NULL) { }

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(inputrec->opts),enerd);

  if (fr->print_force >= 0 && bDoForces)
  {
    print_large_forces(stderr,mdatoms,cr,step,fr->print_force,x,f);
  }
}


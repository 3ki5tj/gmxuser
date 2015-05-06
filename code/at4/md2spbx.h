/* functions missing in md2spb.h */

/* allocate memory and assign default settings for a single spb
 * We assume that spb->bins, spb->min, spb->max
 * have been set already */
#define spb_initdefault(bs, id, flags) spb_initdefault_(bs, id, flags, __FILE__,__LINE__)
static int spb_initdefault_(spbonds_t *bs, int id, int flags, const char *fnsrc, int lineno)
{
  int j;
  spb_t *spb = bs->arr+id;

  spb->id       = id;
  spb->type     = -1;
  spb->exclcnt  = 0;
  spb->exclcap  = 0;
  spb->excl     = NULL;
  spb->bufcnt   = 0;
  spb->bufcap   = 0;
  spb->buf      = NULL;
  spb->mfcom    = 0.0;
  spb->pmfmin   = 0.0;
  spb->pmfmax   = 0.0;
  spb->flags    = SPB_PERIODIC;

  if (flags > 0) {
    spb->atmbuf = NULL;
    spb->atmsiz = 0;
    for (j = 0; j < SPB_MAXNATOMS; j++) {
      spb->atoms[j] = NULL;
    }
  }

  /* copy the default common settings */
  spb->distmin    = bs->distmin;
  spb->distmax    = bs->distmax;
  spb->distpwrscl = bs->distpwrscl;
  spb->mfmax      = bs->mfmax;
  spb->pmfside    = bs->pmfside;
  for (j = 0; j < SPB_PMFSCAL_CNT; j++)
    spb->popscal[j] = bs->popscal[j];

  /* round spb->min and spb->max to multiples of M_PI, if they are close */
  spb_round2pi(spb);
  spb->binw = (spb->max-spb->min)/spb->bins;

  /* allocate array (double), we add size by one for safety */
  #define SPB_ALLOCARR(var, num) \
    if ((spb->var = calloc((num)+1, sizeof spb->var[0])) == NULL) { \
      die_if(1, "cannot allocate %s for spb %d, from %s, line %d.\n",   \
          #var, id, fnsrc, lineno);                                 \
    } else { for (j = 0; j < (num); j++) spb->var[j] = 0.0;       }

  /* allocate spaces for the distribution */
  SPB_ALLOCARR(distref, spb->bins + 1);

  /* assume a uniform distribution by default */
  for (j = 0; j <= spb->bins; j++)
    spb->distref[j] = 1.0;

  SPB_ALLOCARR(mfref,  spb->bins);

  SPB_ALLOCARR(hist,   spb->bins);
  SPB_ALLOCARR(sbf,    spb->bins);
  SPB_ALLOCARR(sdiv,   spb->bins);
  SPB_ALLOCARR(sbl,    spb->bins);
  SPB_ALLOCARR(mf,     spb->bins);
  SPB_ALLOCARR(pmf,    spb->bins+1);
  SPB_ALLOCARR(dbltmp, spb->bins+1);

  return 0;
}

/* allocate bs->arr and assign default setting
 * we assume bs->cnt is set
 * this is for the master only
 * individual spb information is not initialized here
 * */
#define spbs_initdefault(bs,iopref) spbs_initdefault_(bs,iopref,__FILE__,__LINE__)
static int spbs_initdefault_(spbonds_t *bs, const char *iopref, const char *fnsrc, int lineno)
{
  int i;

  die_if(bs->cnt < 0, "# of spbs is negative.\n", bs->cnt);

  /* for safety, we copy default setting even if bs->cnt == 0; */
  bs->arr         = NULL;
  bs->mpi_rank    = 0;
  bs->distmin     = 1e-4;
  bs->distmax     = 1.0;
  bs->distpwrscl  = 1.0;
  bs->mfmax       = 400.0;
  bs->pmfside     = -1; /* now default is -1, i.e., pmf is always nonpositive,
                         * which is natural in separate peaks */
  for (i = 0; i < SPB_PMFSCAL_CNT; i++)
    bs->popscal[i] = 1.0;
  bs->bin_file = ssdup(""); /* cannot know if file names are previously */
  bs->txt_file = ssdup(""); /* managed by ss, so use ssdup instead of sscpy */
  if (iopref) { /* install the prefix iopref first */
    sscpy(bs->bin_file, iopref);
    sscpy(bs->txt_file, iopref);
  }
  sscat(bs->bin_file, "spb.bin");
  sscat(bs->txt_file, "spb.txt");

  /* skip the allocation step */
  if (bs->cnt == 0) return 0;

  /* allocate space for the spb array */
  die_if((bs->arr = calloc(bs->cnt, sizeof bs->arr[0])) == NULL,
    "cannot allocate space for spb array, file: %s, lineno: %d\n",
    fnsrc, lineno);
  return 0;
}


/* initialize an spbonds_t structure from a previous text output,
 * NOTE!!!
 * it will not check its validity against the current settings. */
SPBSTRCLS int spbs_read(spbonds_t *bs, const char *fname, unsigned flags)
{
  int j, id, econ, ver, cntmax, next;
  int verbose, verify = 0;
  spb_t *spb0, *spb;
  FILE *fp;
  char *buf = NULL, *p, *q;
  size_t size;

  verbose = (flags&SPB_VERBOSE);

  if (verify) {
    fprintf(stderr, "spbs_read: will try to preserve the original infomration\n");
  }

  spbs_check(bs, SPB_CHECKNULL|SPB_CHECKMASTER);

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot read from file [%s].\n", fname);
    return -1;
  }

  /* get global information from the first line */
  if (NULL == ssfgets(buf, &size, fp))
    return -1;

  bs->cnt = 0;
  p = buf+1;
  if (1 != sscanf(p, "%d%n", &bs->cnt, &next) || bs->cnt <= 0) {
    fprintf(stderr, "cannot read cnt, p=%s\n", p);
    goto QUIT;
  }
  if (verbose) fprintf(stderr, "%s has %d spbs\n", fname, bs->cnt);
  p += next;
  if (1 != sscanf(p, "%d%n", &ver, &next) ) {
    if (verbose) fprintf(stderr, "spb warning: a version 0 text input\n");
    ver = 0;
  } else {
    p += next;
  }

  if (!verify) {
    if (0 != spbs_initdefault(bs, NULL)) goto QUIT;
  }

  /* read basic information for each spb */
  spb0 = bs->arr;
  for (econ = 1, cntmax = 0, id = 0; id < bs->cnt; id++) {
    int spbbins;
    double spbmin, spbmax;

    spb = bs->arr+id;
    q = strchr(p, '|');
    if (q == NULL) {
      fprintf(stderr, "corrupted information line %s", buf);
      goto QUIT;
    }
    q++;
    if (3 != sscanf(q, "%d%lf%lf%n", &(spbbins), &(spbmin), &(spbmax), &next)) {
      fprintf(stderr, "insufficient basic information for spb %d\n", id);
      goto QUIT;
    }
    q += next;

    /* for version 1 or higher, we can initialize atom information,
     * so we use a stronger cleaning */
    if (!verify) {
      spb->min = spbmin;
      spb->max = spbmax;
      spb->bins = spbbins;
      if (0 != spb_initdefault(bs, id, (ver > 0))) goto QUIT;
    } else {
      /* keep the original values  */
      if (fabs(spbmin-spb->min) > 1e-3 || fabs(spbmax-spb->max) > 1e-3 || spb->bins != spbbins) {
        fprintf(stderr, "[%d] range is wrong got %d (%g,%g) should be %d (%g,%g)\n",
            id, spbbins, spbmin, spbmax, spb->bins, spb->min, spb->max);
        goto QUIT;
      }
    }

    if (cntmax < spb->bins) cntmax = spb->bins;

    if (ver > 0) { /* read atom names */
      char atbuf[128]=""; /* warning limited buffer */

      if (2 != sscanf(q, "%lf%d%n", &(spb->mfmax), &(spb->pmfside), &next)) {
        fprintf(stderr, "insufficient additional information for spb %d\n", id);
        goto QUIT;
      }
      q += next;

      if (1 != sscanf(q, "%d%n", &(spb->natoms), &next)) {
        fprintf(stderr, "cannot have the # of atoms for spb %d\n", id);
        goto QUIT;
      }
      q += next;

      if (spb->natoms > 0) {
        if (1 != sscanf(q, "%127s%n", atbuf, &next)) {
          fprintf(stderr, "cannot have the atoms for spb %d\n", id);
          goto QUIT;
        }
        q += next;

        spb->atmbuf = ssdup(atbuf);
        /* compute spb->atoms, spb->natoms, spb->atmcap */
        spb_parse_atoms_txt(spb, verbose);
      }
    }

    p = q;

    if (id > 0 && (fabs(spb->min-spb0->min) > 1e-8 || fabs(spb->max-spb0->max) > 1e-8 || spb->bins != spb0->bins) ) {
      econ = 0; /* cannot use economic (multiple-column) mode */
    }
  }
  if (verbose) fprintf(stderr, "reading economic mode=%d, lines=%d\n", econ, cntmax);

  if (ver > 0) cntmax++;
  for (j = 0; j < cntmax; j++) {
    if (ssfgets(buf, &size, fp) == NULL) {
      fprintf(stderr, "unable to read line %d\n", j);
      break;
    }
    p = buf;

    for (id = 0; id < bs->cnt; id++) {
      double arr[9];
      spb = bs->arr+id;
      if (id == 0 || !econ) {
        if (1 != sscanf(p, "%lf%n", &arr[0], &next)) {
          fprintf(stderr, "cannot read index. p=%s\n", p);
          goto QUIT;
        }
        p += next;
      }
      if (7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf%n",
           &arr[1], &arr[2], &arr[3], &arr[4],
           &arr[5], &arr[6], &arr[7], &next) ) {
        fprintf(stderr, "insufficient data!, row=%d, id=%d, p=%s\n", j, id, p);
        goto QUIT;
      }
      p += next;
      if (ver >= 2) {
        if (1 != sscanf(p, "%lf%n", &arr[8], &next)) {
          fprintf(stderr, "insufficient data, row=%d, id=%d, p=%s\n", j, id, p);
          goto QUIT;
        }
        p += next;
      } else {
        arr[8] = arr[6]; /* distref0 = distref */
      }
      if (j < spb->bins) {
        spb->hist[j] = arr[1];
        spb->sbl[j] = arr[4]*spb->hist[j];
        spb->sbf[j] = arr[2]*spb->sbl[j];
        spb->sdiv[j] = arr[3]*spb->sbl[j];
        if (!verify) spb->mfref[j] = arr[7];
      }
      if (j <= spb->bins) {
        spb->pmf[j] = arr[5];
        if (!verify) {
          spb->distref[j] = arr[6];
          spb->distref0[j] = arr[8];
        }
      }
    } /* loop over spbs */
  }/* end of loop over bins (rows) */

  for (id = 0; id < bs->cnt; id++) {
    spb = bs->arr+id;
    //fprintf(stderr, "rt sbf: %25.15e %25.15e...\n", spb->sbf[0], spb->sbf[1]);
    //fprintf(stderr, "rt mf: %25.15e %25.15e...\n", spb->mf[0], spb->mf[1]);
  }
  //fprintf(stderr, "[0] %s... [1] %s...\n", bs->arr[0].atoms[0], bs->arr[1].atoms[0]);

  spbs_calcmfpmf(bs, flags&SPB_VERBOSE);
  fclose(fp);
  return 0;
QUIT:
  fclose(fp);
  return -1;
}



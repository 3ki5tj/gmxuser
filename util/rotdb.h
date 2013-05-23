#define RBSIZ 20  /* number of rotamers an integer can encode */
#define RBMAX 8   /* number of blocks */
/* rotamer database */
typedef struct {
  unsigned code[RBMAX]; /* code */
  double w; /* weight */
} rotrec_t;

typedef struct {
  rotrec_t *arr, ref;
  int len, n;
} rotdb_t;

static void rotcode_copy(unsigned *dest, const unsigned *src, int len)
{
  int j, m;

  m = (len+RBSIZ-1)/RBSIZ;
  for (j = 0; j < m; j++)
    dest[j] = src[j];
}

static int rotcode_cmp(const unsigned *dest, const unsigned *src, int len)
{
  int j, m;

  m = (len+RBSIZ-1)/RBSIZ;
  for (j = 0; j < m; j++) {
    if (dest[j] > src[j]) return 1;
    else if (dest[j] < src[j]) return -1;
  }
  return 0;
}

static void rotcode_print(FILE *fp, const unsigned *c, int len)
{
  int j, m;

  m = (len+RBSIZ-1)/RBSIZ;
  for (j = 0; j < m; j++) {
    fprintf(fp, "%11u ", c[j]);
  }
}

static rotdb_t *rotdb_init(int len, unsigned *code)
{
  rotdb_t *rot;

  xnew(rot, 1);
  rot->n = 0;
  rot->len = len;
  if (len > RBMAX*RBSIZ) {
    fprintf(stderr, "code too long len = %d\n", len);
    exit(1);
  }
  rotcode_copy(rot->ref.code, code, len);
  rot->ref.w = 1;
  xnew(rot->arr, 1);
  return rot;
}

static void rotdb_free(rotdb_t *rot)
{
  if (rot->arr) free(rot->arr);
  free(rot);
}

/* add a rotamer configuration */
static int rotdb_add(rotdb_t *rot, unsigned *code, double w)
{
  int i;
  rotrec_t *rr;

  /* search for existing code */
  for (i = 0; i < rot->n; i++) {
    if (rotcode_cmp(rot->arr[i].code, code, rot->len) == 0)
      break;
  }
  if (i >= rot->n) { /* add a new rotamer */
    rot->n += 1;
    xrenew(rot->arr, rot->n);
    rr = rot->arr + rot->n - 1;
    rotcode_copy(rr->code, code, rot->len);
    rr->w = w;
  } else {
    rr = rot->arr + i;
    rr->w += w;
  }
  return 0;
}

static int rotrec_cmp(const void *a, const void *b)
{
  const rotrec_t *ra = (const rotrec_t *)a, *rb = (const rotrec_t *)b;
  //printf("%g %g %p %p\n", ra->w, rb->w, ra, rb); getchar();
  if (ra->w > rb->w) return -1; /* descending */
  else if (ra->w < rb->w) return 1;
  return 0;
}

static void rotdb_trim(rotdb_t *rot)
{
  qsort(rot->arr, rot->n, sizeof(rot->arr[0]), rotrec_cmp);
}

/* return a rotamer string */
static char *rotdb_decode(const unsigned *code, int len)
{
  int i, ip = 0;
  unsigned icode;
  static char s[RBSIZ*RBMAX+1];

  for (i = 0; i < len; i++) {
    if (i % RBSIZ == 0) {
      icode = code[ip];
      ip++;
    }
    s[i] = (char)('0' + (icode % 3));
    icode /= 3;
  }
  return s;
}

static int rotdb_write(const rotdb_t *rot, const char *fn)
{
  FILE *fp;
  int i;
  rotrec_t *rr;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d %d ", rot->len, rot->n);
  rotcode_print(fp, rot->ref.code, rot->len);
  fprintf(fp, "%s\n", rotdb_decode(rot->ref.code, rot->len));
  for (i = 0; i < rot->n; i++) {
    rr = rot->arr + i;
    rotcode_print(fp, rr->code, rot->len);
    fprintf(fp, "%12.3f\t%s\n", rr->w, rotdb_decode(rr->code, rot->len));
  }
  fclose(fp);
  return 0;
}

/* load previous database */
static int rotdb_load(rotdb_t *rot, const char *fn)
{
  FILE *fp;
  char s[1024], *p;
  int len, n, i, m, j, next;
  unsigned code[RBMAX];
  double wt;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  /* the first line */
  fgets(s, sizeof s, fp);
  if (2 != sscanf(s+1, "%d%d", &len, &n)) {
    fprintf(stderr, "broken first line %s: %s\n", fn, s);
    fclose(fp);
    return -1;
  }
  if (len != rot->len) {
    fprintf(stderr, "fn %s, len %d vs. %d mismatch\n",
        fn, len, rot->len);
    fclose(fp);
    return -1;
  }
  for (i = 0; i < n; i++) {
    fgets(s, sizeof s, fp);
    m = (len+RBSIZ-1)/RBSIZ;
    /* scan the code */
    p = s;
    for (j = 0; j < m; j++) {
      if (1 != sscanf(p, "%u%n", &code[j], &next)) {
        fprintf(stderr, "unable to read code %d at block %d, %s\n", i, j, fn);
        fclose(fp);
        return -1;
      }
      p += next;
    }

    if (1 != sscanf(p, "%lf", &wt)) {
      fprintf(stderr, "unable to read weight i = %d, %s\n", i, fn);
      fclose(fp);
      return -1;
    }
    rotdb_add(rot, code, wt);
  }
  fclose(fp);
  return 0;
}


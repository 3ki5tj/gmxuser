/*
  Submit multiple successive PBS files,
  For SUG@R, BIOU clusters in Rice University
  Be careful, the program is powerful but could be dangerous.

    Copyright (C) 2010  Cheng Zhang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 + check existing jobs with the same name
 + automatically append cluster name

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define SMTFILE "PBS_SUBMIT_RESULT"
int verbose = 1;
char cmd[1024]; /* command line */
char buf[1024]; /* line buffer */
char cluster[1024]; /* cluster name */

#define TRUE 1
#define FALSE 0

/* execute command `s' */
static int docmd(const char *s, int id, int silent)
{
  if (!silent) {
    printf("CMD");
    if (id > 0) printf(" %2d", id);
    printf(": %s\n", s);
  }
  return system(s);
}

static int check_after(char *after)
{
  char *p;
  // trim the tail
  for (p = after+strlen(after)-1; p >= after && isspace(*p); p--) {
    *p='\0';
  }
  for (p = after; *p && isspace(*p); p++) ;

  if (*p=='\0') return 1;
  for (; *p && *p !='.'; p++) {
    if (!isdigit(*p)) return 1;
  }
  return 0;
}

static int get_clustername(void)
{
  FILE *fp;
  char *p;

  sprintf(cmd, "qstat -u `whoami` > %s", SMTFILE);
  docmd(cmd, 0, TRUE);

  if ( (fp = fopen(SMTFILE, "r")) == NULL) {
    fprintf(stderr, "Cannot read output file [%s]\n", SMTFILE);
  } else {
    while (fgets(buf, sizeof buf, fp)) {
      if ( (p = strchr(buf, ':')) != NULL) {
        *p='\0';
        strcpy(cluster, buf);
        printf("cluster is [%s]\n", cluster);
        break;
      }
    }
    fclose(fp);
  }
  remove(SMTFILE); // remove the output if it exists
  return 0;
}

// scan the queue for jobs with the same name as `jobname'
// return the latest job in the queue with the name `jobname'
// or -1 if none
static int find_existing(const char *jobname, int after, int *match)
{
  int last = -1;
  FILE *fp;

  sprintf(cmd, "qstat -u `whoami` | grep %s > %s", jobname, SMTFILE);
  docmd(cmd, 0, FALSE);

  if (match != NULL) *match = 0;

  // read the output
  if ( (fp = fopen(SMTFILE, "r")) == NULL) {
    fprintf(stderr, "Cannot read output file [%s]\n", SMTFILE);
    return 1;
  }
  buf[0]='\0';
  while( fgets(buf, sizeof buf, fp) != NULL) { // not empty
    static char jstr[80], usrnm[80], qnm[80], jnm[80];
    int jid = 0;

    sscanf(buf, "%s%s%s%s", jstr, usrnm, qnm, jnm);
    if (strcmp(jnm, jobname) == 0) {
      jid = atoi(jstr);
      if (match != NULL && jid == after) *match = 1;
      if (jid > last) last = jid;
      printf("%s", buf);
    }
  }
  fclose(fp);
  remove(SMTFILE);
  return last;
}

static int intcmp(const void*p, const void*q)
{
  int i, j;
  i = *((const int *) p);
  j = *((const int *) q);
  return i - j;
}

// delete job of the name `jobname'
static int delete_job(const char *jobname, int delete_running)
{
  int *jobs, lines, jcnt, jid, j;
  FILE *fp;

  get_clustername();

  sprintf(cmd, "qstat -u `whoami` | grep %s > %s", jobname, SMTFILE);
  docmd(cmd, 0, FALSE);

  // read the output
  if ((fp = fopen(SMTFILE, "r")) == NULL) {
    fprintf(stderr, "Cannot read output file [%s]\n", SMTFILE);
    return 1;
  }
  buf[0]='\0';
  // first count the # of lines
  for (lines = 0; fgets(buf, sizeof buf, fp) != NULL; ) {
    lines++;
  }
  if (lines <= 0) {
    fprintf(stderr, "no jobs were found\n");
    fclose(fp);
    remove(SMTFILE);
    return 1;
  }

  if ((jobs = calloc(lines, sizeof(int))) == NULL) {
    fprintf(stderr, "cannot allocate spaces for %d jobs\n", lines);
    return -1;
  }

  fseek(fp, 0L, SEEK_SET); // go back to the beginning
  for (jcnt = 0; fgets(buf, sizeof buf, fp) != NULL; ) { // not empty
    char jstr[80], usrnm[80], qnm[80], jnm[80];

    sscanf(buf, "%s%s%s%s", jstr, usrnm, qnm, jnm);
    // make sure the jobname is right,
    // although we did a grep in the command line, it doesn't mean it
    // occurs in the right column
    if (strcmp(jnm, jobname) == 0) {

      jid = atoi(jstr);
      if ('R' == buf[86]) {
        if (!delete_running) continue;
      }
      // put job id's into the array
      jobs[jcnt++] = jid;
      printf("%s", buf);
    }
  }
  fclose(fp);
  remove(SMTFILE);

  // sort jobs
  qsort(jobs, jcnt, sizeof(int), intcmp);
  // remove jobs from the queue in the reverse order
  for (j = jcnt-1; j >= 0; j--) {
    sprintf(cmd, "qdel %d.%s", jobs[j], cluster);
    docmd(cmd, j+1, FALSE);
  }

  return 0;
}

int main(int argc, char *argv[])
{
  char *fname;
  int nrep, i;
  FILE *fp;
  static char option[1024], after[1024], jobname[1024], *p;

  // check if we want to delete a job, option -d or -D
  if (argc >= 2 && argv[1][0]=='-' && strchr("dD", argv[1][1]) != NULL) {
    char *jname = NULL;
    if (argv[1][2]=='\0') {
      if (argc >= 3) jname = argv[2];
    } else {
      jname = argv[1]+2;
    }
    if (jname != NULL) {
      delete_job(jname, (argv[1][1]=='D'));
    } else {
      goto HELP;
    }
    return 0;
  }

  if (argc < 3) {
HELP:
    fprintf(stderr, "%s  Copyright (C) 2010  Cheng Zhang\n"
      "This program comes with ABSOLUTELY NO WARRANTY. "
      "It is free software, and you are welcome to redistribute "
      "it under certain conditions.\n\n", argv[0]);

    fprintf(stderr,
        "Example:\n\n"
        "  %s 3 foo.pbs\n"
        "to submit three successive foo.pbs, or\n\n"
        "  %s 10 abc.pbs 1088.sugarman.rcsg.rice.edu\n"
        "to subimt ten abc.pbs after job 1088\n\n", argv[0], argv[0]);
    fprintf(stderr,
        "To delete jobs from queue\n"
        "  %s -d jobname\n"
        "Replace -d by -D if you want to delete currently-running jobs\n",
        argv[0]);
    return 1;
  }


  for (p = argv[1]; *p; p++) {
    if (!isdigit(*p)) goto BAD_NREP;
  }

  nrep = atoi(argv[1]);
  if (nrep <= 0 || nrep > 1000) {
BAD_NREP:
    fprintf(stderr, "# of repeats [%s] is crazy.\n", argv[1]);
    return 1;
  }

  fname = argv[2];
  // test if the file exists
  if ( (fp = fopen(fname, "r")) == NULL ) {
    fprintf(stderr, "It seems the PBS file [%s] does not exist.\n", fname);
    return 1;
  } else {
    // try to get the jobname
    for (jobname[0]='\0'; fgets(buf, sizeof buf, fp); ) {
      if (strncmp(buf, "#PBS", 4) != 0) continue;
      if ( (p = strstr(buf, "-N")) == NULL) continue;
      p += 3;
      sscanf(p, "%s", jobname);
      printf("Job name is [%s]\n", jobname);
    }
    fclose(fp);
  }
  if (jobname[0] == '\0') {
    fprintf(stderr, "cannot identify job name.\n");
    return 1;
  }

  if (argc > 3) strncpy(after, argv[3], sizeof after);

  // try to get the cluster name
  get_clustername();

  // detect job clashing
  if (after[0] == '\0') {
    int jid_last = find_existing(jobname, -1, NULL);

    if (jid_last >= 0) {
      printf("\nAt least one job with the same or similar name is already running.\n");
      printf("You might want to use\n%s %d %s %d",
        argv[0], nrep, fname, jid_last);
      if (cluster[0]) printf(".%s", cluster);
      printf("\nBut DO check if %d is the latest by\nqstat -f %d | grep beforeok\n", jid_last, jid_last);
      printf("Do you want to continue (y/n)?\n");
      fgets(buf, sizeof buf, stdin);
      if (buf[0]!='y' && buf[0]!='Y') {
        printf("Abort submission.\n");
        return 1;
      }
    }

  } else { // an after-job is specified
    int jid_after, jid_last, match;

    jid_after = (int) atoi(after);
    jid_last = find_existing(jobname, jid_after, &match);

    // trying to check if the job to follow belongs to the
    // current project and if it is the last job.
    if (!match || (match && (jid_last != jid_after))) {
      if (match) {
        printf("Not following the last job (%d) of the project\n", jid_last);
      } else {
        printf("Trying to follow a job of a different or nonexisting project, currrent last=%d\n", jid_last);
      }

      printf("Do you want to continue (y/n)?\n");
      fgets(buf, sizeof buf, stdin);
      if (buf[0]!='y' && buf[0]!='Y') {
        printf("Abort submission.\n");
        return 1;
      }
    }

    if ( (p = strchr(after, '.')) == NULL) {
      for (p = after; *p; p++) {
        if (!isdigit(*p) ) {
          fprintf(stderr, "Not a valid dependence [%s]\n", after);
          return 1;
        }
      }
      if (cluster[0] != '\0') {
        strcat(strcat(after, "."), cluster);
        printf("complete dependence as [%s].\n", after);
      } else {
        printf("Incomplete dependence [%s].\n", after);
        return 1;
      }
    }
  }

  // submit a sequence of jobs one after one
  for (i = 1; i <= nrep; i++) {
    if (after[0] == '\0') {
      option[0]='\0';
    } else {
      if ( check_after(after) != 0) {
        fprintf(stderr, "i=%d: the depended job [%s] is invalid.\n", i, after);
        return 1;
      }
      sprintf(option, "-W depend=afterok:%s ", after);
    }
    sprintf(cmd, "qsub %s%s > %s", option, fname, SMTFILE);
    docmd(cmd, i, !verbose);

    // read the output
    if ( (fp = fopen(SMTFILE, "r")) == NULL) {
      fprintf(stderr, "Cannot read output file [%s]\n", SMTFILE);
      return 1;
    }
    after[0]='\0';
    fgets(after, sizeof after, fp);
    after[strlen(after)-1]='\0';
    if (verbose) printf("JOB %2d: [%s]\n", i, after);
    fclose(fp);
  }

  remove(SMTFILE); // cleanup
  return 0;
}


#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include "include/rcore.h"
#include <stdbool.h>

#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))


struct ContingencyTable {
   char  X[3];
   char  Y[3];
   char  Z[3];
   int   ***n;
   int   **ni;
   int   **nj;
   int   *nk;
};
const int TABLE_SIZE = 20;
struct ContingencyTable GLOBAL_TABLES[TABLE_SIZE];

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP gsquare, SEXP adjusted) {

int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double *res = NULL;
SEXP result;
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  if (isTRUE(adjusted))
    res[0] = c_chisqtest(xx, llx, yy, lly, num, res + 1, MI_ADF, isTRUE(gsquare));
  else
    res[0] = c_chisqtest(xx, llx, yy, lly, num, res + 1, MI, isTRUE(gsquare));

  UNPROTECT(1);

  return result;

}/*MI*/

/* unconditional parametric asymptotic tests for categorical data. */
double c_chisqtest(int *xx, int llx, int *yy, int lly, int num, double *df,
    test_e test, int scale) {

  clock_t start, end, setup, conting, checks, cleanup, stat;
  start = clock();
  int **n = NULL, *ni = NULL, *nj = NULL, ncomplete = 0, adj = IS_ADF(test);
  double res = 0;
  setup = clock();

  /* initialize the contingency table and the marginal frequencies. */
  ncomplete = fill_2d_table(xx, yy, &n, &ni, &nj, llx, lly, num);
  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? df_adjust(ni, llx, nj, lly) : (llx - 1) * (lly - 1);

  /* if there are no complete data points, return independence. */
  if (ncomplete == 0)
    goto free_and_return;

  /* if there are less than 5 observations per cell on average, assume the
   * test does not have enough power and return independence. */
  if (adj)
    if (ncomplete < 5 * llx * lly)
      goto free_and_return;

  checks = clock();
  /* compute the mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = mi_kernel(n, ni, nj, llx, lly, ncomplete) / ncomplete;
  else if ((test == X2) || (test == X2_ADF))
    res = x2_kernel(n, ni, nj, llx, lly, ncomplete);
  /* rescale to match the G^2 test. */
  if (scale)
    res *= 2 * ncomplete;
  stat = clock();

free_and_return:

  Free2D(n, llx);
  Free1D(ni);
  Free1D(nj);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (checks - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - checks)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "%f,%f,%f,%f,%f\n", time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}/*C_CHISQTEST*/


bool use_3d_table_buffer(const char* x, const char* y, const char* sx, int ****n, int ***ni, int ***nj, int **nk) {
  for (int i = 0; i < TABLE_SIZE; i++) {
    if ((strcmp(x,GLOBAL_TABLES[i].X) == 0) && (strcmp(y,GLOBAL_TABLES[i].Y) == 0) && (strcmp(sx,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    } else if ((strcmp(x,GLOBAL_TABLES[i].X) == 0) && (strcmp(sx,GLOBAL_TABLES[i].Y) == 0) && (strcmp(y,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    } else if ((strcmp(y,GLOBAL_TABLES[i].X) == 0) && (strcmp(x,GLOBAL_TABLES[i].Y) == 0) && (strcmp(sx,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    } else if ((strcmp(y,GLOBAL_TABLES[i].X) == 0) && (strcmp(sx,GLOBAL_TABLES[i].Y) == 0) && (strcmp(x,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    } else if ((strcmp(sx,GLOBAL_TABLES[i].X) == 0) && (strcmp(x,GLOBAL_TABLES[i].Y) == 0) && (strcmp(y,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    } else if ((strcmp(sx,GLOBAL_TABLES[i].X) == 0) && (strcmp(y,GLOBAL_TABLES[i].Y) == 0) && (strcmp(x,GLOBAL_TABLES[i].Z) == 0)) {
      return true;
    }
  }
  return false;
}

void load_3d_table_into_buffer(int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz, int num) {
      return;
    }

/* conditional mutual information, to be used in C code.
  xx - all observations of variable x
  llx - num of categories of x
  yy - all observations of variable y
  lly - num of categories of y
  zz - all observations of variable z
  llz - num of categories of z
  num - number of observations
  df - degrees of freedom
  test - symbol for the test to use
  scale - mostly zero...

*/

double c_cchisqtest_better(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, int scale, const char *x, const char *y, const char* sx, int sepset_length) {
  if (sepset_length == 1) {
    Rprintf("%s %s, %s %d\n", x, y, sx, sepset_length);
  }
  clock_t start, end, setup, conting, checks, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0, adj = IS_ADF(test);
  double res = 0;
  bool buffered = false;

  setup = clock();
  //only if there is one conditional variable
  if (sepset_length == 1) {
    buffered = use_3d_table_buffer(x,y,sx, &n, &ni, &nj, &nk);
  }
   /* initialize the contingency table and the marginal frequencies. */
   if (!buffered) {
     ncomplete = fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
     load_3d_table_into_buffer(&n, &ni, &nj, &nk, llx, lly, llz, num);
   }
  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? cdf_adjust(ni, llx, nj, lly, llz) : (llx - 1) * (lly - 1) * llz;

  /* if there are no complete data points, return independence. */
  if (ncomplete == 0)
    goto free_and_return;

  /* if there are less than 5 observations per cell on average, assume the
   * test does not have enough power and return independence. */
  if (adj)
    if (ncomplete < 5 * llx * lly * llz)
      goto free_and_return;
  checks = clock();
  /* compute the conditional mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = cmi_kernel(n, ni, nj, nk, llx, lly, llz) / ncomplete;
  else if ((test == X2) || (test == X2_ADF))
    res = cx2_kernel(n, ni, nj, nk, llx, lly, llz);
  /* rescale to match the G^2 test. */
  if (scale)
    res *= 2 * ncomplete;
  stat = clock();

free_and_return:

  Free3D(n, llz, llx);
  Free2D(ni, llz);
  Free2D(nj, llz);
  Free1D(nk);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (checks - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - checks)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "%f,%f,%f,%f,%f\n", time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}

double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, int scale) {
  clock_t start, end, setup, conting, checks, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0, adj = IS_ADF(test);
  double res = 0;

  setup = clock();

   /* initialize the contingency table and the marginal frequencies. */
  ncomplete = fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? cdf_adjust(ni, llx, nj, lly, llz) : (llx - 1) * (lly - 1) * llz;

  /* if there are no complete data points, return independence. */
  if (ncomplete == 0)
    goto free_and_return;

  /* if there are less than 5 observations per cell on average, assume the
   * test does not have enough power and return independence. */
  if (adj)
    if (ncomplete < 5 * llx * lly * llz)
      goto free_and_return;
  checks = clock();
  /* compute the conditional mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = cmi_kernel(n, ni, nj, nk, llx, lly, llz) / ncomplete;
  else if ((test == X2) || (test == X2_ADF))
    res = cx2_kernel(n, ni, nj, nk, llx, lly, llz);
  /* rescale to match the G^2 test. */
  if (scale)
    res *= 2 * ncomplete;
  stat = clock();

free_and_return:

  Free3D(n, llz, llx);
  Free2D(ni, llz);
  Free2D(nj, llz);
  Free1D(nk);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (checks - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - checks)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "%f,%f,%f,%f,%f\n", time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}



/* compute the mutual information from the joint and marginal frequencies. */
double mi_kernel(int **n, int *nrowt, int *ncolt, int nrow, int ncol,
    int length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < nrow; i++)
    for (j = 0; j < ncol; j++)
      res += MI_PART(n[i][j], nrowt[i], ncolt[j], length);

  return res;

}/*MI_KERNEL*/

/* compute Pearson's X^2 coefficient from the joint and marginal frequencies. */
double x2_kernel(int **n, int *nrowt, int *ncolt, int nrow, int ncol,
    int length) {

int i = 0, j = 0;
double res = 0, expected = 0;

  for (i = 0; i < nrow; i++)
    for (j = 0; j < ncol; j++) {

      expected = nrowt[i] * (double)ncolt[j] / length;

      if (expected != 0)
        res += (n[i][j] - expected) * (n[i][j] - expected) / expected;

    }/*FOR*/

  return res;

}/*X2_KERNEL*/

/* compute the conditional mutual information from the joint and marginal
 * frequencies. */
double cmi_kernel(int ***n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < nl; k++)
    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++)
        res += MI_PART(n[k][i][j], nrowt[k][i], ncolt[k][j], ncond[k]);

  return res;

}/*CMI_KERNEL*/

/* compute the Pearson's conditional X^2 coefficient from the joint and
 * marginal frequencies. */
double cx2_kernel(int ***n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl) {

int i = 0, j = 0, k = 0;
double expected = 0, res = 0;

  for (k = 0; k < nl; k++) {

    if (ncond[k] == 0) continue;

    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++) {

       expected = nrowt[k][i] * (double)ncolt[k][j] / ncond[k];

       if (expected != 0)
          res += (n[k][i][j] - expected) * (n[k][i][j] - expected) / expected;

      }/*FOR*/

  }/*FOR*/

  return res;

}/*CX2_KERNEL*/

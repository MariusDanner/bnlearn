#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include "include/rcore.h"
#include <stdbool.h>

#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))

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

  if (test != X2) {
    Rprintf("This test can't be used in that way/n");
    return -1.0;
  }
  clock_t start, end, setup, conting, degrees, cleanup, stat;
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
  degrees = clock();
  /* if there are no complete data points, return independence. */
  if (ncomplete == 0)
    goto free_and_return;
  res = x2_kernel(n, ni, nj, llx, lly, ncomplete);

  stat = clock();

free_and_return:

  Free2D(n, llx);
  Free1D(ni);
  Free1D(nj);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (degrees - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - degrees)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "0,0,%f,%f,%f,%f,%f\n", time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}/*C_CHISQTEST*/

double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, int scale, int sepset_length) {
  if (test != X2) {
    Rprintf("This test can't be used in that way/n");
    return -1.0;
  }
  clock_t start, end, setup, conting, degrees, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0;
  double res = 0;

  setup = clock();
  ncomplete = fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = (llx - 1) * (lly - 1) * llz;
  degrees = clock();
  if (ncomplete == 0)
    goto free_and_return;
  res = cx2_kernel(n, ni, nj, nk, llx, lly, llz);
  stat = clock();

free_and_return:

  Free3D(n, llz, llx);
  Free2D(ni, llz);
  Free2D(nj, llz);
  Free1D(nk);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (degrees - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - degrees)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "%d,0,%f,%f,%f,%f,%f\n",sepset_length, time1, time2, time3, time4, time5);
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

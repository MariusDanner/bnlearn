#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include <stdbool.h>

struct ContingencyTable {
   char  X[3];
   char  Y[3];
   char  Z[3];
   int   ***n;
   int   **ni;
   int   **nj;
   int   *nk;
};

#define TABLE_SIZE 20

//const int TABLE_SIZE = 20;
struct ContingencyTable GLOBAL_TABLES[TABLE_SIZE];

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
#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include <stdbool.h>
#include "include/hashmap.h"

#define KEY_MAX_LENGTH (256)
#define KEY_PREFIX ("")
#define KEY_COUNT (1024*1024)

map_t conting_hashmap;


struct ContingencyTable {
   char*  X;
   char*  Y;
   char*  Z;
   int   ***n;
   int   **ni;
   int   **nj;
   int   *nk;
};

bool use_3d_table_buffer(const char* x, const char* y, const char* z, int ****n, int ***ni, int ***nj, int **nk) {
  if(conting_hashmap == NULL){
    conting_hashmap = hashmap_new();
  }
  return false;
}

void load_3d_table_into_buffer(int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz, int num) {
      if(conting_hashmap == NULL){
        conting_hashmap = hashmap_new();
      }
      return;
    }

/*  xx pointer to x observations
    llx number of categories of x
    yy pointer to y observations
    lly number of categories of y
    zz pointer to z observations
    llz number of categories of z
    num number of observations?
    df should be degrees of freedom but somehow irrelevant
    test test to use
    scale (test == MI) || (test == MI_ADF)

  */

double c_cchisqtest_better(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, int scale, const char *x, const char *y, const char *z, int sepset_length) {
  clock_t start, end, setup, conting, checks, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0, adj = IS_ADF(test);
  double res = 0;
  bool buffered = false;

  if(sepset_length == 1){
    Rprintf("x: %s, y: %s, z: %s\n", x, y, z);
  }

  setup = clock();
  //only if there is one conditional variable
  if (sepset_length == 1) {
     buffered = use_3d_table_buffer(x, y, z, &n, &ni, &nj, &nk);
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
  fprintf(fp, "%d, %f,%f,%f,%f,%f\n", sepset_length, time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}

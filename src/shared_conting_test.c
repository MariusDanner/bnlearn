#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include <stdbool.h>
#include "include/hashmap.h"

int test_count = 0;

#define KEY_MAX_LENGTH (256)
#define KEY_PREFIX ("")
#define KEY_COUNT (1024*1024)

map_t conting_hashmap;

typedef struct ContingencyTable {
   char X[50];
   char Y[50];
   char Z[50];
   int ***n;
   int **ni;
   int **nj;
   int *nk;
} conting_table_t;

char* order_concat_keys(const char* x, const char* y, const char* z) {
  char* concat = malloc(150 * sizeof(char));
  memset(concat,0,150);
  char x1[50];
  char y1[50];
  char z1[50];
  char temp[50];
  strcpy(x1,x);
  strcpy(y1,y);
  strcpy(z1,z);
  if(strcmp(x1,z1)>0){
    strcpy(temp,x1);
    strcpy(x1,z1);
    strcpy(z1,temp);
  }
  if(strcmp(x1,y1)>0){
    strcpy(temp,x1);
    strcpy(x1,y1);
    strcpy(y1,temp);
  }
  if(strcmp(y1,z1)>0){
    strcpy(temp,y1);
    strcpy(y1,z1);
    strcpy(z1,temp);
  }
  strcpy(concat, x1);
  strcat(concat, y1);
  strcat(concat, z1);
  return concat;
}

int get_permutation_id(const char* x_v, const char* y_v, const char* z_v, const char* x_i, const char* y_i, const char* z_i){
  if(strcmp(x_v, x_i)==0){
    if(strcmp(y_v, y_i)==0){
      return 0;
    } else{
      return 1;
    }
  } else if(strcmp(x_v, y_i)==0){
    if(strcmp(y_v, x_i)==0){
      return 2;
    } else{
      return 3;
    }
  } else{
    if(strcmp(y_v, x_i)==0){
      return 4;
    } else{
      return 5;
    }
  }
}

bool use_3d_table_buffer(const char* x, const char* y, const char* z, int ****n, int ***ni, int ***nj, int **nk, int *llx, int *lly, int *llz, int *perm_id) {
  struct ContingencyTable* value;
  if(conting_hashmap == NULL){
    conting_hashmap = hashmap_new();
  }

  int error = hashmap_get(conting_hashmap,order_concat_keys(x,y,z),(void**)(&value));
  if (error == MAP_OK) {

    *perm_id = get_permutation_id(value->X, value->Y, value->Z, x, y, z);
    *n = value->n;

    return true;
  }
  return false;
}


void load_3d_table_into_buffer(const char* x, const char* y, const char* z,int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz) {
      if(conting_hashmap == NULL){
        conting_hashmap = hashmap_new();
      }
      conting_table_t* value = malloc(sizeof(conting_table_t));
      strcpy(value->X, x);
      strcpy(value->Y, y);
      strcpy(value->Z, z);
      value->n = *n;
      value->ni = *ni;
      value->nj = *nj;
      value->nk = *nk;

      hashmap_put(conting_hashmap,order_concat_keys(x,y,z),value);
    }

/* compute the Pearson's conditional X^2 coefficient from the joint and
 * marginal frequencies. */
double cx2_kernel_better(int ***n, int **ni, int **nj, int *nk, int llx,
    int lly, int llz, int perm_id) {

int i = 0, j = 0, k = 0;
double expected = 0, res = 0;
ni = (int **) Calloc2D(llz, llx, sizeof(int));
nj = (int **) Calloc2D(llz, lly, sizeof(int));
nk = (int *) Calloc1D(llz, sizeof(int));
  switch (perm_id) {
    case 1: {
      for (int i = 0; i < llx; i++)
        for (int j = 0; j < lly; j++)
          for (int k = 0; k < llz; k++) {
            ni[k][i] += n[j][i][k];
            nj[k][j] += n[j][i][k];
            nk[k] += n[j][i][k];
          }
      for (k = 0; k < llz; k++) {
        if (nk[k] == 0) continue;
        for (j = 0; j < lly; j++) {
          for (i = 0; i < llx; i++) {
           expected = ni[k][i] * (double)nj[k][j] / nk[k];
           if (expected != 0)
              res += (n[j][i][k] - expected) * (n[j][i][k] - expected) / expected;
          }/*FOR*/
        }
      }
      break;
    }
    case 2: {
      for (int i = 0; i < llx; i++)
        for (int j = 0; j < lly; j++)
          for (int k = 0; k < llz; k++) {

            ni[k][i] += n[k][j][i];
            nj[k][j] += n[k][j][i];
            nk[k] += n[k][j][i];
          }
      for (k = 0; k < llz; k++) {
        if (nk[k] == 0) continue;
        for (j = 0; j < lly; j++) {
          for (i = 0; i < llx; i++) {
           expected = ni[k][i] * (double)nj[k][j] / nk[k];
           if (expected != 0)
              res += (n[k][j][i] - expected) * (n[k][j][i] - expected) / expected;
          }/*FOR*/
        }
      }
      break;
    }
    case 3: {
      for (int i = 0; i < llx; i++)
        for (int j = 0; j < lly; j++)
          for (int k = 0; k < llz; k++) {

            ni[k][i] += n[i][j][k];
            nj[k][j] += n[i][j][k];
            nk[k] += n[i][j][k];
          }
      for (k = 0; k < llz; k++) {
        if (nk[k] == 0) continue;
        for (j = 0; j < lly; j++) {
          for (i = 0; i < llx; i++) {
           expected = ni[k][i] * (double)nj[k][j] / nk[k];
           if (expected != 0)
              res += (n[i][j][k] - expected) * (n[i][j][k] - expected) / expected;
          }/*FOR*/
        }
      }
      break;
    }
    case 4: {
      for (int i = 0; i < llx; i++)
        for (int j = 0; j < lly; j++)
          for (int k = 0; k < llz; k++) {

            ni[k][i] += n[j][k][i];
            nj[k][j] += n[j][k][i];
            nk[k] += n[j][k][i];
          }
      for (k = 0; k < llz; k++) {
        if (nk[k] == 0) continue;
        for (j = 0; j < lly; j++) {
          for (i = 0; i < llx; i++) {
           expected = ni[k][i] * (double)nj[k][j] / nk[k];
           if (expected != 0)
              res += (n[j][k][i] - expected) * (n[j][k][i] - expected) / expected;
          }/*FOR*/
        }
      }
      break;
    }
    case 5: {
      for (int i = 0; i < llx; i++)
        for (int j = 0; j < lly; j++)
          for (int k = 0; k < llz; k++) {

            ni[k][i] += n[i][k][j];
            nj[k][j] += n[i][k][j];
            nk[k] += n[i][k][j];
          }
      for (k = 0; k < llz; k++) {
        if (nk[k] == 0) continue;
        for (j = 0; j < lly; j++) {
          for (i = 0; i < llx; i++) {
           expected = ni[k][i] * (double)nj[k][j] / nk[k];
           if (expected != 0)
              res += (n[i][k][j] - expected) * (n[i][k][j] - expected) / expected;
          }/*FOR*/
        }
      }
      break;
    }
  }
  int ncomplete = 0;
  for (int k = 0; k < llz; k++)
    ncomplete += nk[k];

  if (ncomplete == 0){
    return -1.0;
  }

  return res;

}/*CX2_KERNEL*/

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
  test_count++;
  Rprintf("%d\n", test_count);
  if (test != X2) {
    Rprintf("This test can't be used in that way/n");
    return -1.0;
  }
  clock_t start, end, setup, conting, checks, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0;
  double res = 0;
  bool buffered = false;

  setup = clock();
  //only if there is one conditional variable
  int perm_id = 0;
  if (sepset_length == 1) {
     buffered = use_3d_table_buffer(x, y, z, &n, &ni, &nj, &nk, &llx, &lly, &llz, &perm_id);
  }

  /* initialize the contingency table and the marginal frequencies. */
  if (!buffered) {
      ncomplete = fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
      if (sepset_length == 1) {
        load_3d_table_into_buffer(x, y, z, &n, &ni, &nj, &nk, llx, lly, llz);
      }
  }

  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = (llx - 1) * (lly - 1) * llz;



  checks = clock();
  /* compute the conditional mutual information or Pearson's X^2. */
  if (!buffered){
    /* if there are no complete data points, return independence. */
    if (ncomplete == 0)
      goto free_and_return;
    res = cx2_kernel(n, ni, nj, nk, llx, lly, llz);
  } else {
    res= cx2_kernel_better(n, ni, nj, nk, llx,lly,llz, perm_id);
  }


  stat = clock();

free_and_return:


  //Free3D(n, llz, llx);
  //Free2D(ni, llz);
  //Free2D(nj, llz);
  //Free1D(nk);
  cleanup = clock();
  double time1 = ((double) (setup - start)) / CLOCKS_PER_SEC;
  double time2 = ((double) (conting - setup)) / CLOCKS_PER_SEC;
  double time3 = ((double) (checks - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - checks)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  FILE *fp = fopen("ci_benchmark.csv", "a");
  fprintf(fp, "%d, %d, %f,%f,%f,%f,%f\n", sepset_length, buffered, time1, time2, time3, time4, time5);
  fclose(fp);
  return res;

}

#include "include/rcore.h"
#include "include/tests.h"
#include "time.h"
#include <stdbool.h>

int test_count = 0;

#define KEY_MAX_LENGTH (256)
#define KEY_PREFIX ("")
#define KEY_COUNT (1024*1024)

SEXP setup_lookup(SEXP n, SEXP nodes) {
  reverse_lookup_hashmap = hashmap_new();
  int* nptr = INTEGER(n);
  var_count = *nptr;
  for (int i = 0; i < *nptr; i++) {
    char* string_x = CHAR(STRING_ELT(nodes,i));
    char* string_y = malloc(strlen(string_x));
    strcpy(string_y, string_x);
    hash_value_t* value = malloc(sizeof(hash_value_t));
    value->id = i;
    hashmap_put(reverse_lookup_hashmap, string_y, value);
  }
  table_buffer = calloc((*nptr) * (*nptr) * (*nptr), sizeof(conting_table_t*));
  return n;
}

SEXP cleanup_lookup(){
  hashmap_free(reverse_lookup_hashmap);
  for (int i = 0; i < var_count * var_count * var_count; i++) {
    if (table_buffer != NULL) {
      free(table_buffer[i]);
    }
  }
  free(table_buffer);
}

int get_key(int x, int y, int z) {
  int swap;
  if (x > z) {
    swap = z;
    z = x;
    x = swap;
  }
  if (x > y) {
    swap = y;
    y = x;
    x = swap;
  }
  if (y > z){
    swap = y;
    y = z;
    z = swap;
  }
  return (var_count * var_count * x + var_count * y + z);

}

int get_permutation_id(int x_v, int y_v, int z_v, int x_i, int y_i, int z_i){
  if(x_v == x_i){
    if(y_v == y_i){
      return 0;
    } else{
      return 1;
    }
  } else if(x_v == y_i){
    if(y_v == x_i){
      return 2;
    } else{
      return 3;
    }
  } else{
    if(y_v == x_i){
      return 4;
    } else{
      return 5;
    }
  }
}

bool use_3d_table_buffer(int x, int y, int z, int ****n, int ***ni, int ***nj, int **nk, int *llx, int *lly, int *llz, int *perm_id) {
  struct ContingencyTable* value;

  value = table_buffer[get_key(x,y,z)];
  if (value != NULL) {
    //Rprintf("buffered\n");
    *perm_id = get_permutation_id(value->x, value->y, value->z, x, y, z);
    *n = value->n;

    return true;
  }
  return false;
}


void load_3d_table_into_buffer(int x, int y, int z,int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz) {
      conting_table_t* value = malloc(sizeof(conting_table_t));
      value->x = x;
      value->y = y;
      value->z = z;
      value->n = *n;
      value->ni = *ni;
      value->nj = *nj;
      value->nk = *nk;

      table_buffer[get_key(x,y,z)] = value;
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
    int num, double *df, test_e test, int scale, char *x, char *y, char *z, int sepset_length) {
  hash_value_t* xid;
  hash_value_t *yid;
  hash_value_t* zid;
  if (sepset_length == 1) {
    hashmap_get(reverse_lookup_hashmap, x, (void**)(&xid));
    hashmap_get(reverse_lookup_hashmap, y, (void**)(&yid));
    hashmap_get(reverse_lookup_hashmap, z, (void**)(&zid));
  }

  if (test != X2) {
    Rprintf("This test can't be used in that way/n");
    return -1.0;
  }
  clock_t start, end, setup, conting, degrees, cleanup, stat;
  start = clock();
  int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
  int ncomplete = 0;
  double res = 0;
  bool buffered = false;
  setup = clock();
  //only if there is one conditional variable
  int perm_id = 0;
  if (sepset_length == 1) {
    buffered = use_3d_table_buffer(xid->id, yid->id, zid->id, &n, &ni, &nj, &nk, &llx, &lly, &llz, &perm_id);
  }
  /* initialize the contingency table and the marginal frequencies. */
  if (!buffered) {
    // Rprintf("%s %s %s\n", x, y, z);
      ncomplete = fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
      if (sepset_length == 1) {
        load_3d_table_into_buffer(xid->id, yid->id, zid->id, &n, &ni, &nj, &nk, llx, lly, llz);
      }
  }
  conting = clock();
  /* compute the degrees of freedom. */
  if (df)
    *df = (llx - 1) * (lly - 1) * llz;



  degrees = clock();
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
  double time3 = ((double) (degrees - conting)) / CLOCKS_PER_SEC;
  double time4 = ((double) (stat - degrees)) / CLOCKS_PER_SEC;
  double time5 = ((double) (cleanup - stat)) / CLOCKS_PER_SEC;
  // FILE *fp = fopen("ci_benchmark.csv", "a");
  // fprintf(fp, "%d, %d, %f,%f,%f,%f,%f\n", sepset_length, buffered, time1, time2, time3, time4, time5);
  // fclose(fp);
  return res;

}

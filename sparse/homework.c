#include "fem.h"
#include <time.h>

typedef struct sparseMatrix {
  int size;
  int nnz;
  int *col;
  int *rptr;
  double *val;
} sparseMatrix;

void sparseMatrixFree(sparseMatrix *sp) {
  free(sp->col);
  free(sp->rptr);
  free(sp->val);
  free(sp);
}

sparseMatrix* to_sparse(double **A, int size) {
  int nnz = 0;
  for (int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      nnz += (A[i][j] != 0);
    }
  }
  int* col = malloc(nnz * sizeof(int));
  int* rptr = malloc((size + 1) * sizeof(int));
  double* val = malloc(nnz * sizeof(double));
  
  nnz = 0;
  for(int i = 0; i < size; i++) {
    rptr[i] = nnz;
    for(int j = 0; j < size; j++) {
      if(A[i][j] != 0) {
        col[nnz] = j;
        val[nnz] = A[i][j];
        nnz++;
      }
    }
  }
  rptr[size] = nnz;

  sparseMatrix* sp = malloc(sizeof(sparseMatrix));
  sp->size = size;
  sp->nnz = nnz;
  sp->col = col;
  sp->rptr = rptr;
  sp->val = val;
  return sp;
}

static inline void spmv(const sparseMatrix* sp, const double* x, double* y) {
  for(int i = 0; i < sp->size; i++) {
    double s = 0;
    for(int j = sp->rptr[i]; j < sp->rptr[i+1]; j++) {
      s += sp->val[j] * x[sp->col[j]];
    }
    y[i] = s;
  }
}

static inline void residual(const sparseMatrix* sp, const double* x, const double* b, double* r) {
  spmv(sp, x, r);
  for(int i = 0; i < sp->size; i++) {
    r[i] = b[i] - r[i];
  }
}

static inline double dot(const double* x, const double* y, int size) {
  double s = 0;
  for(int i = 0; i < size; i++) {
    s += x[i] * y[i];
  }
  return s;
}

static inline void axpy(double* x, const double* y, double a, int size) {
  for(int i = 0; i < size; i++) {
    x[i] += a * y[i];
  }
}

double *solve_cg(femFullSystem *mySystem) {
  int size = mySystem->size;
  double **A = mySystem->A;
  double *B = mySystem->B;
  clock_t start, stop;

  start = clock();
  sparseMatrix* sp = to_sparse(A, size);
  stop = clock();
  printf("Time to create sparse mat: %f ms\n", 1000 * (double)(stop - start) / CLOCKS_PER_SEC);
  start = clock();

  int niter = 0;
  // TODO :-)

  stop = clock();
  sparseMatrixFree(sp);
  printf("Time to solve: %f ms, %d iters (%d nodes)\n", 1000 * (double)(stop - start) / CLOCKS_PER_SEC, niter, size);
  return (mySystem->B);
}


#ifndef SPARSEUTILS
#define SPARSEUTILS

template <typename SCALAR>
double spdot(const double * __restrict v, const SCALAR * __restrict a,
	     const int * __restrict jcol, int nnz)
{
  int k;
  double dot = 0.0;

  for (k = 0;  k < nnz;  k++) {
    dot += v[jcol[k]] * a[k];
  }
  return dot;
}


template <typename SCALAR>
void spaxpy(double alpha, const SCALAR * __restrict x,
	    const int * __restrict jcol,
	    int nnz, double * __restrict y)
{
  int k;
  if (alpha == 1.0) {
    for (k = 0;  k < nnz;  k++) {
      y[jcol[k]] += x[k];
    }
  } else if (alpha == -1.0) {
    for (k = 0;  k < nnz;  k++) {
      y[jcol[k]] -= alpha * x[k];
    }
  } else if (alpha == 0.0) {
    return;
  } else {
    for (k = 0;  k < nnz;  k++) {
      y[jcol[k]] += alpha * x[k];
    }
  }
}

template <typename SCALAR>
void spsyr(int n, double alpha, const SCALAR * __restrict x,
	   const int * __restrict jcol, int nnz,
	   double * __restrict sparse_matrix)
{
  int ii;
  for (ii = 0;  ii < nnz;  ii++) {
    int j, jj;
    int i = jcol[ii];
    double * row = sparse_matrix + i * n;
    double temp = alpha * x[ii];
    for (jj = 0;  jj < nnz && (j = jcol[jj]) <= i;  jj++) {
      row[j] += x[jj] * temp;
    }
  }
}


template <typename SCALAR>
void spspr_1(int n, double alpha, const SCALAR * __restrict x, 
	     const int * __restrict jcol, int nnz,
	     double * __restrict packed_matrix)
{
  int ii;
  for (ii = 0;  ii < nnz;  ii++) {
    int j, jj;
    int i = jcol[ii];
    double * row = packed_matrix + i * (i + 1) / 2;
    double temp = alpha * x[ii];
    for (jj = 0;  jj < nnz && (j = jcol[jj]) < i;  jj++) {
      row[j] += x[jj] * temp;
    }
  }
}

#endif

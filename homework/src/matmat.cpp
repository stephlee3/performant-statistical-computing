#include <Rcpp.h>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// [[Rcpp::export]]
NumericMatrix row_col_dot_matmat(
    const NumericMatrix& A, const NumericMatrix& B
  ) {
  int n_row_A = A.nrow();
  int interior_dim = A.ncol();
  int n_col_B = B.ncol();
  if (interior_dim != B.nrow()) {
    Rcpp::stop("Incompatible dimensions");
  }
  NumericMatrix C(n_row_A, n_col_B);
  // Fill in; remember to use (,) instead of [,] for accessing matrix elements
  for (int i = 0; i < n_row_A; ++i){
    for (int j = 0; j < n_col_B; ++j){
      for(int k = 0; k < interior_dim; ++k){
        C(i,j) += A(i,k) * B(k,j);
      } 
    }
  }
  return C;
}

// [[Rcpp::export]]
NumericMatrix col_oriented_matmat(
    const NumericMatrix& A, const NumericMatrix& B
  ) {
  int n_row_A = A.nrow();
  int interior_dim = A.ncol();
  int n_col_B = B.ncol();
  if (interior_dim != B.nrow()) {
    Rcpp::stop("Incompatible dimensions");
  }
  NumericMatrix C(n_row_A, n_col_B);
  
  for (int j = 0; j < n_col_B; ++j) {
    // Calculate C[, j] = A %*% B[_, j] in column-oriented manner for each j
    // Fill in
    NumericVector tmp(n_row_A);
    for(int k = 0; k < interior_dim; ++k){
      tmp = tmp + A(Rcpp::_ , k) * B(k, j);
    }
    C(Rcpp::_, j) = tmp;
  }
  return C;
}

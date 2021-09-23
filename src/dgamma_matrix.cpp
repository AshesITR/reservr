#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dgamma_matrix(NumericVector x, NumericVector shape, double scale) {
  int n = x.size(), m = shape.size();
  NumericMatrix res(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      res(i, j) = R::dgamma(x[i], shape[j], scale, 0);
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix pgamma_diff_matrix(NumericVector lower, NumericVector upper, NumericVector shape, double scale) {
  int n = lower.size(), m = shape.size();
  NumericMatrix res(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      res(i, j) = R::pgamma(upper[i], shape[j], scale, 1, 0) - R::pgamma(lower[i], shape[j], scale, 1, 0);
    }
  }
  return res;
}

/*** R
dgamma_matrix(1:10, 1:2, 1)
*/

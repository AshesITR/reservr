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
NumericMatrix pgamma_diff_matrix(const NumericVector lower, const NumericVector upper, const NumericVector shape, const NumericVector scale) {
  R_xlen_t n = std::max(lower.size(), scale.size()), m = shape.size();
  NumericMatrix res(n, m);
  R_xlen_t i_scale = 0, d_scale = scale.size() > 1 ? 1 : 0;
  for (R_xlen_t i = 0; i < n; i++) {
    for (R_xlen_t j = 0; j < m; j++) {
      res(i, j) = R::pgamma(upper[i], shape[j], scale[i_scale], 1, 0) - R::pgamma(lower[i], shape[j], scale[i_scale], 1, 0);
    }
    i_scale += d_scale;
  }
  return res;
}

/*** R
dgamma_matrix(1:10, 1:2, 1)
*/

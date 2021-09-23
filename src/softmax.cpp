#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat softmax_mat(arma::mat x) {
  x.each_col() -= max(x, 1);
  mat ex = exp(x);
  vec row_sums = sum(ex, 1);
  ex.each_col() /= row_sums;
  return ex;
}

// [[Rcpp::export]]
std::vector<double> softmax_vec(arma::vec x) {
  x -= max(x);
  vec ex = exp(x);
  return Rcpp::as<std::vector<double>>(Rcpp::wrap(ex / sum(ex)));
}

// [[Rcpp::export]]
arma::mat dsoftmax_vec(arma::vec x) {
  x -= max(x);
  vec ex = exp(x);
  ex /= sum(ex);
  return diagmat(ex) - ex * ex.t();
}

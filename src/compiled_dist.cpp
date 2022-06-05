// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// efficiently aggregate mixture densities / probabilities

template <typename T>
arma::vec aggregate_mixture(const arma::mat compmat, const T probs);

template <>
arma::vec aggregate_mixture<arma::subview_cols<double>>(const arma::mat compmat, const arma::subview_cols<double> probs) {
  arma::vec res = arma::sum(compmat % probs, 1);
  res /= arma::sum(probs, 1);
  return res;
}

template <>
arma::vec aggregate_mixture<arma::vec>(const arma::mat compmat, const arma::vec probs) {
  arma::vec res = compmat * probs;
  res /= arma::accu(probs);
  return res;
}

// efficiently compute mixture densities with possibly fixed probs

template <typename T>
arma::vec dist_mixture_density_impl(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const arma::uvec is_discrete, const T probs) {
  int k = comp_densities.size();
  int n = std::max(x.n_elem, params.n_rows);
  int i = 0;
  arma::mat compdens(n, k);
  SEXP curr_params;
  for (int j = 0; j < k; j++) {
    if (param_sizes[j] > 0) {
      curr_params = Rcpp::wrap(params.cols(i, i + param_sizes[j] - 1));
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }
    Rcpp::Function curr_fun = comp_densities[j];
    compdens.col(j) = Rcpp::as<arma::vec>(curr_fun(x, curr_params, false));
  }
  if (!all(is_discrete) && any(is_discrete)) {
    // mixed type => need to conditionally zero continuous densities
    arma::uvec cont_components = find(is_discrete == 0);
    arma::uvec disc_points  = any(compdens.cols(find(is_discrete)), 1);
    compdens.submat(disc_points, cont_components).fill(0.0);
  }
  arma::vec res = aggregate_mixture(compdens, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_mixture_density_free(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const arma::uvec is_discrete) {
  return dist_mixture_density_impl(x, params, log_p, param_sizes, comp_densities, is_discrete, params.tail_cols(comp_densities.size()));
}

// [[Rcpp::export]]
arma::vec dist_mixture_density_fixed(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const arma::uvec is_discrete, const arma::vec probs) {
  return dist_mixture_density_impl(x, params, log_p, param_sizes, comp_densities, is_discrete, probs);
}

// efficiently compute mixture probabilities with possibly fixed probs

template <typename T>
arma::vec dist_mixture_probability_impl(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_probabilities, const T probs) {
  int k = comp_probabilities.size();
  int n = std::max(q.n_elem, params.n_rows);
  int i = 0;
  arma::mat compprob(n, k);
  SEXP curr_params;
  for (int j = 0; j < k; j++) {
    if (param_sizes[j] > 0) {
      curr_params = Rcpp::wrap(params.cols(i, i + param_sizes[j] - 1));
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }
    Rcpp::Function curr_fun = comp_probabilities[j];
    compprob.col(j) = Rcpp::as<arma::vec>(curr_fun(q, curr_params, lower_tail, false));
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_mixture_probability_free(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_probabilities) {
  return dist_mixture_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_probabilities, params.tail_cols(comp_probabilities.size()));
}

// [[Rcpp::export]]
arma::vec dist_mixture_probability_fixed(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_probabilities, const arma::vec probs) {
  return dist_mixture_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_probabilities, probs);
}

// efficiently compute mixture interval probabilities with possibly fixed probs

template <typename T>
arma::vec dist_mixture_iprobability_impl(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const T probs) {
  int k = comp_iprobabilities.size();
  int n = std::max(qmin.n_elem, params.n_rows);
  int i = 0;
  arma::mat compprob(n, k);
  SEXP curr_params;
  for (int j = 0; j < k; j++) {
    if (param_sizes[j] > 0) {
      curr_params = Rcpp::wrap(params.cols(i, i + param_sizes[j] - 1));
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }
    Rcpp::Function curr_fun = comp_iprobabilities[j];
    compprob.col(j) = Rcpp::as<arma::vec>(curr_fun(qmin, qmax, curr_params, false));
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_mixture_iprobability_free(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities) {
  return dist_mixture_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, params.tail_cols(comp_iprobabilities.size()));
}

// [[Rcpp::export]]
arma::vec dist_mixture_iprobability_fixed(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs) {
  return dist_mixture_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, probs);
}

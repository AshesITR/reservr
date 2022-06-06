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
arma::vec aggregate_mixture<arma::mat>(const arma::mat compmat, const arma::mat probs) {
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

template <typename T>
int num_components(const T probs);

template <>
int num_components<arma::subview_cols<double>>(const arma::subview_cols<double> probs) {
  return probs.n_cols;
}

template <>
int num_components<arma::mat>(const arma::mat probs) {
  return probs.n_cols;
}

template <>
int num_components<arma::vec>(const arma::vec probs) {
  return probs.n_elem;
}

template <typename T>
bool is_matrix(const T x);

template <>
bool is_matrix<arma::vec>(const arma::vec x) {
  return false;
}

template <>
bool is_matrix<arma::mat>(const arma::mat x) {
  return true;
}

template <>
bool is_matrix<arma::subview_cols<double>>(const arma::subview_cols<double> x) {
  return true;
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

// efficiently compute erlang mixture densities with possibly fixed parameters

template <typename TP, typename TS>
arma::vec dist_erlangmix_density_impl(const arma::vec x, bool log_p, const TP probs, const arma::vec scale, const TS shapes) {
  int k = num_components(probs);
  int n = std::max(std::max(x.n_elem, probs.n_rows), std::max(scale.n_elem, shapes.n_rows));
  bool shape_is_matrix = is_matrix(shapes);
  arma::mat compdens(n, k);
  int i_x = 0, d_x = x.n_elem > 1 ? 1 : 0;
  int i_s = 0, d_s = scale.n_elem > 1 ? 1 : 0;
  double curr_shape;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      if (shape_is_matrix) {
        curr_shape = shapes(i, j);
      } else {
        curr_shape = shapes[j];
      }
      compdens(i, j) = R::dgamma(x[i_x], curr_shape, scale[i_s], 0);
    }
    i_x += d_x;
    i_s += d_s;
  }
  arma::vec res = aggregate_mixture(compdens, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_free(const arma::vec x, const arma::mat params, bool log_p) {
  int k = (params.n_cols - 1) / 2;
  return dist_erlangmix_density_impl(x, log_p, params.tail_cols(k), params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_shape(const arma::vec x, const arma::mat params, bool log_p, const arma::vec shapes) {
  int k = shapes.n_elem;
  return dist_erlangmix_density_impl(x, log_p, params.tail_cols(k), params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_scale(const arma::vec x, const arma::mat params, bool log_p, const arma::vec scale) {
  int k = params.n_cols / 2;
  return dist_erlangmix_density_impl(x, log_p, params.tail_cols(k), scale, params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_probs(const arma::vec x, const arma::mat params, bool log_p, const arma::vec probs) {
  int k = probs.n_elem;
  return dist_erlangmix_density_impl(x, log_p, probs, params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_probs_scale(const arma::vec x, const arma::mat params, bool log_p, const arma::vec probs, const arma::vec scale) {
  return dist_erlangmix_density_impl(x, log_p, probs, scale, params);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_probs_shape(const arma::vec x, const arma::mat params, bool log_p, const arma::vec probs, const arma::vec shapes) {
  return dist_erlangmix_density_impl(x, log_p, probs, params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_scale_shape(const arma::vec x, const arma::mat params, bool log_p, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_density_impl(x, log_p, params, scale, shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_density_fixed_probs_scale_shape(const arma::vec x, bool log_p, const arma::vec probs, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_density_impl(x, log_p, probs, scale, shapes);
}

// efficiently compute erlang mixture probabilities with possibly fixed parameters

template <typename TP, typename TS>
arma::vec dist_erlangmix_probability_impl(const arma::vec q, bool lower_tail, bool log_p, const TP probs, const arma::vec scale, const TS shapes) {
  int k = num_components(probs);
  int n = std::max(std::max(q.n_elem, probs.n_rows), std::max(scale.n_elem, shapes.n_rows));
  bool shape_is_matrix = is_matrix(shapes);
  arma::mat compprob(n, k);
  int i_q = 0, d_q = q.n_elem > 1 ? 1 : 0;
  int i_s = 0, d_s = scale.n_elem > 1 ? 1 : 0;
  double curr_shape;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      if (shape_is_matrix) {
        curr_shape = shapes(i, j);
      } else {
        curr_shape = shapes[j];
      }
      compprob(i, j) = R::pgamma(q[i_q], curr_shape, scale[i_s], lower_tail, 0);
    }
    i_q += d_q;
    i_s += d_s;
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_free(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p) {
  int k = (params.n_cols - 1) / 2;
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, params.tail_cols(k), params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_shape(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec shapes) {
  int k = shapes.n_elem;
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, params.tail_cols(k), params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_scale(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec scale) {
  int k = params.n_cols / 2;
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, params.tail_cols(k), scale, params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_probs(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec probs) {
  int k = probs.n_elem;
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, probs, params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_probs_scale(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec probs, const arma::vec scale) {
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, probs, scale, params);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_probs_shape(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec probs, const arma::vec shapes) {
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, probs, params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_scale_shape(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, params, scale, shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_probability_fixed_probs_scale_shape(const arma::vec q, bool lower_tail, bool log_p, const arma::vec probs, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_probability_impl(q, lower_tail, log_p, probs, scale, shapes);
}

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
const unsigned int num_components(const T x);

template <>
const unsigned int num_components<arma::subview_cols<double>>(const arma::subview_cols<double> x) {
  return x.n_cols;
}

template <>
const unsigned int num_components<arma::mat>(const arma::mat x) {
  return x.n_cols;
}

template <>
const unsigned int num_components<arma::vec>(const arma::vec x) {
  return x.n_elem;
}

template <typename T>
const unsigned int num_observations(const T x);

template <>
const unsigned int num_observations<arma::subview_cols<double>>(const arma::subview_cols<double> x) {
  return x.n_rows;
}

template <>
const unsigned int num_observations<arma::mat>(const arma::mat x) {
  return x.n_rows;
}

template <>
const unsigned int num_observations<arma::vec>(const arma::vec x) {
  return 1;
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

template <typename T>
arma::vec column_or_element(const T x, int col);

template <>
arma::vec column_or_element<arma::vec>(const arma::vec x, int col) {
  return x.subvec(col, col);
}

template <>
arma::vec column_or_element<arma::mat>(const arma::mat x, int col) {
  return x.col(col);
}

template <>
arma::vec column_or_element<arma::subview_cols<double>>(const arma::subview_cols<double> x, int col) {
  return x.col(col);
}

arma::uvec find_relevant(const arma::vec xmin, const arma::vec xmax, const arma::vec u_lo, const arma::vec e_lo, const arma::vec u_hi, const arma::vec e_hi) {
  if (u_lo.n_elem > 1 && e_lo.n_elem > 1 && u_hi.n_elem > 1 && e_hi.n_elem > 1) { // all vec
    return arma::find(u_lo - e_lo <= xmax && xmin <= u_hi + e_hi);

  } else if (u_lo.n_elem > 1 && e_lo.n_elem > 1 && u_hi.n_elem > 1) { // 1 const
    return arma::find(u_lo - e_lo <= xmax && xmin <= u_hi + e_hi[0]);
  } else if (u_lo.n_elem > 1 && e_lo.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo - e_lo <= xmax && xmin <= u_hi[0] + e_hi);
  } else if (u_lo.n_elem > 1 && u_hi.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo - e_lo[0] <= xmax && xmin <= u_hi + e_hi);
  } else if (e_lo.n_elem > 1 && u_hi.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo <= xmax && xmin <= u_hi + e_hi);

  } else if (u_lo.n_elem > 1 && e_lo.n_elem > 1) { // 2 const
    return arma::find(u_lo - e_lo <= xmax && xmin <= u_hi[0] + e_hi[0]);
  } else if (u_lo.n_elem > 1 && u_hi.n_elem > 1) {
    return arma::find(u_lo - e_lo[0] <= xmax && xmin <= u_hi + e_hi[0]);
  } else if (e_lo.n_elem > 1 && u_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo <= xmax && xmin <= u_hi + e_hi[0]);
  } else if (u_lo.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo - e_lo[0] <= xmax && xmin <= u_hi[0] + e_hi);
  } else if (e_lo.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo <= xmax && xmin <= u_hi[0] + e_hi);
  } else if (u_hi.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo[0] <= xmax && xmin <= u_hi + e_hi);

  } else if (u_lo.n_elem > 1) { // 3 const
    return arma::find(u_lo - e_lo[0] <= xmax && xmin <= u_hi[0] + e_hi[0]);
  } else if (e_lo.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo <= xmax && xmin <= u_hi[0] + e_hi[0]);
  } else if (u_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo[0] <= xmax && xmin <= u_hi + e_hi[0]);
  } else if (e_hi.n_elem > 1) {
    return arma::find(u_lo[0] - e_lo[0] <= xmax && xmin <= u_hi[0] + e_hi);

  } else { // all const
    return arma::find(u_lo[0] - e_lo[0] <= xmax && xmin <= u_hi[0] + e_hi[0]);
  }
}

arma::uvec find_high(const arma::vec x, const arma::vec u_hi, const arma::vec e_hi) {
  if (u_hi.n_elem > 1 && e_hi.n_elem > 1) {
    return arma::find(x >= u_hi + e_hi);
  } else if (u_hi.n_elem > 1) {
    return arma::find(x >= u_hi + e_hi[0]);
  } else if (e_hi.n_elem > 1) {
    return arma::find(x >= u_hi[0] + e_hi);
  } else {
    return arma::find(x >= u_hi[0] + e_hi[0]);
  }
}

arma::uvec find_low(const arma::vec x, const arma::vec u_lo, const arma::vec e_lo) {
  if (u_lo.n_elem > 1 && e_lo.n_elem > 1) {
    return arma::find(x <= u_lo - e_lo);
  } else if (u_lo.n_elem > 1) {
    return arma::find(x <= u_lo - e_lo[0]);
  } else if (e_lo.n_elem > 1) {
    return arma::find(x <= u_lo[0] - e_lo);
  } else {
    return arma::find(x <= u_lo[0] - e_lo[0]);
  }
}

void blend_transform(arma::vec& x, const arma::vec u_lo, const arma::vec e_lo, const arma::vec u_hi, const arma::vec e_hi) {
  if (e_lo.n_elem > 1 || e_lo[0] != 0.0) {
    // blend lower region, x < u_lo + e_lo
    arma::uvec i;
    if (u_lo.n_elem > 1 && e_lo.n_elem > 1) {
      i = arma::find(u_lo - e_lo < x && x < u_lo + e_lo);
      x(i) = 0.5 * (x(i) + u_lo(i) + e_lo(i)) - e_lo(i) / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_lo(i)) / e_lo(i));
    } else if (u_lo.n_elem > 1) {
      i = arma::find(u_lo - e_lo[0] < x && x < u_lo + e_lo[0]);
      x(i) = 0.5 * (x(i) + u_lo(i) + e_lo[0]) - e_lo[0] / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_lo(i)) / e_lo[0]);
    } else if (e_lo.n_elem > 1) {
      i = arma::find(u_lo[0] - e_lo < x && x < u_lo[0] + e_lo);
      x(i) = 0.5 * (x(i) + u_lo[0] + e_lo(i)) - e_lo(i) / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_lo[0]) / e_lo(i));
    } else {
      i = arma::find(u_lo[0] - e_lo[0] < x && x < u_lo[0] + e_lo[0]);
      x(i) = 0.5 * (x(i) + u_lo[0] + e_lo[0]) - e_lo[0] / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_lo[0]) / e_lo[0]);
    }
  }

  if (e_hi.n_elem > 1 || e_hi[0] != 0.0) {
    arma::uvec i;
    if (u_hi.n_elem > 1 && e_hi.n_elem > 1) {
      i = arma::find(u_hi - e_hi < x && x < u_hi + e_hi);
      x(i) = 0.5 * (x(i) + u_hi(i) - e_hi(i)) + e_hi(i) / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_hi(i)) / e_hi(i));
    } else if (u_hi.n_elem > 1) {
      i = arma::find(u_hi - e_hi[0] < x && x < u_hi + e_hi[0]);
      x(i) = 0.5 * (x(i) + u_hi(i) - e_hi[0]) + e_hi[0] / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_hi(i)) / e_hi[0]);
    } else if (e_hi.n_elem > 1) {
      i = arma::find(u_hi[0] - e_hi < x && x < u_hi[0] + e_hi);
      x(i) = 0.5 * (x(i) + u_hi[0] - e_hi(i)) + e_hi(i) / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_hi[0]) / e_hi(i));
    } else {
      i = arma::find(u_hi[0] - e_hi[0] < x && x < u_hi[0] + e_hi[0]);
      x(i) = 0.5 * (x(i) + u_hi[0] - e_hi[0]) + e_hi[0] / arma::datum::pi * cos(arma::datum::pi * 0.5 * (x(i) - u_hi[0]) / e_hi[0]);
    }
  }
}

arma::vec dblend_transform(const arma::vec x, const arma::vec u_lo, const arma::vec e_lo, const arma::vec u_hi, const arma::vec e_hi) {
  arma::vec dout(arma::size(x), arma::fill::ones);

  if (e_lo.n_elem > 1 || e_lo[0] != 0.0) {
    // blend lower region, x < u_lo + e_lo
    arma::uvec i;
    if (u_lo.n_elem > 1 && e_lo.n_elem > 1) {
      i = arma::find(u_lo - e_lo < x && x < u_lo + e_lo);
      dout(i) = 0.5 + 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_lo(i)) / e_lo(i));
    } else if (u_lo.n_elem > 1) {
      i = arma::find(u_lo - e_lo[0] < x && x < u_lo + e_lo[0]);
      dout(i) = 0.5 + 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_lo(i)) / e_lo[0]);
    } else if (e_lo.n_elem > 1) {
      i = arma::find(u_lo[0] - e_lo < x && x < u_lo[0] + e_lo);
      dout(i) = 0.5 + 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_lo[0]) / e_lo(i));
    } else {
      i = arma::find(u_lo[0] - e_lo[0] < x && x < u_lo[0] + e_lo[0]);
      dout(i) = 0.5 + 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_lo[0]) / e_lo[0]);
    }
  }

  if (e_hi.n_elem > 1 || e_hi[0] != 0.0) {
    arma::uvec i;
    if (u_hi.n_elem > 1 && e_hi.n_elem > 1) {
      i = arma::find(u_hi - e_hi < x && x < u_hi + e_hi);
      dout(i) = 0.5 - 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_hi(i)) / e_hi(i));
    } else if (u_hi.n_elem > 1) {
      i = arma::find(u_hi - e_hi[0] < x && x < u_hi + e_hi[0]);
      dout(i) = 0.5 - 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_hi(i)) / e_hi[0]);
    } else if (e_hi.n_elem > 1) {
      i = arma::find(u_hi[0] - e_hi < x && x < u_hi[0] + e_hi);
      dout(i) = 0.5 - 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_hi[0]) / e_hi(i));
    } else {
      i = arma::find(u_hi[0] - e_hi[0] < x && x < u_hi[0] + e_hi[0]);
      dout(i) = 0.5 - 0.5 * sin(arma::datum::pi * 0.5 * (x(i) - u_hi[0]) / e_hi[0]);
    }
  }

  return dout;
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
  if (arma::any(is_discrete == 1) && arma::any(is_discrete == 0)) {
    // mixed type => need to conditionally zero continuous densities
    arma::uvec cont_components = arma::find(is_discrete == 0);
    arma::uvec disc_points = arma::find(arma::any(compdens.cols(arma::find(is_discrete == 1)), 1));
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
  int n = std::max(std::max(x.n_elem, num_observations(probs)), std::max(scale.n_elem, num_observations(shapes)));
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
  int n = std::max(std::max(q.n_elem, num_observations(probs)), std::max(scale.n_elem, num_observations(shapes)));
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

// efficiently compute erlang mixture interval probabilities with possibly fixed parameters

template <typename TP, typename TS>
arma::vec dist_erlangmix_iprobability_impl(const arma::vec qmin, const arma::vec qmax, bool log_p, const TP probs, const arma::vec scale, const TS shapes) {
  int k = num_components(probs);
  int n = std::max(std::max(qmin.n_elem, qmax.n_elem), std::max(num_observations(probs), std::max(scale.n_elem, num_observations(shapes))));
  bool shape_is_matrix = is_matrix(shapes);
  arma::mat compprob(n, k);
  int i_qmin = 0, d_qmin = qmin.n_elem > 1 ? 1 : 0;
  int i_qmax = 0, d_qmax = qmax.n_elem > 1 ? 1 : 0;
  int i_s = 0, d_s = scale.n_elem > 1 ? 1 : 0;
  double curr_shape;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      if (shape_is_matrix) {
        curr_shape = shapes(i, j);
      } else {
        curr_shape = shapes[j];
      }
      compprob(i, j) = R::pgamma(qmax[i_qmax], curr_shape, scale[i_s], 1, 0) - R::pgamma(qmin[i_qmin], curr_shape, scale[i_s], 1, 0);
    }
    i_qmin += d_qmin;
    i_qmax += d_qmax;
    i_s += d_s;
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_free(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p) {
  int k = (params.n_cols - 1) / 2;
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, params.tail_cols(k), params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_shape(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec shapes) {
  int k = shapes.n_elem;
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, params.tail_cols(k), params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_scale(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec scale) {
  int k = params.n_cols / 2;
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, params.tail_cols(k), scale, params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_probs(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec probs) {
  int k = probs.n_elem;
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, probs, params.col(k), params.head_cols(k));
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_probs_scale(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec probs, const arma::vec scale) {
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, probs, scale, params);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_probs_shape(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec probs, const arma::vec shapes) {
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, probs, params.col(0), shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_scale_shape(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, params, scale, shapes);
}

// [[Rcpp::export]]
arma::vec dist_erlangmix_iprobability_fixed_probs_scale_shape(const arma::vec qmin, const arma::vec qmax, bool log_p, const arma::vec probs, const arma::vec scale, const arma::vec shapes) {
  return dist_erlangmix_iprobability_impl(qmin, qmax, log_p, probs, scale, shapes);
}

// efficiently compute blended densities with possibly fixed parameters

template <typename TP, typename TU, typename TE>
arma::vec dist_blended_density_impl(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const TP probs, const TU breaks, const TE epsilons) {
  int k = comp_densities.size();
  int n = std::max(std::max(x.n_elem, params.n_rows), std::max(num_observations(probs), std::max(num_observations(breaks), num_observations(epsilons))));
  int i = 0;
  arma::mat compdens(n, k, arma::fill::zeros);
  SEXP curr_params;
  bool breaks_is_matrix = is_matrix(breaks);
  bool epsilons_is_matrix = is_matrix(epsilons);
  for (int j = 0; j < k; j++) {
    arma::uvec curr_relevant;
    arma::vec curr_xblend;
    arma::vec curr_dblend;
    arma::vec curr_u_high, curr_e_high, curr_u_low, curr_e_low;
    if (j == 0) {
      curr_u_high = column_or_element(breaks, j);
      curr_e_high = column_or_element(epsilons, j);
      curr_u_low = arma::vec::fixed<1>{-std::numeric_limits<double>::infinity()};
      curr_e_low = arma::zeros<arma::vec>(1);
      curr_relevant = find_relevant(x, x, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) curr_u_high = curr_u_high(curr_relevant);
      if (epsilons_is_matrix) curr_e_high = curr_e_high(curr_relevant);
    } else if (j == k - 1) {
      curr_u_high = arma::vec::fixed<1>{std::numeric_limits<double>::infinity()};
      curr_e_high = arma::zeros<arma::vec>(1);
      curr_u_low = column_or_element(breaks, j - 1);
      curr_e_low = column_or_element(epsilons, j - 1);
      curr_relevant = find_relevant(x, x, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) curr_u_low = curr_u_low(curr_relevant);
      if (epsilons_is_matrix) curr_e_low = curr_e_low(curr_relevant);
    } else {
      curr_u_low = column_or_element(breaks, j - 1);
      curr_e_low = column_or_element(epsilons, j - 1);
      curr_u_high = column_or_element(breaks, j);
      curr_e_high = column_or_element(epsilons, j);
      curr_relevant = find_relevant(x, x, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) {
        curr_u_low = curr_u_low(curr_relevant);
        curr_u_high = curr_u_high(curr_relevant);
      }
      if (epsilons_is_matrix) {
        curr_e_low = curr_e_low(curr_relevant);
        curr_e_high = curr_e_high(curr_relevant);
      }
    }
    if (curr_relevant.n_elem == 0) continue;
    curr_xblend = x.elem(curr_relevant);
    curr_dblend = dblend_transform(curr_xblend, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
    blend_transform(curr_xblend, curr_u_low, curr_e_low, curr_u_high, curr_e_high);

    if (param_sizes[j] > 0) {
      arma::mat curr_params_mat = params.submat(curr_relevant, arma::regspace<arma::uvec>(i, i + param_sizes[j] - 1));
      curr_params = Rcpp::wrap(curr_params_mat);
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }

    Rcpp::Function curr_dfun = comp_densities[j];
    Rcpp::Function curr_ipfun = comp_iprobabilities[j];

    arma::vec curr_dens(n, arma::fill::zeros);
    arma::vec ptrunc = Rcpp::as<arma::vec>(curr_ipfun(curr_u_low, curr_u_high, curr_params, false));
    if (ptrunc.n_elem > 1) {
      curr_dens(curr_relevant) = Rcpp::as<arma::vec>(curr_dfun(curr_xblend, curr_params, false)) % curr_dblend / ptrunc;
    } else {
      curr_dens(curr_relevant) = Rcpp::as<arma::vec>(curr_dfun(curr_xblend, curr_params, false)) % curr_dblend / ptrunc[0];
    }

    compdens.col(j) = curr_dens;
  }
  if (arma::any(is_discrete == 1) && arma::any(is_discrete == 0)) {
    // mixed type => need to conditionally zero continuous densities
    arma::uvec cont_components = arma::find(is_discrete == 0);
    arma::uvec disc_points = arma::find(arma::any(compdens.cols(arma::find(is_discrete == 1)), 1));
    compdens.submat(disc_points, cont_components).fill(0.0);
  }
  arma::vec res = aggregate_mixture(compdens, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_blended_density_free(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, params.tail_cols(k), params.cols(arma::span(params.n_cols - 3 * k + 2, params.n_cols - 2 * k)), params.cols(arma::span(params.n_cols - 2 * k + 1, params.n_cols - k - 1)));
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_probs(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec probs) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, probs, params.cols(arma::span(params.n_cols - 2 * k + 2, params.n_cols - k)), params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_breaks(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec breaks) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, params.tail_cols(k), breaks, params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_eps(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec epsilons) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, params.tail_cols(k), params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_probs_breaks(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec probs, const arma::vec breaks) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, probs, breaks, params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_probs_eps(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec probs, const arma::vec epsilons) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, probs, params.tail_cols(k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_breaks_eps(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec breaks, const arma::vec epsilons) {
  int k = comp_densities.size();
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, params.tail_cols(k), breaks, epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_density_fixed_probs_breaks_eps(const arma::vec x, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_densities, const Rcpp::List comp_iprobabilities, const arma::uvec is_discrete, const arma::vec probs, const arma::vec breaks, const arma::vec epsilons) {
  return dist_blended_density_impl(x, params, log_p, param_sizes, comp_densities, comp_iprobabilities, is_discrete, probs, breaks, epsilons);
}

// efficiently compute blended probabilities with possibly fixed parameters

template <typename TP, typename TU, typename TE>
arma::vec dist_blended_probability_impl(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const TP probs, const TU breaks, const TE epsilons) {
  int k = comp_iprobabilities.size();
  int n = std::max(std::max(q.n_elem, params.n_rows), std::max(num_observations(probs), std::max(num_observations(breaks), num_observations(epsilons))));
  int i = 0;
  arma::mat compprob(n, k, arma::fill::zeros);
  SEXP curr_params;
  for (int j = 0; j < k; j++) {
    arma::vec curr_u_high, curr_e_high, curr_u_low, curr_e_low;
    if (j == 0) {
      curr_u_low = arma::vec::fixed<1>{-std::numeric_limits<double>::infinity()};
      curr_e_low = arma::zeros<arma::vec>(1);
    } else {
      curr_u_low = column_or_element(breaks, j - 1);
      curr_e_low = column_or_element(epsilons, j - 1);
    }
    if (j == k - 1) {
      curr_u_high = arma::vec::fixed<1>{std::numeric_limits<double>::infinity()};
      curr_e_high = arma::zeros<arma::vec>(1);
    } else {
      curr_u_high = column_or_element(breaks, j);
      curr_e_high = column_or_element(epsilons, j);
    }

    arma::vec curr_prob(n, arma::fill::zeros);
    if (lower_tail) {
      curr_prob(find_high(q, curr_u_high, curr_e_high)).fill(1.0);
    } else {
      curr_prob(find_low(q, curr_u_low, curr_e_low)).fill(1.0);
    }
    arma::uvec curr_relevant = find_relevant(q, q, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
    if (curr_relevant.n_elem == 0) {
      compprob.col(j) = curr_prob;
      continue;
    }
    arma::vec curr_xblend = q.elem(curr_relevant);
    if (curr_u_low.n_elem > 1) {
      curr_u_low = curr_u_low(curr_relevant);
    }
    if (curr_e_low.n_elem > 1) {
      curr_e_low = curr_e_low(curr_relevant);
    }
    if (curr_u_high.n_elem > 1) {
      curr_u_high = curr_u_high(curr_relevant);
    }
    if (curr_e_high.n_elem > 1) {
      curr_e_high = curr_e_high(curr_relevant);
    }
    blend_transform(curr_xblend, curr_u_low, curr_e_low, curr_u_high, curr_e_high);

    if (param_sizes[j] > 0) {
      arma::mat curr_params_mat = params.submat(curr_relevant, arma::regspace<arma::uvec>(i, i + param_sizes[j] - 1));
      curr_params = Rcpp::wrap(curr_params_mat);
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }

    Rcpp::Function curr_ipfun = comp_iprobabilities[j];

    arma::vec ptrunc = Rcpp::as<arma::vec>(curr_ipfun(curr_u_low, curr_u_high, curr_params, false));
    arma::vec pobs;
    if (lower_tail) {
      pobs = Rcpp::as<arma::vec>(curr_ipfun(curr_u_low, curr_xblend, curr_params, false));
    } else {
      pobs = Rcpp::as<arma::vec>(curr_ipfun(curr_xblend, curr_u_high, curr_params, false));
    }
    
    if (ptrunc.n_elem > 1) {
      curr_prob(curr_relevant) = pobs / ptrunc;
    } else {
      curr_prob(curr_relevant) = pobs / ptrunc[0];
    }

    compprob.col(j) = curr_prob;
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_free(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), params.cols(arma::span(params.n_cols - 3 * k + 2, params.n_cols - 2 * k)), params.cols(arma::span(params.n_cols - 2 * k + 1, params.n_cols - k - 1)));
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_probs(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, probs, params.cols(arma::span(params.n_cols - 2 * k + 2, params.n_cols - k)), params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_breaks(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec breaks) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), breaks, params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_eps(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_probs_breaks(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec breaks) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, probs, breaks, params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_probs_eps(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, probs, params.tail_cols(k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_breaks_eps(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec breaks, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), breaks, epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_probability_fixed_probs_breaks_eps(const arma::vec q, const arma::mat params, bool lower_tail, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec breaks, const arma::vec epsilons) {
  return dist_blended_probability_impl(q, params, lower_tail, log_p, param_sizes, comp_iprobabilities, probs, breaks, epsilons);
}

// efficiently compute blended interval probabilities with possibly fixed parameters

template <typename TP, typename TU, typename TE>
arma::vec dist_blended_iprobability_impl(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const TP probs, const TU breaks, const TE epsilons) {
  int k = comp_iprobabilities.size();
  int n = std::max(std::max(qmin.n_elem, qmax.n_elem), std::max(std::max(params.n_rows, num_observations(probs)), std::max(num_observations(breaks), num_observations(epsilons))));
  int i = 0;
  arma::mat compprob(n, k, arma::fill::zeros);
  SEXP curr_params;
  bool breaks_is_matrix = is_matrix(breaks);
  bool epsilons_is_matrix = is_matrix(epsilons);
  for (int j = 0; j < k; j++) {
    arma::uvec curr_relevant;
    arma::vec curr_qminblend, curr_qmaxblend;
    arma::vec curr_u_high, curr_e_high, curr_u_low, curr_e_low;
    if (j == 0) {
      curr_u_high = column_or_element(breaks, j);
      curr_e_high = column_or_element(epsilons, j);
      curr_u_low = arma::vec::fixed<1>{-std::numeric_limits<double>::infinity()};
      curr_e_low = arma::zeros<arma::vec>(1);
      curr_relevant = find_relevant(qmin, qmax, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) curr_u_high = curr_u_high(curr_relevant);
      if (epsilons_is_matrix) curr_e_high = curr_e_high(curr_relevant);
    } else if (j == k - 1) {
      curr_u_high = arma::vec::fixed<1>{std::numeric_limits<double>::infinity()};
      curr_e_high = arma::zeros<arma::vec>(1);
      curr_u_low = column_or_element(breaks, j - 1);
      curr_e_low = column_or_element(epsilons, j - 1);
      curr_relevant = find_relevant(qmin, qmax, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) curr_u_low = curr_u_low(curr_relevant);
      if (epsilons_is_matrix) curr_e_low = curr_e_low(curr_relevant);
    } else {
      curr_u_low = column_or_element(breaks, j - 1);
      curr_e_low = column_or_element(epsilons, j - 1);
      curr_u_high = column_or_element(breaks, j);
      curr_e_high = column_or_element(epsilons, j);
      curr_relevant = find_relevant(qmin, qmax, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
      if (breaks_is_matrix) {
        curr_u_low = curr_u_low(curr_relevant);
        curr_u_high = curr_u_high(curr_relevant);
      }
      if (epsilons_is_matrix) {
        curr_e_low = curr_e_low(curr_relevant);
        curr_e_high = curr_e_high(curr_relevant);
      }
    }
    if (curr_relevant.n_elem == 0) continue;
    curr_qminblend = qmin.elem(curr_relevant);
    curr_qmaxblend = qmax.elem(curr_relevant);
    blend_transform(curr_qminblend, curr_u_low, curr_e_low, curr_u_high, curr_e_high);
    blend_transform(curr_qmaxblend, curr_u_low, curr_e_low, curr_u_high, curr_e_high);

    if (curr_u_low.n_elem > 1) {
      curr_qminblend.elem(find_low(qmin.elem(curr_relevant), curr_u_low, curr_e_low)) = curr_u_low;
    } else {
      curr_qminblend.elem(find_low(qmin.elem(curr_relevant), curr_u_low, curr_e_low)).fill(curr_u_low[0]);
    }
    if (curr_u_high.n_elem > 1) {
      curr_qmaxblend.elem(find_high(qmax.elem(curr_relevant), curr_u_high, curr_e_high)) = curr_u_high;
    } else {
      curr_qmaxblend.elem(find_high(qmax.elem(curr_relevant), curr_u_high, curr_e_high)).fill(curr_u_high[0]);
    }

    if (param_sizes[j] > 0) {
      arma::mat curr_params_mat = params.submat(curr_relevant, arma::regspace<arma::uvec>(i, i + param_sizes[j] - 1));
      curr_params = Rcpp::wrap(curr_params_mat);
      i += param_sizes[j];
    } else {
      curr_params = R_NilValue;
    }

    Rcpp::Function curr_ipfun = comp_iprobabilities[j];

    arma::vec curr_prob(n, arma::fill::zeros);
    arma::vec ptrunc = Rcpp::as<arma::vec>(curr_ipfun(curr_u_low, curr_u_high, curr_params, false));
    arma::vec pobs = Rcpp::as<arma::vec>(curr_ipfun(curr_qminblend, curr_qmaxblend, curr_params, false));
    
    if (ptrunc.n_elem > 1) {
      curr_prob(curr_relevant) = pobs / ptrunc;
    } else {
      curr_prob(curr_relevant) = pobs / ptrunc[0];
    }

    compprob.col(j) = curr_prob;
  }
  arma::vec res = aggregate_mixture(compprob, probs);
  if (log_p) res = log(res);
  return res;
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_free(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), params.cols(arma::span(params.n_cols - 3 * k + 2, params.n_cols - 2 * k)), params.cols(arma::span(params.n_cols - 2 * k + 1, params.n_cols - k - 1)));
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_probs(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, probs, params.cols(arma::span(params.n_cols - 2 * k + 2, params.n_cols - k)), params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_breaks(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec breaks) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), breaks, params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_eps(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), params.cols(params.n_cols - 2 * k + 1, params.n_cols - k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_probs_breaks(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec breaks) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, probs, breaks, params.tail_cols(k - 1));
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_probs_eps(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, probs, params.tail_cols(k - 1), epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_breaks_eps(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec breaks, const arma::vec epsilons) {
  int k = comp_iprobabilities.size();
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, params.tail_cols(k), breaks, epsilons);
}

// [[Rcpp::export]]
arma::vec dist_blended_iprobability_fixed_probs_breaks_eps(const arma::vec qmin, const arma::vec qmax, const arma::mat params, bool log_p, const arma::uvec param_sizes, const Rcpp::List comp_iprobabilities, const arma::vec probs, const arma::vec breaks, const arma::vec epsilons) {
  return dist_blended_iprobability_impl(qmin, qmax, params, log_p, param_sizes, comp_iprobabilities, probs, breaks, epsilons);
}

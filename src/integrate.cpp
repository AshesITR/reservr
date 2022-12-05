// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <boost/heap/priority_queue.hpp>
#include <RcppArmadillo.h>
#include <Rcpp.h>

struct Bounds {
  double lower, upper, value, error;

  Bounds() : lower(0), upper(0), value(0), error(0) {}

  Bounds(double lower, double upper, double value, double error) : lower(lower), upper(upper),
    value(value), error(error) {}

  bool operator<(const Bounds b) const {
    return this->error < b.error;
  }
};

static const Rcpp::Function asNamespace = Rcpp::Environment::base_env()["asNamespace"];
static const Rcpp::Environment pkg_namespace = asNamespace("reservr");
static const Rcpp::Function pick_params_at_idx = pkg_namespace["pick_params_at_idx"];

static const arma::mat::fixed<15, 2> gk_weights{
  // GK15 weights
  0.022935322010529, 0.022935322010529, 0.063092092629979, 0.063092092629979, 0.104790010322250, 0.104790010322250,
  0.140653259715525, 0.140653259715525, 0.169004726639267, 0.169004726639267, 0.190350578064785, 0.190350578064785,
  0.204432940075298, 0.204432940075298, 0.209482141084728,
  // GK7 weights
  0, 0, 0.129484966168870, 0.129484966168870, 0, 0, 0.279705391489277, 0.279705391489277, 0, 0, 0.381830050505119,
  0.381830050505119, 0, 0, 0.417959183673469
};
static const arma::rowvec::fixed<15> gk_nodes{
  0.991455371120813, -0.991455371120813, 0.949107912342759, -0.949107912342759, 0.864864423359769, -0.864864423359769,
  0.741531185599394, -0.741531185599394, 0.586087235467691, -0.586087235467691, 0.405845151377397, -0.405845151377397,
  0.207784955007898, -0.207784955007898, 0.000000000000000
};

template <typename T>
void integrate_gk_step(const Rcpp::Function fun, const arma::vec& lower, const arma::vec& upper,
  const T params, const std::vector<unsigned int>& indices, std::vector<boost::heap::priority_queue<Bounds>>& pq) {

  arma::vec midpoint = 0.5 * (lower + upper);
  arma::vec radius = 0.5 * (upper - lower);

  arma::mat f_eval = midpoint * arma::ones(1, 15) + radius * gk_nodes;
  f_eval.each_col( [fun, params, radius](arma::vec& col) { col = Rcpp::as<arma::vec>(fun(col, params)) % radius; });

  arma::mat estimates = f_eval * gk_weights;
  arma::vec values = estimates.col(0);
  arma::vec errors = abs(estimates.col(0) - estimates.col(1));

  if (lower.n_elem < pq.size()) {
    for (std::size_t i = 0; i < lower.n_elem; i++) {
      pq[indices[i]].push(Bounds(lower(i), upper(i), values(i), errors(i)));
    }
  } else {
    for (std::size_t i = 0; i < lower.n_elem; i++) {
      pq[i].push(Bounds(lower(i), upper(i), values(i), errors(i)));
    }
  }
}

template <typename T>
T params_at(const T& params, const std::vector<unsigned int> indexes);

template <>
Rcpp::List params_at(const Rcpp::List& params, const std::vector<unsigned int> indexes) {
  return pick_params_at_idx(params, indexes);
}

template <>
arma::mat params_at(const arma::mat& params, const std::vector<unsigned int> indexes) {
  return params.rows(arma::uvec(indexes));
}

template <typename T>
arma::vec integrate_impl(const Rcpp::Function fun, const arma::vec& lower, const arma::vec& upper,
  const T params, double tolerance, int max_iter, arma::vec& res, arma::cube& info) {

  std::vector<boost::heap::priority_queue<Bounds>> bounds_pq(lower.n_elem);
  bool all_converged = false, some_converged = false;
  std::vector<bool> converged(lower.n_elem);

  int iter = 1;

  // 1. Initial integral estimate
  integrate_gk_step(fun, lower, upper, params, std::vector<unsigned int>(), bounds_pq);
  T curr_params;

  // 2. While not converged (and max_iter not reached) ...
  while (iter < max_iter && !all_converged) {
    iter++;
    all_converged = true;
    some_converged = false;
    std::vector<unsigned int> non_converged;
    // find out which lower[i], upper[i] integrals aren't done yet ...
    for (std::size_t i = 0; i < lower.n_elem; i++) {
      if (!converged[i]) {
        double curr_error = 0; // TODO std::accumulate(bounds_pq[i].begin(), bounds_pq[i].end()).error;
        for (Bounds bd : bounds_pq[i]) {
          curr_error += bd.error;
        }
        if (curr_error < tolerance) {
          converged[i] = true;
          some_converged = true;
        } else {
          non_converged.push_back(i);
          all_converged = false;
        }
      }
    }

    if (!all_converged) {
      if (some_converged) {
        // only re-pick params if some integrals converged in the last iteration
        curr_params = params_at(params, non_converged);
      } else if (non_converged.size() == lower.n_elem) {
        curr_params = params;
      }
      arma::vec curr_lower(non_converged.size());
      arma::vec curr_upper(non_converged.size());

      // Pop the most inaccurate interval of each lower[i], upper[i] ...
      for (std::size_t i = 0; i < non_converged.size(); i++) {
        Bounds curr = bounds_pq[non_converged[i]].top();
        bounds_pq[non_converged[i]].pop();
        curr_lower[i] = curr.lower;
        curr_upper[i] = curr.upper;
      }

      arma::vec curr_midpoint = 0.5 * (curr_lower + curr_upper);

      // ... and bisect it.
      integrate_gk_step(fun, curr_lower, curr_midpoint, curr_params, non_converged, bounds_pq);
      integrate_gk_step(fun, curr_midpoint, curr_upper, curr_params, non_converged, bounds_pq);
    }
  }

  // 3. aggregate all integral values
  for (std::size_t i = 0; i < lower.n_elem; i++) {
    // TODO res[i] = std::accumulate(bounds_pq[i].begin(), bounds_pq[i].end()).value;
    int j = 0;
    for (Bounds bd : bounds_pq[i]) {
      res[i] += bd.value;
      if (info.n_elem > 0) {
        info(i, j, 0) = bd.value;
        info(i, j, 1) = bd.error;
        info(i, j, 2) = bd.lower;
        info(i, j, 3) = bd.upper;
      }
      j++;
    }
  }

  if (!all_converged) {
    Rcpp::warning("`.max_iter` reached. Some integrals did not converge.");
  }

  return res;
}

template <typename T>
Rcpp::List do_integrate_gk(const Rcpp::Function fun, const arma::vec& lower, const arma::vec& upper,
  const T& params, const double tolerance, const int max_iter, bool debug) {
  arma::cube* info;
  if (debug) {
    info = new arma::cube(lower.n_elem, max_iter, 4, arma::fill::zeros);
  } else {
    // empty cube to save memory allocation
    info = new arma::cube();
  }
  arma::vec res(lower.n_elem, arma::fill::zeros);
  integrate_impl(fun, lower, upper, params, tolerance, max_iter, res, *info);

  Rcpp::NumericVector r_res = Rcpp::wrap(res);
  r_res.attr("dim") = R_NilValue;

  // Return as vector instead of colvec.
  if (debug) {
    return Rcpp::List::create(
      Rcpp::Named("value") = r_res,
      Rcpp::Named("info") = *info
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("value") = r_res
    );
  }
}

// [[Rcpp::export]]
Rcpp::List do_integrate_gk_lst(const Rcpp::Function fun, const arma::vec& lower, const arma::vec& upper,
  const Rcpp::List& params, const double tolerance, const int max_iter, bool debug) {
  return do_integrate_gk(fun, lower, upper, params, tolerance, max_iter, debug);
}

// [[Rcpp::export]]
Rcpp::List do_integrate_gk_mat(const Rcpp::Function fun, const arma::vec& lower, const arma::vec& upper,
  const arma::mat& params, const double tolerance, const int max_iter, bool debug) {
  return do_integrate_gk(fun, lower, upper, params, tolerance, max_iter, debug);
}

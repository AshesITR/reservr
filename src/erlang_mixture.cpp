// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>

const double cens_gamma_logdens(const double xmin, const double xmax, const double shape, const double scale) {
  if (xmin <= 0.0) {
    return R::pgamma(xmax, shape, scale, 1, 1);
  } else if (xmax == 1.0 / 0.0) {
    return R::pgamma(xmin, shape, scale, 0, 1);
  } else {
    return log(R::pgamma(xmax, shape, scale, 1, 0) - R::pgamma(xmin, shape, scale, 1, 0));
  }
}

const double cens_gamma_loggrad(const double xmin, const double xmax, const double shape, const double scale) {
  if (xmin <= 0.0 && xmax == 1.0 / 0.0) {
    return 0.0;
  } else if (xmin <= 0.0) {
    return -xmax * R::dgamma(xmax, shape, scale, 0) / scale / R::pgamma(xmax, shape, scale, 1, 0);
  } else if (xmax == 1.0 / 0.0) {
    return xmin * R::dgamma(xmin, shape, scale, 0) / scale / R::pgamma(xmin, shape, scale, 0, 0);
  } else {
    return (xmin * R::dgamma(xmin, shape, scale, 0) - xmax * R::dgamma(xmax, shape, scale, 0)) / scale /
      (R::pgamma(xmax, shape, scale, 1, 0) - R::pgamma(xmin, shape, scale, 1, 0));
  }
}

struct EllikWorker : public RcppParallel::Worker {
  const arma::vec *xmin, *xmax, *tmin, *tmax, *weight, *shapes;
  const double scale;
  const arma::mat *zadj;
  arma::mat *logdens, *logdens_grad_mat;
  arma::vec *logdens_grad_vec;

  EllikWorker(const arma::vec *xmin, const arma::vec *xmax, const arma::vec *tmin, const arma::vec *tmax,
    const arma::vec *weight, const arma::vec *shapes, const double scale, const arma::mat *zadj,
    arma::mat *logdens, arma::mat *logdens_grad_mat, arma::vec *logdens_grad_vec) :
    xmin(xmin), xmax(xmax), tmin(tmin), tmax(tmax), weight(weight), shapes(shapes), scale(scale), zadj(zadj),
    logdens(logdens), logdens_grad_mat(logdens_grad_mat), logdens_grad_vec(logdens_grad_vec) {}

  EllikWorker(EllikWorker& base, RcppParallel::Split) :
    xmin(base.xmin), xmax(base.xmax), tmin(base.tmin), tmax(base.tmax), weight(base.weight), shapes(base.shapes),
    scale(base.scale), zadj(base.zadj), logdens(base.logdens), logdens_grad_mat(base.logdens_grad_mat),
    logdens_grad_vec(base.logdens_grad_vec) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Iterate over xmin, xmax, tmin, tmax from begin to end
    for (std::size_t i = begin; i != end; i++) {
      const bool is_censored = (*xmin)(i) != (*xmax)(i);
      for (std::size_t j = 0; j < shapes->n_elem; j++) {
        if ((*zadj)(i, j) > 0.0) {
          if (is_censored) {
            (*logdens)(i, j) = cens_gamma_logdens((*xmin)(i), (*xmax)(i), (*shapes)(j), scale) -
              cens_gamma_logdens((*tmin)(i), (*tmax)(i), (*shapes)(j), scale);
            (*logdens_grad_mat)(i, j) = cens_gamma_loggrad((*xmin)(i), (*xmax)(i), (*shapes)(j), scale) -
              cens_gamma_loggrad((*tmin)(i), (*tmax)(i), (*shapes)(j), scale);
          } else {
            (*logdens)(i, j) = R::dgamma((*xmin)(i), (*shapes)(j), scale, 1) -
              cens_gamma_logdens((*tmin)(i), (*tmax)(i), (*shapes)(j), scale);
            (*logdens_grad_mat)(i, j) = -(*shapes)(j) / scale -
              cens_gamma_loggrad((*tmin)(i), (*tmax)(i), (*shapes)(j), scale);
          }
        }
        if (!is_censored) {
          (*logdens_grad_vec)(i) = -(*xmin)(i) * (*weight)(i) / scale / scale;
        }
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List trunc_erlangmix_ellik(const arma::vec& xmin, const arma::vec& xmax, const arma::vec& tmin,
  const arma::vec& tmax, const arma::vec& weight, const arma::vec& shapes, double scale, const arma::mat& zadj,
  bool parallel) {

  if (xmin.n_elem == 0) {
    arma::vec empty = arma::zeros(0);
    return Rcpp::List::create(
      Rcpp::Named("objective") = empty,
      Rcpp::Named("gradient") = empty
    );
  }

  arma::mat logdens(arma::size(zadj));
  arma::mat logdens_grad_mat(arma::size(zadj));
  arma::vec logdens_grad_vec(xmin.n_elem);

  EllikWorker worker(&xmin, &xmax, &tmin, &tmax, &weight, &shapes, scale, &zadj,
    &logdens, &logdens_grad_mat, &logdens_grad_vec);

  if (parallel) {
    RcppParallel::parallelFor(0, xmin.n_elem, worker);
  } else {
    worker(0, xmin.size());
  }

  return Rcpp::List::create(
    Rcpp::Named("objective") = accu(logdens % zadj),
    Rcpp::Named("gradient") = accu(logdens_grad_mat % zadj) + sum(logdens_grad_vec)
  );
}

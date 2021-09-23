#' Fit a general distribution to observations
#'
#' The default implementation performs maximum likelihood estimation on all
#' placeholder parameters.
#'
#' For Erlang mixture distributions and for Mixture distributions, an
#' EM-Algorithm is instead used to improve stability.
#'
#' @details
#' `fit()` and `fit_dist()` will chose an optimisation method optimized for the specific distribution given.
#' `fit_dist_direct()` can be used to force direct maximisation of the likelihood.
#'
#' @param dist A `Distribution` object.
#' @param obs Set of observations as produced by [trunc_obs()] or convertible
#' via [as_trunc_obs()].
#' @param start Initial values of all placeholder parameters.
#' If missing, starting values are obtained from [fit_dist_start()].
#' @param ... Distribution-specific arguments for the fitting procedure
#'
#' @return A list with at least the elements
#'
#'   * `params` the fitted parameters in the same structure as `init`.
#'   * `logLik` the final log-likelihood
#'
#' Additional information may be provided depending on `dist`.
#'
#' @export
#' @family distribution fitting functions
#'
#' @examples
#' x <- rexp(100)
#' lambda_hat <- 1 / mean(x)
#' lambda_hat2 <- fit_dist(dist_exponential(), x)$params$rate
#' identical(lambda_hat, lambda_hat2)
fit_dist <- function(dist, obs, start, ...) {
  UseMethod("fit_dist")
}

#' Find starting values for distribution parameters
#'
#' @param dist A `Distribution` object.
#' @param obs Observations to fit to.
#' @param ... Additional arguments for the initialisation procedure
#'
#' @return A list of initial parameters suitable for passing to [fit_dist()].
#' @export
#'
#' @examples
#' fit_dist_start(dist_exponential(), rexp(100))
fit_dist_start <- function(dist, obs, ...) {
  UseMethod("fit_dist_start")
}

#' @param .start_with_default Before directly optimising the likelihood, use an optimised algorithm for finding
#' better starting values?
#'
#' @examples
#' dist <- dist_mixture(list(dist_normal(), dist_genpareto1(u = 6)))
#' params <- list(
#'   dists = list(list(mean = 5, sd = 1), list(sigmau = 1, xi = 0)), probs = list(0.95, 0.05)
#' )
#' u <- runif(100, 3, 10)
#' x <- dist$sample(100, with_params = params)
#' obs <- trunc_obs(x = x[x <= u], tmin = -Inf, tmax = u[x <= u])
#'
#' default_fit <- fit_dist(dist, obs)
#' direct_fit <- fit_dist_direct(dist, obs)
#' # NB: direct optimisation steps with pre-run take a few seconds
#' \donttest{
#' direct_fit_init <- fit_dist_direct(dist, obs, start = default_fit$params)
#' direct_fit_auto_init <- fit_dist_direct(dist, obs, .start_with_default = TRUE)
#'
#' stopifnot(direct_fit_init$logLik == direct_fit_auto_init$logLik)
#'
#' c(default_fit$logLik, direct_fit$logLik, direct_fit_init$logLik)
#' }
#'
#' @export
#' @family distribution fitting functions
#' @rdname fit_dist
fit_dist_direct <- function(dist, obs, start, ..., .start_with_default = FALSE) {
  if (.start_with_default) {
    start <- fit_dist(dist, obs, start, ...)$params
  }
  fit_dist.Distribution(dist, obs, start, ...)
}

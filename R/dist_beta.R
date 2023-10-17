#' Beta Distribution
#'
#' See [stats::Beta]
#'
#' All parameters can be overridden with
#' `with_params = list(shape = ..., scale = ...)`.
#'
#' @param shape1 First scalar shape parameter, or `NULL` as a placeholder.
#' @param shape2 Second scalar shape parameter, or `NULL` as a placeholder.
#' @param ncp Scalar non-centrality parameter, or `NULL` as a placeholder.
#'
#' @return A `BetaDistribution` object.
#' @export
#'
#' @examples
#' d_beta <- dist_beta(shape1 = 2, shape2 = 2, ncp = 0)
#' x <- d_beta$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_beta,
#'   estimated = d_beta,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "beta")$estimate
#'     )
#'   ),
#'   .x = seq(0, 2, length.out = 100)
#' )
#'
#' @family Distributions
dist_beta <- function(shape1 = NULL, shape2 = NULL, ncp = NULL) {
  BetaDistribution$new(shape1 = shape1, shape2 = shape2, ncp = ncp)
}

BetaDistribution <- distribution_class_simple(
  name = "Beta",
  fun_name = "beta",
  params = list(
    shape1 = I_POSITIVE_REALS,
    shape2 = I_POSITIVE_REALS,
    ncp = I_POSITIVE_REALS
  ),
  support = I_UNIT_INTERVAL, # TODO diff_*?
  tf_logdensity = function() function(x, args) { # nolint: brace.
    tf <- tensorflow::tf
    shape1 <- args[["shape1"]]
    shape2 <- args[["shape2"]]
    ncp <- args[["ncp"]]

    xbad <- x < 0.0 | x > 1.0 | ncp != 0.0
    x0 <- x == 0.0
    x_safe <- tf$where(x0 | xbad, K$one, x)
    lterm <- tf$where(shape1 < 1.0, K$inf, tf$where(shape1 > 1.0, K$neg_inf, 0.0))
    x1 <- x == 1.0
    xs_safe <- tf$where(x1 | xbad, K$one, 1.0 - x)
    rterm <- tf$where(shape2 < 1.0, K$inf, tf$where(shape2 > 1.0, K$neg_inf, 0.0))

    tf$where(
      xbad,
      K$neg_inf,
      -tf$math$lbeta(tf$stack(list(shape1, shape2), axis = -1L)) +
        tf$where(x0, lterm, (shape1 - 1.0) * log(x_safe)) +
        tf$where(x1, rterm, (shape2 - 1.0) * log(xs_safe))
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    tf <- tensorflow::tf
    shape1 <- args[["shape1"]]
    shape2 <- args[["shape2"]]
    ncp <- args[["ncp"]]

    qmin_safe <- tf$where(qmin < 0.0, K$zero, tf$where(qmin > 1.0, K$one, qmin))
    qmax_safe <- tf$where(qmax < 0.0, K$zero, tf$where(qmax > 1.0, K$one, qmax))

    nomass <- qmax <= 0.0 | qmin >= 1.0 | ncp != 0.0
    prob <- tf$math$betainc(shape1, shape2, qmax_safe) - tf$math$betainc(shape1, shape2, qmin_safe)
    prob_safe <- tf$where(nomass, K$one, prob)

    tf$where(nomass, K$neg_inf, log(prob_safe))
  }
)

#' @export
fit_dist_start.BetaDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w)
  if (length(ph)) {
    helper <- mom[1L] * (1 - mom[1L]) / mom[2L] - 1.0
  }

  if ("shape1" %in% ph) {
    res$shape1 <- mom[1L] * helper
  }

  if ("shape2" %in% ph) {
    res$shape2 <- (1.0 - mom[1L]) * helper
  }

  if ("ncp" %in% ph) {
    # TODO starting values for ncp?
    res$ncp <- 0.0
  }
  res
}

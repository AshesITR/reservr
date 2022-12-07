#' Negative binomial Distribution
#'
#' See [stats::NegBinomial]
#'
#' Both parameters can be overridden with
#' `with_params = list(size = ..., prob = ...)`.
#'
#' @param size Number of successful trials parameter, or `NULL` as a
#' placeholder. Non-integer values > 0 are allowed.
#' @param mu Mean parameter, or `NULL` as a placeholder.
#'
#' @return A `NegativeBinomialDistribution` object.
#' @export
#'
#' @examples
#' d_nbinom <- dist_negbinomial(size = 3.5, mu = 8.75)
#' x <- d_nbinom$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_nbinom,
#'   estimated = d_nbinom,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "nbinom")$estimate
#'     )
#'   ),
#'   .x = 0:max(x)
#' )
#'
#' @family Distributions
dist_negbinomial <- function(size = NULL, mu = NULL) {
  NegativeBinomialDistribution$new(size = size, mu = mu)
}

NegativeBinomialDistribution <- distribution_class_simple(
  name = "NegativeBinomial",
  fun_name = "nbinom",
  type = "discrete",
  params = list(
    size = I_POSITIVE_REALS,
    mu = I_POSITIVE_REALS
  ),
  support = I_NATURALS,
  diff_density = function(x, vars, log, params) {
    res <- vars
    if ("mu" %in% names(vars)) {
      res$mu <- if (log) {
        (x / params$mu - 1.0) * params$size / (params$size + params$mu)
      } else {
        (x / params$mu - 1.0) * params$size / (params$size + params$mu) *
          dnbinom(x = x, size = params$size, mu = params$mu)
      }
    }

    if ("size" %in% names(vars)) {
      log_diff <- digamma(x + params$size) - digamma(params$size) +
        log(params$size / (params$size + params$mu)) +
        (params$mu - x) / (params$size + params$mu)
      res$size <- if (log) {
        log_diff
      } else {
        log_diff * dnbinom(x, size = params$size, mu = params$mu)
      }
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    mu <- tf$broadcast_to(args[["mu"]], x$shape)
    size <- tf$broadcast_to(args[["size"]], x$shape)

    mu0 <- mu == K$zero
    prob <- mu / (size + mu)

    ok <- x >= K$zero
    x_safe <- tf$where(ok, x, K$zero)

    tf$where(
      mu0,
      tf$where(x == K$zero, K$zero, K$neg_inf),
      tf$where(
        ok,
        size * log1p(-prob) + x_safe * log(prob) -
          tf$math$lbeta(tf$stack(list(K$one + x_safe, size), axis = 1L)) - log(x_safe + size),
        K$neg_inf
      )
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    mu <- tf$broadcast_to(args[["mu"]], qmin$shape)
    size <- tf$broadcast_to(args[["size"]], qmin$shape)

    mu0 <- mu == K$zero
    prob <- mu / (size + mu)

    qmin0 <- qmin <= K$zero
    qmax_ok <- qmax >= K$zero & tf$math$is_finite(qmax)
    qmax_nok <- tf$where(qmax < K$zero, K$neg_inf, K$zero)

    qmin_safe <- tf$math$maximum(K$zero, tf$math$ceil(qmin)) - K$one
    qmax_safe <- tf$math$maximum(K$zero, tf$math$floor(qmax))

    tf$where(
      mu0,
      tf$where(qmin0, qmax_nok, K$neg_inf),
      tf$where(
        qmin0,
        tf$where(
          qmax_ok,
          log(tf$math$betainc(size, qmax_safe + K$one, K$one - prob)),
          qmax_nok
        ),
        tf$where(
          qmax_ok,
          log(tf$math$betainc(qmin_safe + K$one, size, prob) - tf$math$betainc(qmax_safe + K$one, size, prob)),
          log(tf$math$betainc(qmin_safe + K$one, size, prob))
        )
      )
    )
  },
  tf_is_discrete_at = function() function(x, args) { # nolint: brace.
    mu <- tf$broadcast_to(args[["mu"]], x$shape)
    tf$where(
      mu == K$zero,
      x == K$zero,
      tf_is_integerish(x) & x >= K$zero
    )
  }
)

#' @export
fit_dist_start.NegativeBinomialDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0L)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 2L)
  if ("mu" %in% ph) {
    res$mu <- mom[1L]
  }
  if ("size" %in% ph) {
    res$size <- if (mom[2L] > mom[1L]) {
      mom[1L]^2.0 / (mom[2L] - mom[1L])
    } else {
      100.0
    }
  }
  res
}

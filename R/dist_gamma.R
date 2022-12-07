#' Gamma distribution
#'
#' See [stats::GammaDist].
#'
#' Both parameters can be overridden with
#' `with_params = list(shape = ..., rate = ...)`.
#'
#' @param shape Scalar shape parameter, or `NULL` as a placeholder.
#' @param rate Scalar rate parameter, or `NULL` as a placeholder.
#'
#' @return A `GammaDistribution` object.
#' @export
#'
#' @examples
#' alpha <- 2
#' beta <- 2
#'
#' d_gamma <- dist_gamma(shape = alpha, rate = beta)
#' x <- d_gamma$sample(100)
#' d_emp <- dist_empirical(x, positive = TRUE)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_gamma,
#'   estimated = d_gamma,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "gamma")$estimate
#'     )
#'   ),
#'   .x = seq(1e-3, max(x), length.out = 100)
#' )
#'
#' @family Distributions
dist_gamma <- function(shape = NULL, rate = NULL) {
  GammaDistribution$new(shape = shape, rate = rate)
}

GammaDistribution <- distribution_class_simple(
  name = "Gamma",
  fun_name = "gamma",
  params = list(shape = I_POSITIVE_REALS, rate = I_POSITIVE_REALS),
  support = I_POSITIVE_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if ("rate" %in% names(vars)) {
      log_diff_rate <- params$shape / params$rate - x

      res$rate <- if (log) {
        log_diff_rate
      } else {
        log_diff_rate * dgamma(x, shape = params$shape, rate = params$rate)
      }
    }

    if ("shape" %in% names(vars)) {
      log_diff_shape <- log(x) + log(params$rate) - digamma(params$shape)

      res$shape <- if (log) {
        log_diff_shape
      } else {
        log_diff_shape * dgamma(x, shape = params$shape, rate = params$rate)
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if ("rate" %in% names(vars)) {
      diff_rate <- q / params$rate * dgamma(
        q, shape = params$shape, rate = params$rate
      )
      diff_rate[is.infinite(q)] <- 0.0

      res$rate <- if (log.p) {
        diff_rate / pgamma(q, shape = params$shape, rate = params$rate,
                           lower.tail = lower.tail)
      } else {
        diff_rate
      }

      if (!lower.tail) res$rate <- -res$rate
    }

    if ("shape" %in% names(vars)) {
      # Analytic gradient contains Meijer G-function => compute numeric
      # derivative
      res$shape <- gradient(
        func = function(shape) {
          sum(pgamma(
            q, shape = shape, rate = params$rate,
            lower.tail = lower.tail, log.p = log.p
          ))
        },
        x = params$shape
      )
    }

    res
  },
  tf_logdensity = function() function(x, args) {
    shape <- tf$broadcast_to(args[["shape"]], x$shape)
    rate <- tf$broadcast_to(args[["rate"]], x$shape)

    ok <- x > 0
    x_safe <- tf$where(ok, x, 1.0)

    tf$where(
      ok,
      shape * log(rate) - tf$math$lgamma(shape) + (shape - 1.0) * log(x_safe) - rate * x_safe,
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) {
    shape <- tf$broadcast_to(args[["shape"]], qmin$shape)
    rate <- tf$broadcast_to(args[["rate"]], qmax$shape)

    qmin0 <- qmin <= 0.0
    qmax0 <- qmax <= 0.0
    qmaxInf <- !tf$math$is_finite(qmax) & qmax > 0.0
    qmin_safe <- tf$maximum(K$zero, qmin)
    qmax_safe <- tf$maximum(K$zero, tf$where(qmaxInf, K$zero, qmax))
    qmax_nok <- tf$where(qmax0, K$neg_inf, K$zero)

    tf$where(
      qmin0,
      tf$where(
        qmax0 | qmaxInf,
        qmax_nok,
        log(tf$math$igamma(shape, qmax_safe * rate))
      ),
      tf$where(
        qmaxInf,
        log(1.0 - tf$math$igamma(shape, qmin_safe * rate)),
        log(tf$math$igamma(shape, qmax_safe * rate) - tf$math$igamma(shape, qmin_safe * rate))
      )
    )
  }
)

#' @export
fit_dist_start.GammaDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 2L)
  if ("shape" %in% ph && "rate" %in% ph) {
    res$shape <- mom[1L]^2.0 / mom[2L]
    res$rate <- mom[1L] / mom[2L]
  } else if ("shape" %in% ph) {
    res$shape <- dist$get_params()$rate * mom[1L]
  } else { # > "rate" %in% ph
    res$rate <- dist$get_params()$shape / mom[1L]
  }
  res
}

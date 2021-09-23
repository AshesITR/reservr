#' Weibull Distribution
#'
#' See [stats::Weibull]
#'
#' Both parameters can be overridden with
#' `with_params = list(shape = ..., scale = ...)`.
#'
#' @param shape Scalar shape parameter, or `NULL` as a placeholder.
#' @param scale Scalar scale parameter, or `NULL` as a placeholder.
#'
#' @return A `WeibullDistribution` object.
#' @export
#'
#' @examples
#' d_weibull <- dist_weibull(shape = 3, scale = 1)
#' x <- d_weibull$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_weibull,
#'   estimated = d_weibull,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "weibull")$estimate
#'     )
#'   ),
#'   .x = seq(0, 2, length.out = 100)
#' )
#'
#' @family Distributions
dist_weibull <- function(shape = NULL, scale = NULL) {
  WeibullDistribution$new(shape = shape, scale = scale)
}

WeibullDistribution <- distribution_class_simple(
  name = "Weibull",
  fun_name = "weibull",
  params = list(
    shape = I_POSITIVE_REALS,
    scale = I_POSITIVE_REALS
  ),
  support = I_POSITIVE_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if (length(vars)) {
      z <- x / params$scale
    }

    if ("shape" %in% names(vars)) {
      log_diff_shape <- 1.0 / params$shape + (1 - z^(params$shape)) * log(z)

      res$shape <- if (log) {
        log_diff_shape
      } else {
        log_diff_shape * dweibull(x, shape = params$shape, scale = params$scale)
      }
    }

    if ("scale" %in% names(vars)) {
      log_diff_scale <- params$shape / params$scale * (z^(params$shape) - 1.0)

      res$scale <- if (log) {
        log_diff_scale
      } else {
        log_diff_scale * dweibull(x, shape = params$shape, scale = params$scale)
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if (length(vars)) {
      z <- q / params$scale
      z_pow_alpha <- z^(params$shape)
    }

    if ("shape" %in% names(vars)) {
      z[is.infinite(q)] <- 0.0
      diff_shape <- log(z) * z_pow_alpha * exp(-z_pow_alpha)
      diff_shape[is.infinite(q)] <- 0.0

      res$shape <- if (log.p) {
        diff_shape / pweibull(q, shape = params$shape, scale = params$scale,
                              lower.tail = lower.tail)
      } else {
        diff_shape
      }

      if (!lower.tail) res$shape <- -res$shape
    }

    if ("scale" %in% names(vars)) {
      diff_scale <- params$shape / params$scale *
        z_pow_alpha *
        exp(-z_pow_alpha)
      diff_scale[is.infinite(q)] <- 0.0

      res$scale <- if (log.p) {
        diff_scale / pweibull(q, shape = params$shape, scale = params$scale,
                              lower.tail = lower.tail)
      } else {
        diff_scale
      }

      if (lower.tail) res$scale <- -res$scale
    }

    res
  },
  tf_logdensity = function() function(x, args) {
    shape <- tf$broadcast_to(args[["shape"]], x$shape)
    scale <- tf$broadcast_to(args[["scale"]], x$shape)

    ok <- x >= 0.0 & tf$math$is_finite(x)
    x_safe <- tf$where(ok, x, 0.0)

    tf$where(
      ok,
      log(shape) + (shape - 1.0) * log(x_safe) - shape * log(scale) - tf$math$pow(x_safe / scale, shape),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) {
    shape <- tf$broadcast_to(args[["shape"]], qmin$shape)
    scale <- tf$broadcast_to(args[["scale"]], qmin$shape)

    qmin0 <- qmin <= 0.0
    qmin_safe <- tf$math$maximum(0.0, qmin / scale)
    qmax0 <- qmax > 0.0
    qmax_ok <- tf$math$is_finite(qmax) & qmax > 0.0
    qmax_safe <- tf$where(qmax_ok, qmax / scale, qmin_safe + 1.0)
    qmax_nok <- tf$where(qmax0, K$neg_inf, K$zero)

    tf$where(
      qmin0,
      tf$where(
        qmax_ok,
        log(1.0 - exp(-tf$math$pow(qmax_safe, shape))),
        qmax_nok
      ),
      tf$where(
        qmax_ok,
        log(exp(-tf$math$pow(qmin_safe, shape)) - exp(-tf$math$pow(qmax_safe, shape))),
        -tf$math$pow(qmin_safe, shape)
      )
    )
  }
)

#' @export
fit_dist_start.WeibullDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  logmom <- weighted_moments(log(x), obs$w, n = 2L)
  if ("shape" %in% ph) {
    shape <- 1.2 / sqrt(logmom[2L])
    res$shape <- shape
  } else {
    shape <- dist$get_params()$shape
  }
  if ("scale" %in% ph) {
    res$scale <- exp(logmom[1L] + 0.572 / shape)
  }
  res
}

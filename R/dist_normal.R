#' Normal distribution
#'
#' See [stats::Normal].
#'
#' Both parameters can be overridden with
#' `with_params = list(mean = ..., sd = ...)`.
#'
#' @param mean Scalar mean parameter, or `NULL` as a placeholder.
#' @param sd Scalar standard deviation parameter, or `NULL` as a placeholder.
#'
#' @return A `NormalDistribution` object.
#' @export
#'
#' @examples
#' mu <- 0
#' sigma <- 1
#'
#' d_norm <- dist_normal(mean = mu, sd = sigma)
#' x <- d_norm$sample(20)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_norm,
#'   estimated = d_norm,
#'   with_params = list(
#'     estimated = list(mean = mean(x), sd = sd(x))
#'   ),
#'   .x = seq(-3, 3, length.out = 100)
#' )
#'
#' @family Distributions
dist_normal <- function(mean = NULL, sd = NULL) {
  NormalDistribution$new(mean = mean, sd = sd)
}

NormalDistribution <- distribution_class_simple(
  name = "normal",
  fun_name = "norm",
  params = list(mean = I_REALS, sd = I_POSITIVE_REALS),
  support = I_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if (length(vars)) {
      z <- (x - params$mean) / params$sd
    }

    if ("mean" %in% names(vars)) {
      res$mean <- if (log) {
        z / params$sd
      } else {
        z / params$sd * dnorm(x, mean = params$mean, sd = params$sd)
      }
    }

    if ("sd" %in% names(vars)) {
      res$sd <- if (log) {
        (z^2 - 1) / params$sd
      } else {
        (z^2 - 1) / params$sd * dnorm(x, mean = params$mean, sd = params$sd)
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if ("mean" %in% names(vars)) {
      res$mean <- if (log.p) {
        -dnorm(q, mean = params$mean, sd = params$sd) /
          pnorm(q, mean = params$mean, sd = params$sd, lower.tail = lower.tail)
      } else {
        -dnorm(q, mean = params$mean, sd = params$sd)
      }
      if (!lower.tail) res$mean <- -res$mean
    }

    if ("sd" %in% names(vars)) {
      z <- (q - params$mean) / params$sd
      res$sd <- if (log.p) {
        -z * dnorm(q, mean = params$mean, sd = params$sd) /
          pnorm(q, mean = params$mean, sd = params$sd, lower.tail = lower.tail)
      } else {
        ifelse(
          is.infinite(q),
          0.0,
          -z * dnorm(q, mean = params$mean, sd = params$sd)
        )
      }

      if (!lower.tail) res$sd <- -res$sd
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    mu <- tf$broadcast_to(tf$squeeze(args[["mean"]]), x$shape)
    sigma <- tf$broadcast_to(tf$squeeze(args[["sd"]]), x$shape)

    z <- (x - mu) / sigma

    -log(sigma) - K$log_sqrt_2pi - K$one_half * z * z
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    mu <- tf$broadcast_to(tf$squeeze(args[["mean"]]), qmin$shape)
    sigma <- tf$broadcast_to(tf$squeeze(args[["sd"]]), qmin$shape)

    qmin_finite <- tf$math$is_finite(qmin)
    qmin_safe <- tf$where(qmin_finite, qmin, K$zero)
    qmax_finite <- tf$math$is_finite(qmax)
    qmax_safe <- tf$where(qmax_finite, qmax, qmin_safe + K$one)

    zmin <- (qmin_safe - mu) / sigma / K$sqrt_2
    zmax <- (qmax_safe - mu) / sigma / K$sqrt_2

    tf$where(
      qmin_finite,
      tf$where(
        qmax_finite,
        log(tf$math$erf(zmax) - tf$math$erf(zmin)) - K$log_2,
        log(K$one - tf$math$erf(zmin)) - K$log_2
      ),
      tf$where(
        qmax_finite,
        log(K$one + tf$math$erf(zmax)) - K$log_2,
        K$zero
      )
    )
  }
)

#' @export
fit_dist_start.NormalDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 2L)
  if ("mean" %in% ph) {
    res$mean <- mom[1L]
  }
  if ("sd" %in% ph) {
    res$sd <- sqrt(mom[2L])
  }
  res
}

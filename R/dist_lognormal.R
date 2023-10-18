#' Log Normal distribution
#'
#' See [stats::Lognormal].
#'
#' Both parameters can be overridden with
#' `with_params = list(meanlog = ..., sdlog = ...)`.
#'
#' @param meanlog Scalar mean parameter on the log scale,
#' or `NULL` as a placeholder.
#' @param sdlog Scalar standard deviation parameter on the log scale,
#' or `NULL` as a placeholder.
#'
#' @return A `LognormalDistribution` object.
#' @export
#'
#' @examples
#' mu <- 0
#' sigma <- 1
#'
#' d_lnorm <- dist_lognormal(meanlog = mu, sdlog = sigma)
#' x <- d_lnorm$sample(20)
#' d_emp <- dist_empirical(x, positive = TRUE)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_lnorm,
#'   estimated = d_lnorm,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "lnorm")$estimate
#'     )
#'   ),
#'   .x = seq(1e-3, 5, length.out = 100)
#' )
#'
#' @family Distributions
dist_lognormal <- function(meanlog = NULL, sdlog = NULL) {
  LogNormalDistribution$new(meanlog = meanlog, sdlog = sdlog)
}

LogNormalDistribution <- distribution_class_simple(
  name = "LogNormal",
  fun_name = "lnorm",
  params = list(meanlog = I_REALS, sdlog = I_POSITIVE_REALS),
  support = I_POSITIVE_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if (length(vars)) {
      z <- (log(pmax(0, x)) - params$meanlog) / params$sdlog
    }

    if ("meanlog" %in% names(vars)) {
      res$meanlog <- if (log) {
        z / params$sdlog
      } else {
        z / params$sdlog * dlnorm(
          x, meanlog = params$meanlog, sdlog = params$sdlog
        )
      }
    }

    if ("sdlog" %in% names(vars)) {
      res$sdlog <- if (log) {
        (z^2 - 1) / params$sdlog
      } else {
        (z^2 - 1) / params$sdlog *
          dlnorm(x, meanlog = params$meanlog, sdlog = params$sdlog)
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if ("meanlog" %in% names(vars)) {
      diff_meanlog <- -q *
        dlnorm(q, meanlog = params$meanlog, sdlog = params$sdlog)
      diff_meanlog[is.infinite(q) | q <= 0.0] <- 0.0

      res$meanlog <- if (log.p) {
        diff_meanlog / plnorm(q, meanlog = params$meanlog, sdlog = params$sdlog,
                              lower.tail = lower.tail)
      } else {
        diff_meanlog
      }

      if (!lower.tail) res$meanlog <- -res$meanlog
    }

    if ("sdlog" %in% names(vars)) {
      z <- (log(pmax(0.0, q)) - params$meanlog) / params$sdlog
      diff_sdlog <- -z * q *
        dlnorm(q, meanlog = params$meanlog, sdlog = params$sdlog)
      diff_sdlog[is.infinite(q) | q <= 0.0] <- 0.0

      res$sdlog <- if (log.p) {
        diff_sdlog / plnorm(q, meanlog = params$meanlog, sdlog = params$sdlog,
                            lower.tail = lower.tail)
      } else {
        diff_sdlog
      }

      if (!lower.tail) res$sdlog <- -res$sdlog
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    meanlog <- args[["meanlog"]]
    sdlog <- args[["sdlog"]]

    ok <- x > K$zero
    x_safe <- tf$where(ok, x, K$one)

    z <- (log(x_safe) - meanlog) / sdlog

    tf$where(
      ok,
      -log(x_safe) - log(sdlog) - K$log_sqrt_2pi - K$one_half * z * z,
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    meanlog <- args[["meanlog"]]
    sdlog <- args[["sdlog"]]

    qmin_ok <- qmin > K$zero
    qmin_safe <- tf$where(qmin_ok, qmin, K$one)
    qmax_ok <- qmax > K$zero & tf$math$is_finite(qmax)
    qmax_safe <- tf$where(qmax_ok, qmax, qmin_safe + K$one)

    zmin <- (log(qmin_safe) - meanlog) / sdlog / K$sqrt_2
    zmax <- (log(qmax_safe) - meanlog) / sdlog / K$sqrt_2

    qmax_nok <- tf$where(
      qmax <= K$zero,
      K$neg_inf,
      K$zero
    )

    tf$where(
      qmin_ok,
      tf$where(
        qmax_ok,
        log(tf$math$erf(zmax) - tf$math$erf(zmin)) - K$log_2,
        log(K$one - tf$math$erf(zmin)) - K$log_2
      ),
      tf$where(
        qmax_ok,
        log(K$one + tf$math$erf(zmax)) - K$log_2,
        qmax_nok
      )
    )
  }
)

#' @export
fit_dist_start.LogNormalDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  logx <- log(x)
  logmom <- weighted_moments(logx, obs$w, n = 2L)
  if ("meanlog" %in% ph) {
    res$meanlog <- logmom[1L]
  }
  if ("sdlog" %in% ph) {
    res$sdlog <- sqrt(logmom[2L])
  }
  res
}

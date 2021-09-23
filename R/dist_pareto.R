#' Pareto Distribution
#'
#' See [actuar::Pareto]
#'
#' Both parameters can be overridden with
#' `with_params = list(shape = ..., scale = ...)`.
#'
#' @param shape Scalar shape parameter, or `NULL` as a placeholder.
#' @param scale Scalar scale parameter, or `NULL` as a placeholder.
#'
#' @return A `ParetoDistribution` object.
#' @export
#'
#' @examples
#' library(actuar) # so fitdistrplus finds it
#'
#' d_pareto <- dist_pareto(shape = 3, scale = 1)
#' x <- d_pareto$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_pareto,
#'   estimated = d_pareto,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "pareto")$estimate
#'     )
#'   ),
#'   .x = seq(0, 2, length.out = 100)
#' )
#'
#' @family Distributions
dist_pareto <- function(shape = NULL, scale = NULL) {
  ParetoDistribution$new(shape = shape, scale = scale)
}

ParetoDistribution <- distribution_class_simple(
  name = "Pareto",
  fun_name = "pareto",
  params = list(
    shape = I_POSITIVE_REALS,
    scale = I_POSITIVE_REALS
  ),
  support = I_POSITIVE_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if ("shape" %in% names(vars)) {
      res$shape <- if (log) {
        1.0 / params$shape + log(params$scale) - log(x + params$scale)
      } else {
        (1.0 / params$shape + log(params$scale) - log(x + params$scale)) *
          actuar::dpareto(x, shape = params$shape, scale = params$scale)
      }
    }

    if ("scale" %in% names(vars)) {
      res$scale <- if (log) {
        params$shape / params$scale - (params$shape + 1.0) / (x + params$scale)
      } else {
        (params$shape / params$scale -
          (params$shape + 1.0) / (x + params$scale)) *
          actuar::dpareto(x, shape = params$shape, scale = params$scale)
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    # Avoid NaNs
    q <- pmax(0.0, q)

    if ("shape" %in% names(vars)) {
      res$shape <- if (log.p) {
        if (lower.tail) {
          (log(params$scale) - log(q + params$scale)) *
            actuar::ppareto(q, shape = params$shape, scale = params$scale,
                            lower.tail = FALSE) /
            actuar::ppareto(q, shape = params$shape, scale = params$scale)
        } else {
          log(params$scale) - log(q + params$scale)
        }
      } else {
        (log(params$scale) - log(q + params$scale)) *
          actuar::ppareto(
            q, shape = params$shape, scale = params$scale, lower.tail = FALSE
          )
      }

      if (!log.p) res$shape[is.infinite(q)] <- 0.0

      if (lower.tail) res$shape <- -res$shape
    }

    if ("scale" %in% names(vars)) {
      res$scale <- if (log.p) {
        if (lower.tail) {
          params$shape / params$scale * q / (q + params$scale) *
            actuar::ppareto(
              q, shape = params$shape, scale = params$scale, lower.tail = FALSE
            ) /
            actuar::ppareto(q, shape = params$shape, scale = params$scale)
        } else {
          params$shape / params$scale * q / (q + params$scale)
        }
      } else {
        (params$shape / params$scale * q / (q + params$scale)) *
          actuar::ppareto(
            q, shape = params$shape, scale = params$scale, lower.tail = FALSE
          )
      }

      if (!log.p) res$scale[is.infinite(q)] <- 0.0

      if (lower.tail) res$scale <- -res$scale
    }

    res
  },
  tf_logdensity = function() function(x, args) {
    shape <- tf$broadcast_to(args[["shape"]], x$shape)
    scale <- tf$broadcast_to(args[["scale"]], x$shape)

    ok <- x >= K$zero & tf$math$is_finite(x)
    x_safe <- tf$where(ok, x, K$zero)

    tf$where(
      ok,
      log(shape) - (shape + K$one) * log(K$one + x_safe / scale) - log(scale),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) {
    shape <- tf$broadcast_to(args[["shape"]], qmin$shape)
    scale <- tf$broadcast_to(args[["scale"]], qmin$shape)

    qmin0 <- qmin <= K$zero
    qmin_safe <- tf$math$maximum(K$zero, qmin)
    qmax0 <- qmax > K$zero
    qmax_ok <- tf$math$is_finite(qmax) & qmax > K$zero
    qmax_safe <- tf$where(qmax_ok, qmax, qmin_safe + K$one)
    qmax_nok <- tf$where(qmax0, K$neg_inf, K$zero)

    qmin_sc <- K$one + qmin_safe / scale
    qmax_sc <- K$one + qmax_safe / scale

    tf$where(
      qmin0,
      tf$where(
        qmax_ok,
        log(K$one - tf$math$pow(qmax_sc, -shape)),
        qmax_nok
      ),
      tf$where(
        qmax_ok,
        log(tf$math$pow(qmin_sc, -shape) - tf$math$pow(qmax_sc, -shape)),
        -shape * log(qmin_sc)
      )
    )
  },
  envir = getNamespace("actuar")
)

#' @export
fit_dist_start.ParetoDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 2L, center = FALSE)
  if ("scale" %in% ph & "shape" %in% ph) {
    sc <- (mom[1L] * mom[2L]) / (mom[2L] - 2.0 * mom[1L]^2.0)

    if (sc < 0.0) { # illegal estimates produced by moment matching
      med <- weighted_median(x, obs$w)
      sc <- med
      for (iter in seq_len(5L)) {
        shpinv <- weighted.mean(log(x + sc), obs$w) - log(sc)
        sc <- med / (2.0^shpinv - 1.0)
      }
      res$scale <- sc
      res$shape <- 1.0 / shpinv
    } else {
      res$scale <- sc
      res$shape <- 1.0 + sc / mom[1L]
    }
  } else if ("shape" %in% ph) {
    sc <- dist$get_params()$scale
    res$shape <- 1 / weighted.mean(log(x + sc) - log(sc), obs$w)
  } else { # > "scale" %in% ph
    shp <- dist$get_params()$shape
    if (shp <= 1.0) {
      res$scale <- weighted_median(x, obs$w) / (2.0^(1.0 / shp) - 1.0)
    } else {
      res$scale <- mom[1L] * (shp - 1.0)
    }
  }
  res
}

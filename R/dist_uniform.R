#' Uniform distribution
#'
#' See [stats::Uniform]
#'
#' Both parameters can be overridden with
#' `with_params = list(min = ..., max = ...)`.
#'
#' @param min Lower limit, or `NULL` as a placeholder.
#' @param max Upper limit, or `NULL` as a placeholder.
#'
#' @return A `UniformDistribution` object.
#' @export
#'
#' @examples
#' d_unif <- dist_uniform(min = 0, max = 1)
#' x <- d_unif$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_unif,
#'   estimated = d_unif,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "unif")$estimate
#'     )
#'   ),
#'   .x = seq(0, 1, length.out = 100)
#' )
#'
#' @family Distributions
dist_uniform <- function(min = NULL, max = NULL) {
  UniformDistribution$new(min = min, max = max)
}

UniformDistribution <- distribution_class_simple(
  name = "Uniform",
  fun_name = "unif",
  params = list(min = I_REALS, max = I_REALS),
  support = function(params) {
    mapply(
      make_interval_union,
      lowest = params$min,
      highest = params$max,
      MoreArgs = list(include_lowest = 1.0, include_highest = 1.0, integer = 0.0),
      SIMPLIFY = FALSE
    )
  },
  is_in_support = function(x, params) {
    params$min <= x & x <= params$max
  },
  diff_density = function(x, vars, log, params) {
    res <- vars

    if (length(vars)) {
      ddens <- self$density(x, with_params = params)
      if (!log) {
        ddens <- ddens^2.0
      } else {
        ddens[ddens == 0.0] <- NaN
      }
    }

    if ("min" %in% names(vars)) {
      res$min <- ddens
    }

    if ("max" %in% names(vars)) {
      res$max <- -ddens
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    if (length(vars)) {
      dens <- self$density(q, with_params = params)
      # (q - m) / (M - m) or (M - q) / (M - m)
    }

    res <- vars

    if ("min" %in% names(vars)) {
      res$min <- if (lower.tail) {
        if (log.p) {
          (q - params$max) / (q - params$min) * dens
        } else {
          (q - params$max) * dens^2
        }
      } else {
        if (log.p) {
          dens[dens == 0.0] <- NaN
          dens
        } else {
          -(q - params$max) * dens^2 # nolint (lintr bug)
        }
      }
    }

    if ("max" %in% names(vars)) {
      res$max <- if (lower.tail) {
        if (log.p) {
          -dens
        } else {
          -(q - params$min) * dens^2 # nolint (lintr bug)
        }
      } else {
        if (log.p) {
          dens[dens == 0.0] <- NaN
          -(q - params$min) / (q - params$max) * dens # nolint (lintr bug)
        } else {
          (q - params$min) * dens^2
        }
      }
    }

    res
  },
  tf_logdensity = function() function(x, args) {
    bmin <- args[["min"]]
    bmax <- args[["max"]]

    tf$where(
      x >= tf$minimum(bmin, bmax) &
        x <= tf$maximum(bmin, bmax),
      -log(abs(bmax - bmin)),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) {
    amin <- args[["min"]]
    amax <- args[["max"]]

    bmin <- tf$minimum(amin, amax)
    bmax <- tf$maximum(amin, amax)

    supp_width <- log(bmax - bmin)
    qmin_clamp <- tf$minimum(bmax, tf$maximum(bmin, qmin))
    qmax_clamp <- tf$minimum(bmax, tf$maximum(bmin, qmax))

    log(qmax_clamp - qmin_clamp) - supp_width
  }
)

#' @export
fit_dist_start.UniformDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs)
  res <- dist$get_placeholders()
  ph <- names(res)
  if ("min" %in% ph) {
    res$min <- min(x)
  }
  if ("max" %in% ph) {
    res$max <- max(x)
  }
  res
}

#' @export
fit_dist.UniformDistribution <- function(dist, obs, start, ...) {
  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start)

  res <- list()
  if ("min" %in% names(start)) {
    res$min <- min(obs$xmin)
  }

  if ("max" %in% names(start)) {
    res$max <- max(obs$xmax)
  }

  list(params = res)
}

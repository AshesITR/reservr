#' Generalized Pareto Distribution
#'
#' See [evmix::gpd]
#'
#' All parameters can be overridden with
#' `with_params = list(u = ..., sigmau = ..., xi = ...)`.
#'
#' @param u Scalar location parameter, or `NULL` as a placeholder.
#' @param sigmau Scalar scale parameter, or `NULL` as a placeholder.
#' @param xi Scalar shape parameter, or `NULL` as a placeholder.
#'
#' @return A `GeneralizedParetoDistribution` object.
#' @export
#'
#' @examples
#' d_genpareto <- dist_genpareto(u = 0, sigmau = 1, xi = 1)
#' x <- d_genpareto$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' d_genpareto$export_functions("gpd") # so fitdistrplus finds it
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_genpareto,
#'   estimated = d_genpareto,
#'   with_params = list(
#'     estimated = fit(dist_genpareto(), x)$params
#'   ),
#'   .x = seq(0, 5, length.out = 100)
#' )
#'
#' @family Distributions
dist_genpareto <- function(u = NULL, sigmau = NULL, xi = NULL) {
  GeneralizedParetoDistribution$new(u = u, sigmau = sigmau, xi = xi)
}

#' @rdname dist_genpareto
#' @details `dist_genpareto1` is equivalent to `dist_genpareto` but enforces
#' bound constraints on `xi` to `[0, 1]`.
#' This ensures unboundedness and finite expected value unless `xi == 1.0`.
#' @export
dist_genpareto1 <- function(u = NULL, sigmau = NULL, xi = NULL) {
  out <- GeneralizedParetoDistribution$new(u = u, sigmau = sigmau, xi = xi)
  out$param_bounds$xi <- I_UNIT_INTERVAL
  out
}

#' @include gpd.R
GeneralizedParetoDistribution <- distribution_class_simple(
  name = "GeneralizedPareto",
  fun_name = "gpd",
  params = list(
    u = I_REALS,
    sigmau = I_POSITIVE_REALS,
    xi = I_REALS
  ),
  support = function(x, params) {
    xi_neg <- params$xi < 0.0
    res <- x >= params$u
    res[xi_neg] <- res[xi_neg] &
      (x <= params$u[xi_neg] - params$sigmau[xi_neg] / params$xi[xi_neg])
    res
  },
  diff_density = function(x, vars, log, params) {
    res <- vars

    if (length(vars)) {
      z <- (x - params$u) / params$sigmau
      zxip1 <- 1.0 + params$xi * z
      oob <- z < 0.0 | zxip1 < 0.0
      oob[is.na(oob)] <- FALSE

      if (!log) {
        dens <- dgpd(
          x = x,
          u = params$u,
          sigmau = params$sigmau,
          xi = params$xi
        )
      }
    }

    if ("u" %in% names(vars)) {
      res$u <- (1.0 + params$xi) / (params$sigmau * zxip1)
      if (!log) res$u <- res$u * dens
      res$u[oob] <- if (log) NaN else 0.0
    }

    if ("sigmau" %in% names(vars)) {
      xi0 <- params$xi %in% 0.0
      res$sigmau <- z * (1.0 + params$xi) /
        (params$sigmau * zxip1) -
        1.0 / params$sigmau
      res$sigmau[xi0] <- (z[xi0] - 1.0) / params$sigmau[xi0]
      if (!log) res$sigmau <- res$sigmau * dens
      res$sigmau[oob] <- if (log) NaN else 0.0
    }

    if ("xi" %in% names(vars)) {
      xi0 <- params$xi %in% 0.0
      res$xi <- numeric(length(z))
      res$xi[!xi0 & !oob] <- log(zxip1[!xi0 & !oob]) /
        params$xi[!xi0 & !oob]^2.0 -
        ((1.0 + params$xi[!xi0 & !oob]) * z[!xi0 & !oob]) /
          (params$xi[!xi0 & !oob] * zxip1[!xi0 & !oob])
      res$xi[xi0] <- z[xi0]^2.0 / 2.0 - z[xi0]
      if (!log) res$xi <- res$xi * dens
      res$xi[oob] <- if (log) NaN else 0.0
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if (length(vars)) {
      z <- (q - params$u) / params$sigmau
      zxip1 <- 1.0 + params$xi * z

      oob <- z < 0.0 | zxip1 < 0.0
      oob[is.na(z)] <- TRUE
      xi0 <- params$xi == 0.0

      dens <- dgpd(
        x = q,
        u = params$u,
        sigmau = params$sigmau,
        xi = params$xi
      )

      prob <- pgpd(
        q = q,
        u = params$u,
        sigmau = params$sigmau,
        xi = params$xi,
        lower.tail = lower.tail
      )
    }

    if ("u" %in% names(vars)) {
      res$u <- if (lower.tail) -dens else dens
      if (log.p) res$u <- res$u / prob
    }

    if ("sigmau" %in% names(vars)) {
      res$sigmau <- z * dens
      res$sigmau[is.infinite(z)] <- 0.0
      if (lower.tail) res$sigmau <- -res$sigmau
      if (log.p) res$sigmau <- res$sigmau / prob
    }

    if ("xi" %in% names(vars)) {
      res$xi <- numeric(length(z))
      res$xi[!xi0 & !oob] <- log(zxip1[!xi0 & !oob]) /
        params$xi[!xi0 & !oob]^2.0 -
        z[!xi0 & !oob] / (zxip1[!xi0 & !oob] * params$xi[!xi0 & !oob])
      res$xi[xi0] <- 0.5 * z[xi0]^2.0
      res$xi[oob] <- 0.0
      res$xi[is.na(z)] <- z[is.na(z)]
      res$xi[is.infinite(z)] <- 0.0

      if (lower.tail) {
        prob_ut <- pgpd(
          q = q,
          u = params$u,
          sigmau = params$sigmau,
          xi = params$xi,
          lower.tail = FALSE
        )

        if (log.p) {
          res$xi <- -res$xi * (prob_ut / prob)
        } else {
          res$xi <- -res$xi * prob_ut
        }
      } else {
        if (!log.p) {
          res$xi <- res$xi * prob
        } else {
          res$xi[oob & zxip1 < 0.0] <- NaN
        }
      }
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    u <- tf$squeeze(args[["u"]])
    sigmau <- tf$squeeze(args[["sigmau"]])
    xi <- tf$squeeze(args[["xi"]])

    z <- (x - u) / sigmau
    zxip1 <- K$one + xi * z

    z_safe <- tf$where(z >= K$zero, z, K$zero)
    zxip1_safe <- tf$where(zxip1 > K$zero, zxip1, K$one)

    xi_nonzero <- tf$where(xi == K$zero, K$one, xi)

    tf$where(
      z >= K$zero & zxip1 > K$zero,
      tf$where(
        xi == K$zero,
        -z_safe +
          # Add artificial derivative w.r.t. xi
          xi * (K$one_half * tf$math$square(z_safe) - z_safe),
        (K$neg_one / xi_nonzero - K$one) * log(zxip1_safe)
      ) -
        log(sigmau),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    tf <- tensorflow::tf
    u <- tf$squeeze(args[["u"]])
    sigmau <- tf$squeeze(args[["sigmau"]])
    xi <- tf$squeeze(args[["xi"]])

    qmin_finite <- tf$maximum(u, tf$where(qmin < K$inf, qmin, K$one))
    qmax_finite <- tf$maximum(u, tf$where(qmax < K$inf, qmax, K$one))

    zmin <- tf$math$maximum(
      K$zero,
      tf$where(qmin < K$inf, (qmin_finite - u) / sigmau, K$one)
    )

    zmax <- tf$math$maximum(
      K$zero,
      tf$where(qmax < K$inf, (qmax_finite - u) / sigmau, K$one)
    )

    zmin_xip1 <- tf$math$maximum(K$zero, K$one + xi * zmin)
    zmax_xip1 <- tf$math$maximum(K$zero, K$one + xi * zmax)

    zmin_safe <- tf$where(zmax == zmin, K$zero, zmin)
    zmax_safe <- tf$where(zmax == zmin, K$one, zmax)
    zmin_xip1_safe <- tf$where(zmax_xip1 == zmin_xip1, K$one, zmin_xip1)
    zmax_xip1_safe <- tf$where(zmax_xip1 == zmin_xip1, K$two, zmax_xip1)

    # Safe version of xi (guaranteed nonzero) for use within actual tf$where
    # The 1.0 branch will never actually be selected, but prevents a division
    # by zero in the off-branch
    # cf. https://github.com/tensorflow/tensorflow/issues/38349
    xi_nonzero <- tf$where(
      tf$math$equal(xi, K$zero),
      K$one,
      xi
    )

    log(
      tf$where(
        xi == K$zero,
        tf$where(
          zmin == zmax,
          K$zero,
          exp(-zmin_safe) -
            tf$where(qmax < K$inf, exp(-zmax_safe), K$zero)
        ),
        tf$where(
          zmin_xip1 == zmax_xip1,
          K$zero,
          zmin_xip1_safe^(K$neg_one / xi_nonzero) -
            tf$where(
              qmax < K$inf,
              zmax_xip1_safe^(K$neg_one / xi_nonzero),
              K$zero
            )
        )
      )
    )
  }
)

#' @export
fit_dist_start.GeneralizedParetoDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)

  res <- dist$get_placeholders()
  ph <- names(res)

  if ("u" %in% ph) {
    u <- min(obs$xmin[is.finite(obs$xmin)])
    res$u <- u
  } else {
    u <- dist$get_params()$u
  }

  x <- .get_init_x(obs, .min = u)

  if ("xi" %in% ph) {
    hillord <- ceiling(sum(x > u) / 2.0)
    if (hillord > 0L) {
      xs <- sort(x[x > u], decreasing = TRUE) - u
      xi <- mean(log(xs[seq_len(hillord - 1L)])) - log(xs[hillord])
    } else {
      # all obs are <= u
      xi <- 0.0
    }

    # Handle dist_genpareto1 initialisation
    xi_bounds <- dist$get_param_bounds()$xi
    if (!xi_bounds$contains(xi)) {
      xi <- pmin(pmax(xi_bounds$range[1L], xi), xi_bounds$range[2L])
    }

    res$xi <- xi
  } else {
    xi <- dist$get_params()$xi
  }

  if ("sigmau" %in% ph) {
    res$sigmau <- if (xi < 1.0) {
      (mean(x) - u) * (1.0 - xi) # moment matching
    } else {
      (median(x) - u) * xi / (2.0^xi - 1.0) # median matching
    }
  }

  res
}

#' @export
fit_dist.GeneralizedParetoDistribution <- function(dist, obs, start, ...) {
  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start)

  if ("u" %in% names(start) &&
      identical(dist$get_param_bounds()[["u"]], I_REALS)) {
    # MLE of u is unstable for gradient-based ML because the right derivative at
    # the optimum (= min(obs$x[obs$w > 0.0])) doesn't exist. To fix this, we add
    # artifical box constraints on u. If the initial value was manually
    # specified and exceeds this constraint, it will instead be considered
    # fixed. This can happen if the GPD is a component of some mixture and the
    # other w > 0.0 observations are guaranteed to stem from a different
    # component.

    u_max <- min(obs$xmax[obs$w > 0.0])
    if (start$u > u_max) {
      # Temporarily fix u
      dist$default_params$u <- start$u
      on.exit(dist$default_params$u <- NULL, add = TRUE)
      res <- fit_dist(dist = dist, obs = obs, start = start, ...)
      # TODO transform history if trace is TRUE
      res$params$u <- start$u
      res
    } else {
      # Temporarily bound u above
      dist$param_bounds$u <- interval(
        range = c(-Inf, u_max),
        include_highest = TRUE
      )
      on.exit(dist$param_bounds$u <- I_REALS, add = TRUE)

      fit_dist(dist = dist, obs = obs, start = start, ...)
    }
  } else {
    NextMethod("fit_dist")
  }
}

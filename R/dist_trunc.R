#' Truncated distribution
#'
#' @param dist An underlying distribution, or `NULL` as a placeholder.
#' @param min Minimum value to truncate at (exclusive),
#' or `NULL` as a placeholder.
#' @param max Maxmimum value to truncate at (inclusive),
#' or `NULL` as a placeholder.
#' @param offset Offset to be added to each observation after truncation,
#' or `NULL` as a placeholder. Truncation of `dist` will occur to (min, max].
#' The offset is then added deterministically.
#' @param max_retry Maximum number of resample attempts when trying to sample
#' with rejection.
#'
#' @return A `TruncatedDistribution` object.
#' @export
#'
#' @examples
#' d_norm <- dist_normal(mean = 0, sd = 1)
#' d_tnorm <- dist_trunc(dist = d_norm, min = -2, max = 2, offset = 1)
#' plot_distributions(d_norm, d_tnorm, .x = seq(-2, 3, length.out = 100))
#'
#' @family Distributions
dist_trunc <- function(dist = NULL, min = NULL, max = NULL,
                       offset = 0, max_retry = 100) {
  TruncatedDistribution$new(dist = dist, min = min, max = max,
                            offset = offset, max_retry = max_retry)
}

TruncatedDistribution <- distribution_class(
  name = "Truncated",
  params = list(
    dist = list(),
    min = I_REALS,
    max = I_REALS,
    offset = I_REALS,
    max_retry = I_NATURALS
  ),
  sample = function(n, params) {
    if (all(params$dist$dist$has_capability(c("quantile", "probability")))) {
      p <- runif(n, min = 0.0, max = 1.0)

      wp <- params
      wp$dist <- wp$dist$params

      self$quantile(
        p = p,
        with_params = wp
      )
    } else {
      # Rejection sampling
      n_retry <- 0
      x <- params$dist$dist$sample(
        n = n,
        with_params = params$dist$params
      )

      min <- params$min
      max <- params$max

      while (all(n_retry < params$max_retry[min > x | x > max]) && any(min > x | x > max)) {
        not_ok <- min > x | x > max
        n_nok <- sum(not_ok)
        x[not_ok] <- params$dist$dist$sample(
          n = n_nok,
          with_params = pick_params_at(params$dist$params, not_ok)
        )
        n_retry <- n_retry + 1
      }

      if (any(min > x | x > max)) {
        stop("Rejection sampling exceeded max_retry.") # nocov
      }

      x + params$offset
    }
  },
  density = function(x, log = FALSE, params) {
    params$dist$dist$require_capability(c("density", "probability"),
                                        fun_name = "dist_trunc$density()")

    x_trans <- x - params$offset
    x_okay <- !is.na(x_trans) &
      params$min <= x_trans &
      x_trans <= params$max

    res <- numeric(length(x))

    res[is.na(x_trans)] <- x_trans[is.na(x_trans)] # Handle NA / NaN

    cdf_min <- params$dist$dist$probability(
      q = params$min,
      with_params = params$dist$params
    )

    if (!params$dist$dist$is_continuous()) {
      disc <- params$dist$dist$is_discrete_at(params$min, with_params = params$dist$params)
      cdf_min[disc] <- cdf_min[disc] - params$dist$dist$density(
        x = params$min[disc],
        with_params = params$dist$params
      )
    }

    cdf_max <- params$dist$dist$probability(
      q = params$max,
      with_params = params$dist$params
    )

    if (log) {
      res[x_okay] <- params$dist$dist$density(
        x = x_trans[x_okay],
        log = TRUE,
        with_params = pick_params_at(params$dist$params, x_okay)
      ) - log(cdf_max[x_okay] - cdf_min[x_okay])

      res[x_trans < params$min] <- -Inf
      res[x_trans > params$max] <- -Inf
    } else {
      res[x_okay] <- params$dist$dist$density(
        x = x_trans[x_okay],
        with_params = pick_params_at(params$dist$params, x_okay)
      ) / (cdf_max[x_okay] - cdf_min[x_okay])

      res[x_trans < params$min] <- 0
      res[x_trans > params$max] <- 0
    }

    res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    params$dist$dist$require_capability("probability",
                                        fun_name = "dist_trunc$probability()")
    cdf_min <- params$dist$dist$probability(
      q = params$min,
      with_params = params$dist$params
    )

    if (!params$dist$dist$is_continuous()) {
      disc <- params$dist$dist$is_discrete_at(params$min, with_params = params$dist$params)
      cdf_min[disc] <- cdf_min[disc] - params$dist$dist$density(
        x = params$min[disc],
        with_params = params$dist$params
      )
    }

    cdf_max <- params$dist$dist$probability(
      q = params$max,
      with_params = params$dist$params
    )

    q_trans <- q - params$offset
    q_okay <- !is.na(q_trans) &
      params$min <= q_trans &
      q_trans <= params$max

    res <- numeric(length(q))

    res[is.na(q_trans)] <- q_trans[is.na(q_trans)]

    dist_prob_okay <- params$dist$dist$probability(
      q = q_trans[q_okay],
      with_params = pick_params_at(params$dist$params, q_okay)
    )

    if (log.p) {
      if (lower.tail) {
        res[q_okay] <- log(dist_prob_okay - cdf_min[q_okay]) -
          log(cdf_max[q_okay] - cdf_min[q_okay])

        res[q_trans < params$min] <- -Inf
        res[q_trans > params$max] <- 0.0
      } else {
        res[q_okay] <- log(cdf_max[q_okay] - dist_prob_okay) -
          log(cdf_max[q_okay] - cdf_min[q_okay])

        res[q_trans < params$min] <- 0.0
        res[q_trans > params$max] <- -Inf
      }
    } else {
      if (lower.tail) {
        res[q_okay] <- (dist_prob_okay - cdf_min[q_okay]) /
          (cdf_max[q_okay] - cdf_min[q_okay])

        res[q_trans < params$min] <- 0.0
        res[q_trans > params$max] <- 1.0
      } else {
        res[q_okay] <- (cdf_max[q_okay] - dist_prob_okay) /
          (cdf_max[q_okay] - cdf_min[q_okay])

        res[q_trans < params$min] <- 1.0
        res[q_trans > params$max] <- 0.0
      }
    }

    res
  },
  quantile = function(p, lower.tail = TRUE, log.p = FALSE, params) {
    params$dist$dist$require_capability(c("quantile", "probability"),
                                        fun_name = "dist_trunc$quantile()")

    cdf_min <- params$dist$dist$probability(
      q = params$min,
      with_params = params$dist$params
    )

    if (!params$dist$dist$is_continuous()) {
      disc <- params$dist$dist$is_discrete_at(params$min, with_params = params$dist$params)
      cdf_min[disc] <- cdf_min[disc] - params$dist$dist$density(
        x = params$min[disc],
        with_params = params$dist$params
      )
    }

    cdf_max <- params$dist$dist$probability(
      q = params$max,
      with_params = params$dist$params
    )

    if (log.p) p <- exp(p)
    if (lower.tail) {
      p_dist <- p * (cdf_max - cdf_min) + cdf_min
    } else {
      p_dist <- cdf_max - p * (cdf_max - cdf_min)
    }

    params$offset + params$dist$dist$quantile(
      p = p_dist,
      with_params = params$dist$params
    )
  },
  support = function(x, params) {
    x_trans <- x - params$offset
    x_okay <- !is.na(x_trans) &
      params$min <= x_trans &
      x_trans <= params$max
    res <- logical(length(x))
    res[!x_okay] <- FALSE
    res[x_okay] <- params$dist$dist$is_in_support(
      x = x_trans[x_okay],
      with_params = params$dist$params
    )
    res
  },
  get_components = function() {
    # Interface chosen identical to dist_mixture()$get_components()
    list(private$.default_params$dist)
  },
  has_capability = function(caps) {
    self$get_components()[[1L]]$has_capability(caps) &
      caps %in% private$.caps
  },
  tf_logdensity = function() {
    dist_logdens <- self$get_components()[[1L]]$tf_logdensity()
    dist_logprob <- self$get_components()[[1L]]$tf_logprobability()

    function(x, args) {
      min <- tf$broadcast_to(args[["min"]], x$shape)
      max <- tf$broadcast_to(args[["max"]], x$shape)
      offset <- tf$broadcast_to(args[["offset"]], x$shape)
      dist_args <- args[["dist"]]

      x_trans <- x - offset
      ok <- x_trans >= min & x_trans <= max
      xt_safe <- tf$where(ok, x_trans, min)

      tf$where(ok, dist_logdens(xt_safe, dist_args) - dist_logprob(min, max, dist_args), K$neg_inf)
    }
  },
  tf_logprobability = function() {
    dist_logprob <- self$get_components()[[1L]]$tf_logprobability()

    function(qmin, qmax, args) {
      min <- tf$broadcast_to(args[["min"]], qmin$shape)
      max <- tf$broadcast_to(args[["max"]], qmin$shape)
      offset <- tf$broadcast_to(args[["offset"]], qmin$shape)
      dist_args <- args[["dist"]]

      xt_min <- qmin - offset
      xt_max <- qmax - offset
      xt_min_safe <- tf$minimum(max, tf$maximum(min, xt_min))
      xt_max_safe <- tf$minimum(max, tf$maximum(min, xt_max))

      zero <- xt_max < min | xt_min > max

      tf$where(zero, K$neg_inf, dist_logprob(xt_min_safe, xt_max_safe, dist_args) - dist_logprob(min, max, dist_args))
    }
  },
  tf_is_discrete_at = function() {
    if (self$is_continuous()) return(super$tf_is_discrete_at())

    dist_disc <- self$get_components()[[1L]]$tf_is_discrete_at()

    function(x, args) {
      min <- tf$broadcast_to(args[["min"]], x$shape)
      max <- tf$broadcast_to(args[["max"]], x$shape)
      offset <- tf$broadcast_to(args[["offset"]], x$shape)
      dist_args <- args[["dist"]]

      x_trans <- x - offset
      ok <- x_trans >= min & x_trans <= max
      xt_safe <- tf$where(ok, x_trans, min)

      tf$where(ok, dist_disc(xt_safe, dist_args), FALSE)
    }
  },
  compile_sample = function() {
    ph <- names(self$get_placeholders())
    ph_min <- "min" %in% ph
    ph_max <- "max" %in% ph
    ph_offset <- "offset" %in% ph
    ph_max_retry <- "max_retry" %in% ph

    dist <- self$get_components()[[1L]]
    dist_quantile <- dist$compile_quantile()
    dist_prob <- dist$compile_probability()
    dist_prob_int <- dist$compile_probability_interval()

    n_params_dist <- attr(dist_quantile, "n_params")
    n_params <- as.integer(ph_min) + as.integer(ph_max) + as.integer(ph_offset) + as.integer(ph_max_retry) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    min_expr <- if (ph_min) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$min
    }

    max_expr <- if (ph_max) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + 1L)])
    } else {
      self$default_params$max
    }

    offset_expr <- if (ph_offset) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + as.integer(ph_max) + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(n, param_matrix) {
        p_upper <- dist_prob(.(max_expr), .(dist_param_expr))
        p_width <- dist_prob_int(.(min_expr), .(max_expr), .(dist_param_expr))
        q <- runif(n = n, min = p_upper - p_width, max = p_upper)
        .(offset_expr) + dist_quantile(q, .(dist_param_expr))
      })),
      n_params
    )
  },
  compile_density = function() {
    ph <- names(self$get_placeholders())
    ph_min <- "min" %in% ph
    ph_max <- "max" %in% ph
    ph_offset <- "offset" %in% ph
    ph_max_retry <- "max_retry" %in% ph

    dist <- self$get_components()[[1L]]
    dist_dens <- dist$compile_density()
    dist_prob_int <- dist$compile_probability_interval()

    n_params_dist <- attr(dist_dens, "n_params")
    n_params <- as.integer(ph_min) + as.integer(ph_max) + as.integer(ph_offset) + as.integer(ph_max_retry) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    min_expr <- if (ph_min) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$min
    }

    max_expr <- if (ph_max) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + 1L)])
    } else {
      self$default_params$max
    }

    offset_expr <- if (ph_offset) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + as.integer(ph_max) + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(x, param_matrix, log = FALSE) {
        xtrans <- x - .(offset_expr)
        xdens <- dist_dens(xtrans, .(dist_param_expr), log = log)
        ptrunc <- dist_prob_int(.(min_expr), .(max_expr), .(dist_param_expr), log.p = log)

        if (log) {
          xdens[xtrans < .(min_expr)] <- -Inf
          xdens[xtrans > .(max_expr)] <- -Inf
          xdens - ptrunc
        } else {
          xdens[xtrans < .(min_expr)] <- 0.0
          xdens[xtrans > .(max_expr)] <- 0.0
          xdens / ptrunc
        }
      })),
      n_params
    )
  },
  compile_probability = function() {
    ph <- names(self$get_placeholders())
    ph_min <- "min" %in% ph
    ph_max <- "max" %in% ph
    ph_offset <- "offset" %in% ph
    ph_max_retry <- "max_retry" %in% ph

    dist <- self$get_components()[[1L]]
    dist_prob_int <- dist$compile_probability_interval()

    n_params_dist <- attr(dist_prob_int, "n_params")
    n_params <- as.integer(ph_min) + as.integer(ph_max) + as.integer(ph_offset) + as.integer(ph_max_retry) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    min_expr <- if (ph_min) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$min
    }

    max_expr <- if (ph_max) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + 1L)])
    } else {
      self$default_params$max
    }

    offset_expr <- if (ph_offset) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + as.integer(ph_max) + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(q, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        qtrans <- q - .(offset_expr)
        if (lower.tail) {
          qprob <- dist_prob_int(.(min_expr), qtrans, .(dist_param_expr), log.p = log.p)
        } else {
          qprob <- dist_prob_int(qtrans, .(max_expr), .(dist_param_expr), log.p = log.p)
        }
        ptrunc <- dist_prob_int(.(min_expr), .(max_expr), .(dist_param_expr), log.p = log.p)

        if (log.p) {
          qprob[qtrans < .(min_expr)] <- if (lower.tail) -Inf else 0.0
          qprob[qtrans > .(max_expr)] <- if (lower.tail) 0.0 else -Inf
          qprob - ptrunc
        } else {
          qprob[qtrans < .(min_expr)] <- if (lower.tail) 0.0 else 1.0
          qprob[qtrans > .(max_expr)] <- if (lower.tail) 1.0 else 0.0
          qprob / ptrunc
        }
      })),
      n_params
    )
  },
  compile_probability_interval = function() {
    ph <- names(self$get_placeholders())
    ph_min <- "min" %in% ph
    ph_max <- "max" %in% ph
    ph_offset <- "offset" %in% ph
    ph_max_retry <- "max_retry" %in% ph

    dist <- self$get_components()[[1L]]
    dist_prob_int <- dist$compile_probability_interval()

    n_params_dist <- attr(dist_dens, "n_params")
    n_params <- as.integer(ph_min) + as.integer(ph_max) + as.integer(ph_offset) + as.integer(ph_max_retry) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    min_expr <- if (ph_min) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$min
    }

    max_expr <- if (ph_max) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + 1L)])
    } else {
      self$default_params$max
    }

    offset_expr <- if (ph_offset) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + as.integer(ph_max) + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
        qtrans_min <- pmin(pmax(qmin - .(offset_expr), .(min_expr)), .(max_expr))
        qtrans_max <- pmin(pmax(qmax - .(offset_expr), .(min_expr)), .(max_expr))

        qprob <- dist_prob_int(qtrans_min, qtrans_max, .(dist_param_expr), log.p = log.p)
        ptrunc <- dist_prob_int(.(min_expr), .(max_expr), .(dist_param_expr), log.p = log.p)

        if (log.p) {
          qprob - ptrunc
        } else {
          qprob / ptrunc
        }
      })),
      n_params
    )
  },
  compile_quantile = function() {
    ph <- names(self$get_placeholders())
    ph_min <- "min" %in% ph
    ph_max <- "max" %in% ph
    ph_offset <- "offset" %in% ph
    ph_max_retry <- "max_retry" %in% ph

    dist <- self$get_components()[[1L]]
    dist_quantile <- dist$compile_quantile()
    dist_prob <- dist$compile_probability()
    dist_prob_int <- dist$compile_probability_interval()

    n_params_dist <- attr(dist_quantile, "n_params")
    n_params <- as.integer(ph_min) + as.integer(ph_max) + as.integer(ph_offset) + as.integer(ph_max_retry) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    min_expr <- if (ph_min) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$min
    }

    max_expr <- if (ph_max) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + 1L)])
    } else {
      self$default_params$max
    }

    offset_expr <- if (ph_offset) {
      quote(param_matrix[, .(n_params_dist + as.integer(ph_min) + as.integer(ph_max) + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(p, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        if (log.p) p <- exp(p)
        p_upper <- dist_prob(.(max_expr), .(dist_param_expr))
        p_width <- dist_prob_int(.(min_expr), .(max_expr), .(dist_param_expr))

        q <- if (lower.tail) {
          p_upper - p_width * (1.0 - p)
        } else {
          p_upper - p_width * p
        }

        .(offset_expr) + dist_quantile(q, .(dist_param_expr))
      })),
      n_params
    )
  }
)

#' @export
fit_dist_start.TruncatedDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  res <- dist$get_placeholders()
  ph <- names(res)
  ph_dist <- length(res$dist) > 0L
  rng <- c(min(obs$xmin), max(obs$xmax))

  .assert_set(ph, c("offset", "max_retry"), "Truncated")

  if ("min" %in% ph) {
    res$min <- rng[1L]
  } else {
    rng[1L] <- dist$get_params()$min
  }
  if ("max" %in% ph) {
    res$max <- rng[2L]
  } else {
    rng[2L] <- dist$get_params()$max
  }
  if (ph_dist) {
    obs_tr <- truncate_obs(obs, tmin_new = rng[1L], tmax_new = rng[2L])

    # shift all data by offset (weights are unmodified)
    obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")] <- lapply(
      obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")],
      "-", dist$get_params()$offset
    )

    res$dist <- fit_dist_start(
      dist = dist$get_components()[[1L]],
      obs = obs_tr,
      ...
    )
  } else {
    res$dist <- NULL
  }
  res
}

#' Tranlsated distribution
#'
#' @param dist An underlying distribution, or `NULL` as a placeholder.
#' @param offset Offset to be added to each observation, or `NULL` as a placeholder.
#' @param multiplier Factor to multiply each observation by, or `NULL` as a placeholder.
#'
#' @return A `TranslatedDistribution` object.
#' @export
#'
#' @examples
#' d_norm <- dist_normal(mean = 0, sd = 1)
#' d_tnorm <- dist_translate(dist = d_norm, offset = 1)
#' plot_distributions(d_norm, d_tnorm, .x = seq(-2, 3, length.out = 100))
#'
#' @family Distributions
dist_translate <- function(dist = NULL, offset = NULL, multiplier = 1.0) {
  TranslatedDistribution$new(dist = dist, offset = offset, multiplier = multiplier)
}

TranslatedDistribution <- distribution_class(
  name = "Translated",
  params = list(
    dist = list(),
    offset = I_REALS,
    multiplier = I_REALS
  ),
  sample = function(n, params) {
    x <- params$dist$dist$sample(
      n = n,
      with_params = params$dist$params
    )
    x * params$multiplier + params$offset
  },
  density = function(x, log = FALSE, params) {
    params$dist$dist$require_capability(
      "density",
      fun_name = "dist_translate$density()"
    )
    x_trans <- (x - params$offset) / params$multiplier
    dens <- params$dist$dist$density(x_trans, log = log, with_params = params$dist$params)
    if (log) dens - log(params$multiplier) else dens / params$multiplier
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    params$dist$dist$require_capability(
      "probability",
      fun_name = "dist_translate$probability()"
    )
    q_trans <- (q - params$offset) / params$multiplier
    params$dist$dist$probability(
      q_trans, lower.tail = lower.tail, log.p = log.p,
      with_params = params$dist$params
    )
  },
  quantile = function(p, lower.tail = TRUE, log.p = FALSE, params) {
    params$offset + params$multiplier * params$dist$dist$quantile(
      p = p,
      lower.tail = lower.tail,
      log.p = log.p,
      with_params = params$dist$params
    )
  },
  support = function(x, params) {
    x_trans <- (x - params$offset) / params$multiplier
    params$dist$dist$is_in_support(
      x = x_trans,
      with_params = params$dist$params
    )
  },
  get_components = function() {
    # Interface chosen identical to dist_mixture()$get_components()
    list(private$.default_params$dist)
  },
  diff_density = function(x, vars, log, params) {
    res <- vars

    x_trans <- (x - params$offset) / params$multiplier

    if (length(vars$dist)) {
      res$dist <- params$dist$dist$diff_density(
        x_trans, log = log, with_params = params$dist$params
      )
    } else {
      res$dist <- NULL
    }

    if ("offset" %in% names(vars)) {
      stop("dist_translate$diff_density()$offset is not implemented yet.")
    }

    if ("multiplier" %in% names(vars)) {
      stop("dist_translate$diff_density()$multiplier is not implemented yet.")
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    q_trans <- (q - params$offset) / params$multiplier

    if (length(vars$dist)) {
      res$dist <- params$dist$dist$diff_probability(
        q_trans,
        lower.tail = lower.tail,
        log.p = log.p,
        with_params = params$dist$params
      )
    } else {
      res$dist <- NULL
    }

    if ("offset" %in% names(vars)) {
      res$offset <- params$dist$dist$density(
        q_trans,
        with_params = params$dist$params
      ) / params$multiplier

      if (log.p) {
        res$offset <- res$offset / params$dist$dist$probability(
          q_trans,
          lower.tail = lower.tail,
          with_params = params$dist$params
        )
      }

      if (lower.tail) {
        res$offset <- -res$offset
      }
    }

    if ("multiplier" %in% names(vars)) {
      res$multiplier <- params$dist$dist$density(
        q_trans,
        with_params = params$dist$params
      ) * q_trans / params$multiplier

      if (log.p) {
        res$multiplier <- res$multiplier / params$dist$dist$probability(
          q_trans,
          lower.tail = lower.tail,
          with_params = params$dist$params
        )
      }

      if (lower.tail) {
        res$multiplier <- -res$multiplier
      }
    }

    res
  },
  has_capability = function(caps) {
    self$get_components()[[1L]]$has_capability(caps)
  },
  tf_logdensity = function() {
    logd_comp <- self$get_components()[[1L]]$tf_logdensity()

    function(x, args) {
      tf <- tensorflow::tf
      multiplier <- tf$squeeze(args[["multiplier"]])
      invert <- multiplier < 0
      x_trans <- tf$where(tf$math$is_finite(x), x * multiplier - tf$squeeze(args[["offset"]]), tf$where(invert, -x, x))
      logd_comp(x_trans, args[["dist"]]) - log(abs(multiplier))
    }
  },
  tf_logprobability = function() {
    logp_comp <- self$get_components()[[1L]]$tf_logprobability()

    function(qmin, qmax, args) {
      tf <- tensorflow::tf
      multiplier <- tf$squeeze(args[["multiplier"]])
      offset <- tf$squeeze(args[["offset"]])
      invert <- multiplier < 0

      qmin_finite <- tf$where(tf$math$is_finite(qmin), qmin, K$zero)
      qmax_finite <- tf$where(tf$math$is_finite(qmax), qmax, K$zero)

      qmin_trans <- tf$where(
        tf$math$is_finite(qmin),
        (qmin_finite - offset) / multiplier,
        tf$where(invert, -qmin, qmin)
      )

      qmax_trans <- tf$where(
        tf$math$is_finite(qmax),
        (qmax_finite - offset) / multiplier,
        tf$where(invert, -qmax, qmax)
      )

      qmin_trans2 <- tf$where(invert, qmax_trans, qmin_trans)
      qmax_trans2 <- tf$where(invert, qmin_trans, qmax_trans)

      logp_comp(qmin_trans2, qmax_trans2, args[["dist"]])
    }
  },
  tf_compile_params = function(input, name_prefix = "") {
    ph <- self$get_placeholders()
    comp <- self$get_components()[[1L]]

    if ("offset" %in% names(ph)) {
      out <- list(
        offset = I_REALS$tf_make_layer(
          input = input,
          name = paste0(name_prefix, "offset")
        )
      )
      has_offset <- TRUE
    } else {
      out <- list()
      has_offset <- FALSE
    }

    if ("multiplier" %in% names(ph)) {
      out <- c(
        out,
        list(
          multiplier = I_REALS$tf_make_layer(
            input = input,
            name = paste0(name_prefix, "multiplier")
          )
        )
      )
      has_multiplier <- TRUE
    } else {
      has_multiplier <- FALSE
    }

    comp_inflater <- NULL

    if (length(ph$dist)) {
      dist_compiled <- comp$tf_compile_params(
        input = input,
        name_prefix = paste0(name_prefix, "dist_")
      )
      comp_inflater <- dist_compiled$output_inflater
      out <- c(out, dist_compiled$outputs)
    }

    list(
      outputs = out,
      output_inflater = eval(bquote(function(outputs) {
        if (!is.list(outputs)) outputs <- list(outputs)
        has_offset <- .(has_offset)
        has_multiplier <- .(has_multiplier)
        comp_inflater <- .(comp_inflater)

        out <- list()
        if (has_offset) {
          out$offset <- outputs[[1L]]
          outputs <- outputs[-1L]
        }
        if (has_multiplier) {
          out$multiplier <- outputs[[1L]]
          outputs <- outputs[-1L]
        }
        if (!is.null(comp_inflater)) {
          out$dist <- comp_inflater(outputs)
        } else {
          out$dist <- list()
        }
        out
      }))
    )
  },
  compile_sample = function() {
    ph <- names(self$get_placeholders())
    ph_offset <- "offset" %in% ph
    ph_multiplier <- "multiplier" %in% ph

    dist_sample <- self$get_components()[[1L]]$compile_sample()
    n_params_dist <- attr(dist_sample, "n_params")
    n_params <- as.integer(ph_offset) + as.integer(ph_multiplier) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    multiplier_expr <- if (ph_multiplier) {
      bquote(param_matrix[, .(n_params)])
    } else {
      self$default_params$multiplier
    }

    offset_expr <- if (ph_offset) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(n, param_matrix) {
        dist_sample(n, .(dist_param_expr)) * .(multiplier_expr) + .(offset_expr)
      })),
      n_params
    )
  },
  compile_density = function() {
    ph <- names(self$get_placeholders())
    ph_offset <- "offset" %in% ph
    ph_multiplier <- "multiplier" %in% ph

    dist_density <- self$get_components()[[1L]]$compile_density()
    n_params_dist <- attr(dist_density, "n_params")
    n_params <- as.integer(ph_offset) + as.integer(ph_multiplier) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    multiplier_expr <- if (ph_multiplier) {
      bquote(param_matrix[, .(n_params)])
    } else {
      self$default_params$multiplier
    }

    offset_expr <- if (ph_offset) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(x, param_matrix, log = FALSE) {
        mult <- .(multiplier_expr)
        dens <- dist_density((x - .(offset_expr)) / mult, .(dist_param_expr), log = log)
        if (log) dens - log(mult) else dens / mult
      })),
      n_params
    )
  },
  compile_probability = function() {
    ph <- names(self$get_placeholders())
    ph_offset <- "offset" %in% ph
    ph_multiplier <- "multiplier" %in% ph

    dist_probability <- self$get_components()[[1L]]$compile_probability()
    n_params_dist <- attr(dist_probability, "n_params")
    n_params <- as.integer(ph_offset) + as.integer(ph_multiplier) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    multiplier_expr <- if (ph_multiplier) {
      bquote(param_matrix[, .(n_params)])
    } else {
      self$default_params$multiplier
    }

    offset_expr <- if (ph_offset) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(p, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        dist_probability(
          (p - .(offset_expr)) / .(multiplier_expr),
          .(dist_param_expr),
          lower.tail = lower.tail,
          log.p = log.p
        )
      })),
      n_params
    )
  },
  compile_probability_interval = function() {
    ph <- names(self$get_placeholders())
    ph_offset <- "offset" %in% ph
    ph_multiplier <- "multiplier" %in% ph

    dist_probability <- self$get_components()[[1L]]$compile_probability_interval()
    n_params_dist <- attr(dist_probability, "n_params")
    n_params <- as.integer(ph_offset) + as.integer(ph_multiplier) + n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    multiplier_expr <- if (ph_multiplier) {
      bquote(param_matrix[, .(n_params)])
    } else {
      self$default_params$multiplier
    }

    offset_expr <- if (ph_offset) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
        dist_probability(
          (qmin - .(offset_expr)) / .(multiplier_expr),
          (qmax - .(offset_expr)) / .(multiplier_expr),
          .(dist_param_expr),
          log.p = log.p
        )
      })),
      n_params
    )
  },
  compile_quantile = function() {
    ph <- names(self$get_placeholders())
    ph_offset <- "offset" %in% ph
    ph_multiplier <- "multiplier" %in% ph

    dist_quantile <- self$get_components()[[1L]]$compile_quantile()
    n_params_dist <- attr(dist_quantile, "n_params")
    n_params <- as.integer(ph_offset) + as.integer(ph_multiplier) +
      n_params_dist

    dist_param_expr <- if (n_params_dist > 0L) {
      bquote(param_matrix[, 1L:.(n_params_dist), drop = FALSE])
    } else {
      NULL
    }

    multiplier_expr <- if (ph_multiplier) {
      bquote(param_matrix[, .(n_params)])
    } else {
      self$default_params$multiplier
    }

    offset_expr <- if (ph_offset) {
      bquote(param_matrix[, .(n_params_dist + 1L)])
    } else {
      self$default_params$offset
    }

    as_compiled_distribution_function(
      eval(bquote(function(p, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        dist_quantile(p, .(dist_param_expr), lower.tail = lower.tail, log.p = log.p) * .(multiplier_expr) +
          .(offset_expr)
      })),
      n_params
    )
  }
)

#' @export
fit_dist_start.TranslatedDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  res <- dist$get_placeholders()
  ph <- names(res)
  ph_dist <- length(res$dist) > 0L

  .assert_set(ph, "offset", "Translated")

  if (ph_dist) {
    obs_tr <- obs

    pp <- dist$get_params()

    # shift all data by offset (weights are unmodified)
    obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")] <- lapply(
      obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")],
      function(x) {
        (x - pp$offset) / pp$multiplier
      }
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

#' @export
fit_dist.TranslatedDistribution <- function(dist, obs, start, ...) {
  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start)

  if ("offset" %in% names(start) || "multiplier" %in% names(start)) {
    return(NextMethod("fit_dist"))
  }

  # If offset is fixed, transform observations and delegate to fit_dist() of
  # child Distribution.
  obs_tr <- obs

  pp <- dist$get_params()

  # shift all data by offset (weights are unmodified)
  obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")] <- lapply(
    obs_tr[, c("x", "xmin", "xmax", "tmin", "tmax")],
    function(x) {
      (x - pp$offset) / pp$multiplier
    }
  )

  res <- fit_dist(dist$get_components()[[1L]], obs_tr, start$dist, ...)
  # Wrap dist params
  # TODO transform history if trace is TRUE
  res$params <- list(dist = res$params)
  res
}

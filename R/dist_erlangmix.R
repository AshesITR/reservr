#' Erlang Mixture distribution
#'
#' @param shapes Shape parameters, a trunc_erlangmix fit, or `NULL` as a
#' placeholder.
#' @param scale Common scale parameter, or `NULL` as a placeholder.
#' @param probs Mixing probabilities, or `NULL` as a placeholder.
#'
#' @return An `ErlangMixtureDistribution` object.
#'
#' @examples
#' params <- list(scale = 1.0, probs = list(0.5, 0.3, 0.2), shapes = list(1L, 2L, 3L))
#' dist <- dist_erlangmix(vector("list", 3L))
#' x <- dist$sample(20, with_params = params)
#' d_emp <- dist_empirical(x, positive = TRUE)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = dist,
#'   with_params = list(
#'     theoretical = params
#'   ),
#'   .x = seq(1e-4, 5, length.out = 100)
#' )
#'
#' @export
#'
#' @family Distributions
dist_erlangmix <- function(shapes, scale = NULL, probs = NULL) {
  if (is.null(probs)) probs <- vector("list", length(shapes))

  ErlangMixtureDistribution$new(shapes = shapes, scale = scale, probs = probs)
}

ErlangMixtureDistribution <- distribution_class(
  name = "ErlangMixture",
  params = list(
    shapes = list(I_POSITIVE_REALS),
    scale = I_POSITIVE_REALS,
    probs = list(I_UNIT_INTERVAL)
  ),
  sample = function(n, params) {
    slot <- runif(n)
    k <- length(params$shapes)
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]
    slot <- rowSums(slot > probmat) + 1

    shapemat <- do.call(cbind, params$shapes)

    rgamma(
      n = n,
      shape = shapemat[cbind(seq_len(n), slot)],
      scale = params$scale
    )
  },
  density = function(x, log = FALSE, params) {
    params$probs <- lapply(params$probs, rep_len, length(x))
    probmat <- do.call(cbind, params$probs)
    probmat <- probmat / rowSums(probmat)

    densmat <- map_dbl_matrix(
      params$shapes,
      dgamma,
      length(x),
      x = x,
      scale = params$scale
    )

    res <- rowSums(densmat * probmat)
    if (log) res <- log(res)
    res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    params$probs <- lapply(params$probs, rep_len, length(q))
    probmat <- do.call(cbind, params$probs)

    cdfmat <- map_dbl_matrix(
      params$shapes,
      pgamma,
      length(q),
      q = q,
      lower.tail = lower.tail,
      scale = params$scale
    )

    res <- rowSums(cdfmat * probmat) / rowSums(probmat)
    if (log.p) res <- log(res)
    res
  },
  support = I_POSITIVE_REALS,
  get_dof = function() {
    super$get_dof() - 1L
  },
  get_param_constraints = function() {
    ph <- self$get_placeholders()
    if (length(ph$probs)) {
      function(params) {
        prob_mat <- do.call(cbind, params$probs)

        list(
          constraints = rowSums(prob_mat) - 1.0,
          jacobian = cbind(
            if (hasName(ph, "scale")) rep_len(0.0, nrow(prob_mat)) else NULL,
            array(1.0, dim = dim(prob_mat))
          )
        )
      }
    } else {
      NULL
    }
  },
  get_components = function() {
    scale <- private$.default_params$scale
    rate <- if (is.null(scale)) NULL else 1.0 / scale
    lapply(private$.default_params$shapes, function(shape) {
      dist_gamma(shape = shape, rate = rate)
    })
  },
  diff_density = function(x, vars, log, params) {
    res <- vars

    k <- length(private$.default_params$shapes)
    if (log) {
      dens <- self$density(x, with_params = params)
    }

    dgammas <- lapply(
      seq_len(k),
      function(k) {
        dgamma(x, scale = params$scale, shape = params$shapes[[k]])
      }
    )

    if (length(vars$probs)) {
      dprobs <- rowSums(
        map_dbl_matrix(
          seq_len(k),
          function(k) {
            -params$probs[[k]] * dgammas[[k]]
          },
          length(x)
        )
      )

      res$probs <- lapply(
        dgammas, `+`, dprobs
      )

      if (log) {
        res$probs <- lapply(res$probs, `/`, dens)
      }
    } else {
      res$probs <- NULL
    }

    if ("scale" %in% names(vars)) {
      res$scale <- rowSums(
        map_dbl_matrix(
          seq_len(k),
          function(k) {
            params$probs[[k]] *
              dgammas[[k]] *
              (x / params$scale - params$shapes[[k]]) / params$scale
          },
          length(x)
        )
      )

      if (log) res$scale <- res$scale / dens
    }

    if (length(vars$shapes)) {
      res$shapes <- lapply(
        seq_len(k),
        function(k) {
          out <- params$probs[[k]] *
            dgammas[[k]] * (
              log(x) -
                log(params$scale) -
                digamma(params$shapes[[k]])
            )

          if (log) out / dens else out
        }
      )
    } else {
      res$shapes <- NULL
    }

    res
  },
  tf_logdensity = function() {
    k <- length(self$get_components())
    function(x, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      shapes <- args[["shapes"]]
      shapes <- tf$reshape(shapes, keras::k_constant(as.integer(c(1, -1)), dtype = "int32"))
      scale <- args[["scale"]]
      scale <- tf$reshape(scale, keras::k_constant(as.integer(c(-1, 1)), dtype = "int32"))

      x_ok <- tf$math$is_finite(x) & x > K$zero
      x_safe <- tf$where(x_ok, x, K$one)

      tf$where(
        x_ok,
        tf$math$reduce_logsumexp(
          tf$where(
            x_ok[, tf$newaxis],
            log(probs) -
              shapes * log(scale) +
              (shapes - K$one) * log(x_safe[, tf$newaxis]) -
              lgamma(shapes),
            K$neg_inf
          ),
          axis = 1L
        ) - x_safe / scale[, 1L],
        K$neg_inf
      )
    }
  },
  tf_logprobability = function() {
    k <- length(self$get_components())
    function(qmin, qmax, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      shapes <- args[["shapes"]]
      shapes <- tf$reshape(shapes, keras::k_constant(as.integer(c(1, -1)), dtype = "int32"))
      scale <- args[["scale"]]
      scale <- tf$reshape(scale, keras::k_constant(as.integer(c(-1, 1)), dtype = "int32"))

      qmax_ok <- tf$math$is_finite(qmax) & qmax > K$zero
      qmin_ok <- qmin > K$zero
      qmax_finite <- tf$where(qmax_ok, qmax, K$two)
      qmin_finite <- tf$where(qmin_ok, qmin, qmax_finite / K$two)

      # When tmax is small compared to shapes, the gradient w.r.t. scale becomes
      # unstable (NaN). We catch and replace unstable shapes by 1. This should
      # stabilize the gradient; The resulting total contribution of the unstable
      # entries is effectively log(0) = -Inf.
      qmax_unstable <- (qmax <= K$zero)[, tf$newaxis] | (tf$math$igamma(
        shapes,
        qmax_finite[, tf$newaxis] / scale
      ) < keras::k_epsilon())

      qmin_unstable <- (qmin <= K$zero)[, tf$newaxis] | (tf$math$igamma(
        shapes,
        qmin_finite[, tf$newaxis] / scale
      ) < keras::k_epsilon())

      tf$where(
        qmax > K$zero,
        tf$math$reduce_logsumexp(
          tf$where(
            qmax_unstable & tf$math$is_finite(qmax)[, tf$newaxis],
            K$neg_inf,
            log(probs) +
              log(
                # igamma is the lower incomplete gamma function regularized with
                # factor 1 / Gamma(shape) i.e.
                # > tf$math$igamma(shape, x / scale) =
                # >  pgamma(x, shape = shape, scale = scale)
                tf$where(
                  qmax_ok[, tf$newaxis],
                  tf$math$igamma(
                    tf$where(qmax_unstable, K$one, shapes),
                    tf$where(qmax_unstable, K$inf, qmax_finite[, tf$newaxis] / scale)
                  ),
                  K$one
                ) -
                  tf$where(
                    qmin_ok[, tf$newaxis],
                    tf$math$igamma(
                      tf$where(qmin_unstable, K$one, shapes),
                      tf$where(qmin_unstable, K$zero, qmin_finite[, tf$newaxis] / scale)
                    ),
                    K$zero
                  )
              )
          ),
          axis = 1L
        ),
        K$neg_inf
      )
    }
  },
  tf_make_constants = function(with_params = list()) {
    check_installed("keras")
    params <- private$.make_params(with_params, 1)
    out <- list()
    if (length(params$probs) && !is.null(params$probs[[1L]])) {
      probs <- as.numeric(params$probs)
      out$probs <- keras::k_constant(probs / sum(probs), shape = c(1L, length(probs)))
    }
    if (length(params$shapes) && !is.null(params$shapes[[1L]])) {
      out$shapes <- keras::k_constant(as.numeric(params$shapes))
    }
    if (!is.null(params$scale)) {
      out$scale <- keras::k_constant(params$scale)
    }

    out
  },
  tf_compile_params = function(input, name_prefix = "") {
    check_installed("keras")
    ph <- self$get_placeholders()
    k <- length(self$get_components())
    out <- list()

    if (length(ph$shapes)) {
      out$shapes <- I_POSITIVE_REALS$tf_make_layer(
        input = input,
        name = paste0(name_prefix, "shapes"),
        size = k
      )
    }

    if (length(ph$probs)) {
      out$probs <- keras::layer_dense(
        input, units = k, activation = "softmax",
        name = paste0(name_prefix, "probs")
      )
    }

    if ("scale" %in% names(ph)) {
      out$scale <- I_POSITIVE_REALS$tf_make_layer(
        input = input,
        name = paste0(name_prefix, "scale")
      )
    }

    out_names <- names(out)

    list(
      outputs = out,
      output_inflater = eval(bquote(function(outputs) {
        if (!is.list(outputs)) outputs <- list(outputs)
        names(outputs) <- .(out_names)
        outputs
      }))
    )
  },
  compile_sample = function() {
    k <- length(self$get_params()$shapes)

    ph <- self$get_placeholders()
    ph_shapes <- length(ph$shapes) > 0L
    ph_scale <- "scale" %in% names(ph)
    ph_probs <- length(ph$probs) > 0L

    n_params <- (as.integer(ph_shapes) + as.integer(ph_probs)) * k + as.integer(ph_scale)

    scale_code <- if (ph_scale) substitute(param_matrix[, i_scale], list(
      i_scale = 1L + if (ph_shapes) k else 0L
    )) else self$default_params$scale

    probs_expr <- if (ph_shapes || ph_scale) {
      substitute(param_matrix[, i_prob:j_prob, drop = FALSE], list(
        i_prob = 1L + if (ph_shapes) k else 0L + as.integer(ph_scale),
        j_prob = k + if (ph_shapes) k else 0L + as.integer(ph_scale)
      ))
    } else {
      quote(param_matrix)
    }

    slot_code <- if (ph_probs) {
      substitute({
        slot <- runif(n)
        probmat <- matrixStats::rowCumsums(probs_expr)
        probmat <- probmat / probmat[, k]
        slot <- rowSums(slot > probmat) + 1
      }, list(k = k, probs_expr = probs_expr))
    } else {
      probs <- as.numeric(self$get_params()$probs)
      probs <- probs / sum(probs)
      substitute(
        slot <- sample.int(n = k, size = n, replace = TRUE, prob = probs),
        list(k = k, probs = probs)
      )
    }

    shape_code <- if (ph_shapes) {
      substitute(param_matrix[seq_len(n) + (slot - 1L) * n])
    } else {
      substitute(shapes[slot], list(shapes = as.numeric(self$get_params()$shapes)))
    }

    as_compiled_distribution_function(
      eval(substitute(function(n, param_matrix) {
        slot_code

        rgamma(
          n = n,
          shape = shape_code,
          scale = scale_code
        )
      }, list(slot_code = slot_code, shape_code = shape_code, scale_code = scale_code))),
      n_params = n_params
    )
  },
  compile_density = function() {
    k <- length(self$get_params()$shapes)

    ph <- self$get_placeholders()
    ph_shapes <- length(ph$shapes) > 0L
    ph_scale <- "scale" %in% names(ph)
    ph_probs <- length(ph$probs) > 0L

    n_params <- (as.integer(ph_shapes) + as.integer(ph_probs)) * k + as.integer(ph_scale)

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_scale) {
      scale_expr <- self$default_params$scale
    }

    if (!ph_shapes) {
      shapes_expr <- as.numeric(self$default_params$shapes)
    }

    dens_code <- if (!ph_probs && !ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_density_fixed_probs_scale_shape(
        x, log, .(prob_expr), .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs && !ph_scale) {
      bquote(drop(dist_erlangmix_density_fixed_probs_scale(
        x, param_matrix, log, .(prob_expr), .(scale_expr)
      )))
    } else if (!ph_probs && !ph_shapes) {
      bquote(drop(dist_erlangmix_density_fixed_probs_shape(
        x, param_matrix, log, .(prob_expr), .(shapes_expr)
      )))
    } else if (!ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_density_fixed_scale_shape(
        x, param_matrix, log, .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_erlangmix_density_fixed_probs(
        x, param_matrix, log, .(prob_expr)
      )))
    } else if (!ph_scale) {
      bquote(drop(dist_erlangmix_density_fixed_scale(
        x, param_matrix, log, .(scale_expr)
      )))
    } else if (!ph_shapes) {
      bquote(drop(dist_erlangmix_density_fixed_shape(
        x, param_matrix, log, .(shapes_expr)
      )))
    } else {
      bquote(drop(dist_erlangmix_density_free(
        x, param_matrix, log
      )))
    }

    as_compiled_distribution_function(
      eval(bquote(function(x, param_matrix, log = FALSE) {
        .(dens_code)
      })),
      n_params
    )
  },
  compile_probability = function() {
    k <- length(self$get_params()$shapes)

    ph <- self$get_placeholders()
    ph_shapes <- length(ph$shapes) > 0L
    ph_scale <- "scale" %in% names(ph)
    ph_probs <- length(ph$probs) > 0L

    n_params <- (as.integer(ph_shapes) + as.integer(ph_probs)) * k + as.integer(ph_scale)

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_scale) {
      scale_expr <- self$default_params$scale
    }

    if (!ph_shapes) {
      shapes_expr <- as.numeric(self$default_params$shapes)
    }

    prob_code <- if (!ph_probs && !ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_probability_fixed_probs_scale_shape(
        q, lower.tail, log.p, .(prob_expr), .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs && !ph_scale) {
      bquote(drop(dist_erlangmix_probability_fixed_probs_scale(
        q, param_matrix, lower.tail, log.p, .(prob_expr), .(scale_expr)
      )))
    } else if (!ph_probs && !ph_shapes) {
      bquote(drop(dist_erlangmix_probability_fixed_probs_shape(
        q, param_matrix, lower.tail, log.p, .(prob_expr), .(shapes_expr)
      )))
    } else if (!ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_probability_fixed_scale_shape(
        q, param_matrix, lower.tail, log.p, .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_erlangmix_probability_fixed_probs(
        q, param_matrix, lower.tail, log.p, .(prob_expr)
      )))
    } else if (!ph_scale) {
      bquote(drop(dist_erlangmix_probability_fixed_scale(
        q, param_matrix, lower.tail, log.p, .(scale_expr)
      )))
    } else if (!ph_shapes) {
      bquote(drop(dist_erlangmix_probability_fixed_shape(
        q, param_matrix, lower.tail, log.p, .(shapes_expr)
      )))
    } else {
      bquote(drop(dist_erlangmix_probability_free(
        q, param_matrix, lower.tail, log.p
      )))
    }

    as_compiled_distribution_function(
      eval(bquote(function(q, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        .(prob_code)
      })),
      n_params
    )
  },
  compile_probability_interval = function() {
    k <- length(self$get_params()$shapes)

    ph <- self$get_placeholders()
    ph_shapes <- length(ph$shapes) > 0L
    ph_scale <- "scale" %in% names(ph)
    ph_probs <- length(ph$probs) > 0L

    n_params <- (as.integer(ph_shapes) + as.integer(ph_probs)) * k + as.integer(ph_scale)

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_scale) {
      scale_expr <- self$default_params$scale
    }

    if (!ph_shapes) {
      shapes_expr <- as.numeric(self$default_params$shapes)
    }

    prob_code <- if (!ph_probs && !ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_iprobability_fixed_probs_scale_shape(
        qmin, qmax, log.p, .(prob_expr), .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs && !ph_scale) {
      bquote(drop(dist_erlangmix_iprobability_fixed_probs_scale(
        qmin, qmax, param_matrix, log.p, .(prob_expr), .(scale_expr)
      )))
    } else if (!ph_probs && !ph_shapes) {
      bquote(drop(dist_erlangmix_iprobability_fixed_probs_shape(
        qmin, qmax, param_matrix, log.p, .(prob_expr), .(shapes_expr)
      )))
    } else if (!ph_scale && !ph_shapes) {
      bquote(drop(dist_erlangmix_iprobability_fixed_scale_shape(
        qmin, qmax, param_matrix, log.p, .(scale_expr), .(shapes_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_erlangmix_iprobability_fixed_probs(
        qmin, qmax, param_matrix, log.p, .(prob_expr)
      )))
    } else if (!ph_scale) {
      bquote(drop(dist_erlangmix_iprobability_fixed_scale(
        qmin, qmax, param_matrix, log.p, .(scale_expr)
      )))
    } else if (!ph_shapes) {
      bquote(drop(dist_erlangmix_iprobability_fixed_shape(
        qmin, qmax, param_matrix, log.p, .(shapes_expr)
      )))
    } else {
      bquote(drop(dist_erlangmix_iprobability_free(
        qmin, qmax, param_matrix, log.p
      )))
    }

    as_compiled_distribution_function(
      eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
        .(prob_code)
      })),
      n_params
    )
  }
)

#' @export
fit_dist_start.ErlangMixtureDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph_scale <- "scale" %in% names(res)
  ph_shapes <- length(res$shapes) > 0L
  ph_probs <- length(res$probs) > 0L
  if (!ph_shapes) res$shapes <- NULL
  if (!ph_probs) res$probs <- NULL

  k <- length(dist$get_components())
  if (ph_probs && ph_shapes && ph_scale) {
    # TODO add weights for .trunc_erlangmix_init
    res <- .trunc_erlangmix_init(x, num_components = k, ...)
    res$shapes <- as.list(res$shapes)
    res$probs <- as.list(res$probs)
  } else if (ph_probs && ph_scale) {
    shapes <- as.integer(dist$get_params()$shapes)
    scale <- max(x) / shapes[k]
    bin <- .bincode(x, c(0, scale * shapes))
    # fix cases where numerically max(x) / shapes[m] * shapes[m] < max(x)
    bin[is.na(bin)] <- k

    probs <- weighted_tabulate(bin, obs$w, nbins = k) / sum(obs$w)

    res$probs <- as.list(probs)
    # Better starting value for scale via method of moments
    res$scale <- weighted.mean(x, obs$w) / weighted.mean(shapes, probs)
  } else if (ph_probs && ph_shapes) {
    full_init <- .trunc_erlangmix_init(x, num_components = k, ...)
    res$probs <- as.list(full_init$probs)
    res$shapes <- as.list(full_init$shapes)
  } else if (ph_shapes && ph_scale) {
    full_init <- .trunc_erlangmix_init(x, num_components = k, ...)
    res$shapes <- as.list(full_init$shapes)
    res$scale <- full_init$scale
  } else if (ph_probs) {
    pp <- dist$get_params()
    shapes <- as.integer(pp$shapes)
    scale <- pp$scale

    # FIXME respect truncation properly
    # cf fit_dist_start for Mixtures
    dmat <- dgamma_matrix(x, shapes, scale)
    dmat <- dmat / rowSums(dmat)

    res$probs <- as.list(as.numeric(t(dmat) %*% obs$w) / sum(obs$w))
  } else if (ph_shapes) {
    full_init <- .trunc_erlangmix_init(x, num_components = k, ...)
    res$shapes <- as.list(full_init$shapes)
  } else { # > ph_scale
    pp <- dist$get_params()
    shapes <- as.integer(pp$shapes)
    probs <- as.numeric(pp$probs)

    res$scale <- weighted.mean(x, obs$w) / weighted.mean(shapes, probs)
  }

  res
}

# TODO fit_erlang_mixture assumes probs to be tunable
# If probs is fixed, we should instead just fit with appropriately weighted
# observations
#' @include fit_erlang_mixture.R
#' @export
fit_dist.ErlangMixtureDistribution <- fit_erlang_mixture

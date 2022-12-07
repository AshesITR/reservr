#' Blended distribution
#'
#' @param dists A list of k >= 2 component Distributions.
#' @param probs k Mixture weight parameters
#' @param breaks k - 1 Centers of the blending zones.
#' `dists[i]` will blend into `dists[i + 1]` around `breaks[i]`.
#' @param bandwidths k - 1 Radii of the blending zones.
#' The i-th blending zone will begin at `breaks[i] - bandwidths[i]` and end at
#' `breaks[i] + bandwidths[i]`. A bandwidth of 0 corresponds to a hard cut-off,
#' i.e. a jump discontinuity in the density of the blended Distribution.
#'
#' @return A `BlendedDistribution` object.
#' @export
#'
#' @examples
#' bd <- dist_blended(
#'   list(
#'     dist_normal(mean = 0.0, sd = 1.0),
#'     dist_genpareto(u = 3.0, sigmau = 1.0, xi = 3.0)
#'   ),
#'   breaks = list(3.0),
#'   bandwidths = list(0.5),
#'   probs = list(0.9, 0.1)
#' )
#'
#' plot_distributions(
#'   bd,
#'   .x = seq(-3, 10, length.out = 100),
#'   plots = c("d", "p")
#' )
#'
#' @family Distributions
dist_blended <- function(dists, probs = NULL, breaks = NULL,
                         bandwidths = NULL) {
  breaks <- breaks %||% vector("list", length(dists) - 1L)
  bandwidths <- bandwidths %||% vector("list", length(dists) - 1L)
  probs <- probs %||% vector("list", length(dists))
  BlendedDistribution$new(dists = dists, probs = probs, breaks = breaks,
                          bandwidths = bandwidths)
}

BlendedDistribution <- distribution_class(
  name = "Blended",
  params = list(
    dists = list(),
    breaks = list(I_REALS),
    bandwidths = list(I_POSITIVE_REALS),
    probs = list(I_UNIT_INTERVAL)
  ),
  sample = function(n, params) {
    comps <- self$get_components()
    k <- length(comps)
    slot <- runif(n)
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]
    slot <- rowSums(slot > probmat) + 1

    u_mat <- do.call(cbind, params$breaks)
    eps_mat <- do.call(cbind, params$bandwidths)

    res <- numeric(n)
    for (i in seq_len(k)) {
      c_min <- if (i == 1L) -Inf else params$breaks[[i - 1L]]
      c_max <- if (i == k) Inf else params$breaks[[i]]
      comp_dist_trunc <- dist_trunc(comps[[i]], min = c_min, max = c_max)
      idx <- which(slot == i)
      trunc_smps <- comp_dist_trunc$sample(
        n = length(idx),
        with_params = list(dist = pick_params_at_idx(params$dists[[i]]$params, idx - 1L))
      )
      res[idx] <- blended_transition_inv(trunc_smps, u_mat[idx, , drop = FALSE], eps_mat[idx, , drop = FALSE], i)
    }

    res
  },
  density = function(x, log = FALSE, params) {
    comps <- self$get_components()
    k <- length(comps)

    u_mat <- do.call(cbind, params$breaks)
    eps_mat <- do.call(cbind, params$bandwidths)

    params$probs <- lapply(params$probs, rep_len, length(x))
    probmat <- do.call(cbind, params$probs)
    probmat <- probmat / rowSums(probmat)

    xtrans <- blended_transition(x, u_mat, eps_mat, .gradient = TRUE, .extend_na = TRUE)
    dblend <- attr(xtrans, "gradient")

    densmat <- map_dbl_matrix(
      seq_len(k),
      function(i) {
        in_range <- !is.na(xtrans[, i])

        if (i == k) {
          pmax <- 1.0
        } else {
          pmax <- comps[[i]]$probability(
            params$breaks[[i]][in_range], with_params = pick_params_at(params$dists[[i]]$params, in_range)
          )
        }

        if (i == 1L) {
          pmin <- 0.0
        } else {
          pmin <- comps[[i]]$probability(
            params$breaks[[i - 1L]][in_range], with_params = pick_params_at(params$dists[[i]]$params, in_range)
          )
        }

        res <- numeric(length(x))
        res[in_range] <- comps[[i]]$density(
          xtrans[in_range, i], with_params = pick_params_at(params$dists[[i]]$params, in_range)
        ) * dblend[in_range, i] / (pmax - pmin)

        res
      },
      length(x)
    )

    if (self$get_type() == "mixed") {
      is_disc <- map_lgl_matrix(
        seq_len(k),
        function(i) {
          in_range <- !is.na(xtrans[, i])

          res <- logical(length(x))
          res[in_range] <- comps[[i]]$is_discrete_at(
            xtrans[in_range, i], with_params = pick_params_at(params$dists[[i]]$params, in_range)
          )
          res
        },
        length(x)
      )

      has_disc <- matrixStats::rowAnys(is_disc)
      # Zero all non-discrete components if any component is discrete at that x
      densmat[!is_disc & has_disc] <- 0.0
    }

    res <- rowSums(probmat * densmat)
    if (log) res <- log(res)
    res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    comps <- self$get_components()
    k <- length(comps)

    params$probs <- lapply(params$probs, rep_len, length(q))
    probmat <- do.call(cbind, params$probs)
    probmat <- probmat / rowSums(probmat)

    u_mat <- do.call(cbind, params$breaks)
    eps_mat <- do.call(cbind, params$bandwidths)

    xtrans <- blended_transition(q, u_mat, eps_mat)

    cdfmat <- map_dbl_matrix(
      seq_len(k),
      function(i) {
        c_min <- if (i == 1L) rep_len(-Inf, length(q)) else params$breaks[[i - 1L]]
        c_max <- if (i == k) rep_len(Inf, length(q)) else params$breaks[[i]]

        pmax <- comps[[i]]$probability(
          c_max, with_params = params$dists[[i]]$params
        )
        pmin <- comps[[i]]$probability(
          c_min, with_params = params$dists[[i]]$params
        )

        comps[[i]]$probability(
          xtrans[, i], with_params = params$dists[[i]]$params
        ) / (pmax - pmin)
      },
      length(q)
    )

    res <- rowSums(probmat * cdfmat)
    if (!lower.tail) res <- 1.0 - res
    if (log.p) res <- log(res)
    res
  },
  support = function(x, params) {
    comps <- self$get_components()
    k <- length(comps)

    u_mat <- do.call(cbind, params$breaks)
    eps_mat <- do.call(cbind, params$bandwidths)

    xtrans <- blended_transition(x, u_mat, eps_mat, .extend_na = TRUE)

    matrixStats::rowAnys(map_lgl_matrix(
      seq_len(k),
      function(i) {
        res <- logical(length(x))
        in_range <- !is.na(xtrans[, i])
        res[in_range] <- comps[[i]]$is_in_support(
          x[in_range], with_params = pick_params_at(params$dists[[i]]$params, in_range)
        )
        res
      },
      length(x)
    ))
  },
  get_components = function() {
    private$.default_params$dists
  },
  get_param_constraints = function() {
    ph <- self$get_placeholders()
    constrs <- list()
    if (length(ph$probs)) {
      # Constrain colSums(do.call(cbind, probs)) == 1.0
      constrs <- c(constrs, function(params) {
        prob_mat <- do.call(cbind, params$probs)
        nms <- names(flatten_params(params))
        jac_full <- matrix(0, nrow = nrow(prob_mat), ncol = length(nms))
        jac_full[, grepl("^probs", nms)] <- 1.0

        list(
          constraints = rowSums(prob_mat) - 1.0,
          jacobian = jac_full
        )
      })
    }

    for (i in seq_along(ph$dists)) {
      if (length(ph$dists[[i]])) {
        # Any free parameters in component i
        comp_constr <- self$default_params$dists[[i]]$get_param_constraints()
        if (!is.null(comp_constr)) {
          constrs <- c(constrs, function(params) {
            comp_res <- comp_constr(params$dists[[i]])
            if (is.list(comp_res)) {
              nms <- names(flatten_params(params))
              jac_full <- matrix(0, nrow = length(comp_res$constraints), ncol = length(nms))
              jac_full[, grepl(paste0("^dists\\[", i, "\\]"), nms)] <- comp_res$jacobian

              list(
                constraints = comp_res$constraints,
                jacobian = jac_full
              )
            } else {
              comp_res
            }
          })
        }
      }
    }

    if (length(constrs) > 1L) {
      function(params) {
        all_constraints <- lapply(constrs, function(constr) constr(params))
        are_lists <- vapply(all_constraints, is.list, logical(1L))
        if (all(are_lists)) {
          # Jacobians for all constraints
          list(
            constraints = do.call(
              c,
              lapply(all_constraints, `[[`, i = "constraints")
            ),
            jacobian = do.call(
              rbind,
              lapply(all_constraints, `[[`, i = "jacobian")
            )
          )
        } else {
          all_constraints[are_lists] <- lapply(
            all_constraints[are_lists], `[[`, i = "constraints"
          )
          # No jacobians for some constraints
          do.call(c, all_constraints)
        }
      }
    } else if (length(constrs) == 1L) {
      constrs[[1L]]
    } else {
      # No constraints necessary
      NULL
    }
  },
  get_dof = function() {
    sdof <- super$get_dof()
    if (length(self$get_placeholders()$probs)) {
      sdof <- sdof - 1L
    }
    sdof
  },
  get_type = function() {
    if (all(purrr::map_lgl(self$get_components(), ~.$is_discrete()))) {
      "discrete"
    } else if (all(purrr::map_lgl(self$get_components(), ~.$is_continuous()))) {
      "continuous"
    } else {
      "mixed"
    }
  },
  has_capability = function(caps) {
    super$has_capability(caps) &
      matrixStats::rowAlls(
        map_lgl_matrix(
          self$get_components(),
          function(comp) comp$has_capability(caps),
          length(caps)
        )
      )
  },
  is_discrete = function(x, params) {
    if (!self$is_continuous()) {
      comps <- self$get_components()
      k_types <- purrr::map_chr(comps, ~.$get_type())

      u_mat <- do.call(cbind, params$breaks)
      eps_mat <- do.call(cbind, params$bandwidths)

      xtrans <- blended_transition(x, u_mat, eps_mat, .extend_na = TRUE)

      matrixStats::rowAnys(map_lgl_matrix(
        which(k_types != "continuous"),
        function(i) {
          res <- logical(length(x))
          in_range <- !is.na(xtrans[, i])
          res[in_range] <- comps[[i]]$is_discrete_at(
            x[in_range], with_params = pick_params_at(params$dists[[i]]$params, in_range)
          )
          res
        },
        length(x)
      ))
    } else {
      super$is_discrete_at(x, params)
    }
  },
  tf_make_constants = function(with_params = list()) {
    check_installed("keras")
    params <- private$.make_params(with_params, 1)
    out <- list()
    if (length(params$probs) && !is.null(params$probs[[1L]])) {
      probs <- as.numeric(params$probs)
      out$probs <- keras::k_constant(probs / sum(probs), shape = list(1L, length(params$probs)))
    }
    if (length(params$breaks) && !is.null(params$breaks[[1L]])) {
      out$breaks <- keras::k_constant(
        as.numeric(params$breaks), shape = list(length(params$breaks))
      )
    }
    if (length(params$bandwidths) && !is.null(params$bandwidths[[1L]])) {
      out$bandwidths <- keras::k_constant(
        as.numeric(params$bandwidths), shape = list(length(params$bandwidths))
      )
    }
    out$dists <- lapply(
      params$dists,
      function(d) d$dist$tf_make_constants(d$params)
    )
    out
  },
  tf_compile_params = function(input, name_prefix = "") {
    ph <- self$get_placeholders()
    comps <- self$get_components()
    k <- length(comps)
    if (length(ph$probs)) {
      out <- list(
        probs = keras::layer_dense(
          input, units = k, activation = "softmax",
          name = paste0(name_prefix, "probs")
        )
      )
      out_indices <- 0L
    } else {
      out <- list()
      out_indices <- integer()
    }

    if (length(ph$breaks)) {
      out$breaks <- keras::layer_dense(
        input, units = k - 1L, activation = "linear",
        name = paste0(name_prefix, "breaks")
      )
      out_indices <- c(out_indices, -1L)
    }

    if (length(ph$bandwidths)) {
      out$bandwidths <- keras::layer_dense(
        input, units = k - 1L, activation = "softplus",
        name = paste0(name_prefix, "bandwidths")
      )
      out_indices <- c(out_indices, -2L)
    }

    comp_inflaters <- vector("list", k)
    for (i in seq_len(k)) {
      if (length(ph$dists[[i]])) {
        comp_compiled <- comps[[i]]$tf_compile_params(
          input = input,
          name_prefix = paste0(name_prefix, "dists_", i, "_")
        )
        comp_inflaters[[i]] <- comp_compiled$output_inflater
        out_indices <- c(out_indices, rep_len(i, length(comp_compiled$outputs)))
        out <- c(out, comp_compiled$outputs)
      }
    }

    list(
      outputs = out,
      output_inflater = eval(bquote(function(outputs) {
        if (!is.list(outputs)) outputs <- list(outputs)
        idx <- .(out_indices)
        inflaters <- .(comp_inflaters)
        out <- list()
        if (0L %in% idx) {
          out$probs <- outputs[[which(idx == 0L)]]
        }
        if (-1L %in% idx) {
          out$breaks <- outputs[[which(idx == -1L)]]
        }
        if (-2L %in% idx) {
          out$bandwidths <- outputs[[which(idx == -2L)]]
        }
        out$dists <- lapply(seq_along(inflaters), function(i) {
          if (i %in% idx) {
            inflaters[[i]](outputs[idx == i])
          } else {
            list()
          }
        })
        out
      }))
    )
  },
  tf_is_discrete_at = function() {
    if (self$is_continuous()) return(super$tf_is_discrete_at())

    comp_discretes <- lapply(self$get_components(), function(comp) {
      comp$tf_is_discrete_at()
    })

    function(x, args) {
      # TODO correctly handle discreteness in blending intervals
      # and out-of-interval discrete points
      are_discrete <- tf$stack(
        lapply(seq_along(comp_discretes),
               function(i) comp_discretes[[i]](x, args[["dists"]][[i]]))
      )
      tf$math$reduce_any(are_discrete, axis = 1L)
    }
  },
  tf_logdensity = function() {
    comps <- self$get_components()
    comps_logdens <- lapply(comps, function(comp) comp$tf_logdensity())
    comps_logprob <- lapply(comps, function(comp) comp$tf_logprobability())
    k <- length(comps)

    if (k != 2L) {
      stop("Currenty, only 2-component blended distributions ",
           "are supported by tf_logdensity.", call. = FALSE)
    }

    function(x, args) {
      tf <- tensorflow::tf
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      dist_args <- args[["dists"]]
      u <- args[["breaks"]]
      eps <- args[["bandwidths"]]

      # TODO make compatible with u trainable
      imin <- tf$concat(list(keras::k_constant(-Inf, shape = 1L), u), axis = 0L)
      imax <- tf$concat(list(u, keras::k_constant(Inf, shape = 1L)), axis = 0L)

      log_probs <- log(probs)

      left_tail <- x < u - eps
      right_tail <- x > u + eps

      x_blend_safe <- tf$where(left_tail | right_tail, u, x)

      blend_angle <- K$one_half * (x_blend_safe - u) / eps
      blend_delta <- eps * cospi(blend_angle) / K$pi

      diff_blend_delta <- K$one_half * sinpi(blend_angle)

      x_left <- tf$where(
        left_tail,
        x,
        tf$where(
          right_tail,
          x + eps,
          K$one_half * (x + u - eps) + blend_delta
        )
      )

      x_right <- tf$where(
        left_tail,
        x - eps,
        tf$where(
          right_tail,
          x,
          K$one_half * (x + u + eps) - blend_delta
        )
      )

      imin_left <- tf$where(
        right_tail[, tf$newaxis],
        K$neg_inf,
        imin
      )
      imin_right <- tf$where(
        left_tail[, tf$newaxis],
        K$neg_inf,
        imin
      )
      imax_left <- tf$where(
        right_tail[, tf$newaxis],
        K$inf,
        imax
      )
      imax_right <- tf$where(
        left_tail[, tf$newaxis],
        K$inf,
        imax
      )

      logdiff_left <- tf$where(
        left_tail, K$zero,
        tf$where(
          right_tail, K$neg_inf,
          log(K$one_half - diff_blend_delta)
        )
      )
      logdiff_right <- tf$where(
        left_tail, K$neg_inf,
        tf$where(
          right_tail, K$zero,
          log(K$one_half + diff_blend_delta)
        )
      )

      compdens_left <- comps_logdens[[1L]](x_left, dist_args[[1L]])
      compdens_right <- comps_logdens[[2L]](x_right, dist_args[[2L]])
      compprob_left <- comps_logprob[[1L]](imin_left[, 1L], imax_left[, 1L], dist_args[[1L]])
      compprob_right <- comps_logprob[[2L]](imin_right[, 2L], imax_right[, 2L], dist_args[[2L]])

      dleft_ok <- compdens_left > K$neg_inf
      dleft_safe <- tf$where(dleft_ok, compdens_left, K$zero) - tf$where(dleft_ok, compprob_left, K$zero)
      dright_ok <- compdens_right > K$neg_inf
      dright_safe <- tf$where(dright_ok, compdens_right, K$zero) - tf$where(dright_ok, compprob_right, K$zero)

      logdens_comps <- tf$stack(list(
        tf$where(dleft_ok, log_probs[, 1L] + logdiff_left + dleft_safe, K$neg_inf),
        tf$where(dright_ok, log_probs[, 2L] + logdiff_right + dright_safe, K$neg_inf)
      ), axis = 1L)

      tf$where(
        left_tail,
        tf$where(dleft_ok, log_probs[, 1L] + dleft_safe, K$neg_inf),
        tf$where(
          right_tail,
          tf$where(dright_ok, log_probs[, 2L] + dright_safe, K$neg_inf),
          tf$reduce_logsumexp(logdens_comps, axis = 1L)
        )
      )
    }
  },
  tf_logprobability = function() {
    comps_logprob <- lapply(
      self$get_components(),
      function(comp) comp$tf_logprobability()
    )

    if (length(comps_logprob) != 2L) {
      stop("Currenty, only 2-component blended distributions ",
           "are supported by tf_logdensity.", call. = FALSE)
    }

    tf_blend <- function(x, imin, epsmin, imax, epsmax) {
      tf$where(
        x <= imin - epsmin,
        imin,
        tf$where(
          x < imin + epsmin,
          K$one_half * (x + imin + epsmin) -
            (epsmin * cospi(K$one_half * (x - imin) / epsmin)) / K$pi,
          tf$where(
            x <= imax - epsmax,
            x,
            tf$where(
              x < imax + epsmax,
              K$one_half * (x + imax - epsmax) +
                (epsmax * cospi(K$one_half * (x - imax) / epsmax)) / K$pi,
              imax
            )
          )
        )
      )
    }

    function(qmin, qmax, args) {
      k <- length(comps_logprob)
      tf <- tensorflow::tf

      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      dist_args <- args[["dists"]]
      u <- args[["breaks"]]
      eps <- args[["bandwidths"]]

      # TODO make compatible with u trainable
      imin <- tf$concat(list(keras::k_constant(-Inf, shape = 1L), u), axis = 0L)
      epsmin <- tf$concat(list(keras::k_constant(0.0, shape = 1L), eps), axis = 0L)
      imax <- tf$concat(list(u, keras::k_constant(Inf, shape = 1L)), axis = 0L)
      epsmax <- tf$concat(list(eps, keras::k_constant(0.0, shape = 1L)), axis = 0L)

      log_probs <- log(probs)

      qmin_trans <- tf$stack(lapply(seq_len(k), function(i) {
        tf_blend(
          qmin,
          imin[i], epsmin[i],
          imax[i], epsmax[i]
        )
      }), axis = 1L)

      qmax_trans <- tf$stack(lapply(seq_len(k), function(i) {
        tf_blend(
          qmax,
          imin[i], epsmin[i],
          imax[i], epsmax[i]
        )
      }), axis = 1L)

      comp_logprobs <- tf$stack(lapply(seq_len(k), function(i) {
        prob <- comps_logprob[[i]](
          qmin_trans[, i],
          qmax_trans[, i],
          dist_args[[i]]
        )

        sc <- tf$where(
          tf$math$is_finite(prob),
          comps_logprob[[i]](
            tf$broadcast_to(imin[i], tf$shape(qmin)),
            tf$broadcast_to(imax[i], tf$shape(qmax)),
            dist_args[[i]]
          ),
          K$zero
        )

        prob - sc
      }), axis = 1L)

      tf$reduce_logsumexp(tf$where(comp_logprobs > K$neg_inf, log_probs + comp_logprobs, K$neg_inf), axis = 1L)
    }
  },
  compile_sample = function() {
    comps <- self$get_components()
    k <- length(comps)
    ph <- self$get_placeholders()
    ph_probs <- length(ph$probs) > 0L
    ph_u <- length(ph$breaks) > 0L
    ph_eps <- length(ph$bandwidths) > 0L

    comp_prob <- lapply(comps, function(comp) comp$compile_probability())
    comp_prob_int <- lapply(comps, function(comp) comp$compile_probability_interval())
    comp_quantile <- lapply(comps, function(comp) comp$compile_quantile())

    comp_param_counts <- vapply(comp_prob, function(fun) attr(fun, "n_params"), integer(1L))
    comp_param_ends <- cumsum(comp_param_counts)
    comp_param_starts <- comp_param_ends - comp_param_counts + 1L

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L + if (ph_u) k - 1L else 0L + if (ph_eps) k - 1L else 0L

    slot_code <- if (ph_probs) {
      bquote({
        slot <- runif(n)
        probmat <- matrixStats::rowCumsums(
          param_matrix[, .(sum(comp_param_counts) + 1L):.(sum(comp_param_counts) + k), drop = FALSE]
        )
        probmat <- probmat / probmat[, .(k)]
        slot <- rowSums(slot > probmat) + 1
      })
    } else {
      probs <- as.numeric(self$get_params()$probs)
      probs <- probs / sum(probs)
      bquote(slot <- sample.int(n = .(k), size = n, replace = TRUE, prob = .(probs)))
    }

    sampling_code <- bquote({
      res <- numeric(n)
      num_samples <- tabulate(slot, .(k))
    })
    for (i in seq_len(k)) {
      comp_param_expr <- if (comp_param_counts[i] > 0L) {
        bquote(param_matrix[idx, .(comp_param_starts[i]):.(comp_param_ends[i]), drop = FALSE])
      } else {
        NULL
      }

      break_min_expr <- if (i == 1L) {
        -Inf
      } else if (ph_u) {
        bquote(param_matrix[idx, .(sum(comp_param_counts) + i - 1L)])
      } else {
        self$default_params$breaks[[i - 1L]]
      }

      epsilon_min_expr <- if (i == 1L) {
        0.0
      } else if (ph_eps) {
        bquote(param_matrix[idx, .(sum(comp_param_counts) + k + i - 2L)])
      } else {
        self$default_params$bandwidths[[i - 1L]]
      }

      break_max_expr <- if (i == k) {
        Inf
      } else if (ph_u) {
        bquote(param_matrix[idx, .(sum(comp_param_counts) + i)])
      } else {
        self$default_params$breaks[[i]]
      }

      epsilon_max_expr <- if (i == k) {
        0.0
      } else if (ph_eps) {
        bquote(param_matrix[idx, .(sum(comp_param_counts) + k + i - 1L)])
      } else {
        self$default_params$bandwiths[[i]]
      }

      sampling_code[[i + 3L]] <- bquote({
        idx <- slot == .(i)
        q_max <- comp_prob[[.(i)]](.(break_max_expr), .(comp_param_expr))
        q_min <- q_max - comp_prob_int[[.(i)]](.(break_min_expr), .(break_max_expr), .(comp_param_expr))
        q <- runif(num_samples[.(i)], min = q_min, max = q_max)
        trunc_smps <- comp_quantile[[.(i)]](q, .(comp_param_expr))
        res[idx] <- blended_transition_finv(
          trunc_smps,
          .(break_min_expr),
          .(break_max_expr),
          .(epsilon_min_expr),
          .(epsilon_max_expr),
          blend_left = .(i > 1L),
          blend_right = .(i < k)
        )
      })
    }

    as_compiled_distribution_function(
      eval(bquote(function(n, param_matrix) {
        .(slot_code)
        .(sampling_code)
        res
      })),
      n_params
    )
  },
  compile_density = function() {
    comps <- self$get_components()
    k <- length(comps)
    ph <- self$get_placeholders()
    ph_probs <- length(ph$probs) > 0L
    ph_u <- length(ph$breaks) > 0L
    ph_eps <- length(ph$bandwidths) > 0L

    comp_dens <- lapply(comps, function(comp) comp$compile_density())
    comp_prob_int <- lapply(comps, function(comp) comp$compile_probability_interval())

    comp_param_counts <- vapply(comp_dens, function(fun) attr(fun, "n_params"), integer(1L))

    comp_types <- vapply(comps, function(comp) comp$get_type(), character(1L))
    # TODO implement mixed type too
    stopifnot(all(comp_types %in% c("discrete", "continuous")))
    comp_discrete <- as.integer(comp_types == "discrete")

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L + if (ph_u) k - 1L else 0L + if (ph_eps) k - 1L else 0L

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_u) {
      u_expr <- as.numeric(self$default_params$breaks)
    }

    if (!ph_eps) {
      eps_expr <- as.numeric(self$default_params$bandwidths)
    }

    dens_code <- if (!ph_probs && !ph_u && !ph_eps) {
      bquote(drop(dist_blended_density_fixed_probs_breaks_eps(
        x, param_matrix, log, .(comp_param_counts),
        comp_dens, comp_prob_int, .(comp_discrete), .(prob_expr), .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs && !ph_u) {
      bquote(drop(dist_blended_density_fixed_probs_breaks(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(prob_expr), .(u_expr)
      )))
    } else if (!ph_probs && !ph_eps) {
      bquote(drop(dist_blended_density_fixed_probs_eps(
        x, param_matrix, log,
        .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(prob_expr), .(eps_expr)
      )))
    } else if (!ph_u && !ph_eps) {
      bquote(drop(dist_blended_density_fixed_breaks_eps(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_blended_density_fixed_probs(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(prob_expr)
      )))
    } else if (!ph_u) {
      bquote(drop(dist_blended_density_fixed_breaks(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(u_expr)
      )))
    } else if (!ph_eps) {
      bquote(drop(dist_blended_density_fixed_eps(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete), .(eps_expr)
      )))
    } else {
      bquote(drop(dist_blended_density_free(
        x, param_matrix, log, .(comp_param_counts), comp_dens, comp_prob_int, .(comp_discrete)
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
    comps <- self$get_components()
    k <- length(comps)
    ph <- self$get_placeholders()
    ph_probs <- length(ph$probs) > 0L
    ph_u <- length(ph$breaks) > 0L
    ph_eps <- length(ph$bandwidths) > 0L

    comp_prob_int <- lapply(comps, function(comp) comp$compile_probability_interval())

    comp_param_counts <- vapply(comp_prob_int, function(fun) attr(fun, "n_params"), integer(1L))

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L + if (ph_u) k - 1L else 0L + if (ph_eps) k - 1L else 0L

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_u) {
      u_expr <- as.numeric(self$default_params$breaks)
    }

    if (!ph_eps) {
      eps_expr <- as.numeric(self$default_params$bandwidths)
    }

    prob_code <- if (!ph_probs && !ph_u && !ph_eps) {
      bquote(drop(dist_blended_probability_fixed_probs_breaks_eps(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs && !ph_u) {
      bquote(drop(dist_blended_probability_fixed_probs_breaks(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(u_expr)
      )))
    } else if (!ph_probs && !ph_eps) {
      bquote(drop(dist_blended_probability_fixed_probs_eps(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(eps_expr)
      )))
    } else if (!ph_u && !ph_eps) {
      bquote(drop(dist_blended_probability_fixed_breaks_eps(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_blended_probability_fixed_probs(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr)
      )))
    } else if (!ph_u) {
      bquote(drop(dist_blended_probability_fixed_breaks(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(u_expr)
      )))
    } else if (!ph_eps) {
      bquote(drop(dist_blended_probability_fixed_eps(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int, .(eps_expr)
      )))
    } else {
      bquote(drop(dist_blended_probability_free(
        q, param_matrix, lower.tail, log.p, .(comp_param_counts), comp_prob_int
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
    comps <- self$get_components()
    k <- length(comps)
    ph <- self$get_placeholders()
    ph_probs <- length(ph$probs) > 0L
    ph_u <- length(ph$breaks) > 0L
    ph_eps <- length(ph$bandwidths) > 0L

    comp_prob_int <- lapply(comps, function(comp) comp$compile_probability_interval())

    comp_param_counts <- vapply(comp_prob_int, function(fun) attr(fun, "n_params"), integer(1L))

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L + if (ph_u) k - 1L else 0L + if (ph_eps) k - 1L else 0L

    if (!ph_probs) {
      prob_expr <- as.numeric(self$default_params$probs)
    }

    if (!ph_u) {
      u_expr <- as.numeric(self$default_params$breaks)
    }

    if (!ph_eps) {
      eps_expr <- as.numeric(self$default_params$bandwidths)
    }

    prob_code <- if (!ph_probs && !ph_u && !ph_eps) {
      bquote(drop(dist_blended_iprobability_fixed_probs_breaks_eps(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs && !ph_u) {
      bquote(drop(dist_blended_iprobability_fixed_probs_breaks(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(u_expr)
      )))
    } else if (!ph_probs && !ph_eps) {
      bquote(drop(dist_blended_iprobability_fixed_probs_eps(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr), .(eps_expr)
      )))
    } else if (!ph_u && !ph_eps) {
      bquote(drop(dist_blended_iprobability_fixed_breaks_eps(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(u_expr), .(eps_expr)
      )))
    } else if (!ph_probs) {
      bquote(drop(dist_blended_iprobability_fixed_probs(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(prob_expr)
      )))
    } else if (!ph_u) {
      bquote(drop(dist_blended_iprobability_fixed_breaks(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(u_expr)
      )))
    } else if (!ph_eps) {
      bquote(drop(dist_blended_iprobability_fixed_eps(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int, .(eps_expr)
      )))
    } else {
      bquote(drop(dist_blended_iprobability_free(
        qmin, qmax, param_matrix, log.p, .(comp_param_counts), comp_prob_int
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

# TODO fit_blended assumes weights to be tunable
# If weights is fixed, we should instead just fit with appropriately weighted
# observations
#' @include fit_blended.R
#' @export
fit_dist.BlendedDistribution <- fit_blended

# TODO fit_dist_start assumes breaks and bandwidths to be known
# If they are missing, heuristics should be used as starting values to not
# throw an error.

#' @export
fit_dist_start.BlendedDistribution <- function(dist, obs, dists_start = NULL,
                                               ...) {
  obs <- as_trunc_obs(obs)
  res <- dist$get_placeholders()
  ph_dists <- lengths(res$dists) > 0L
  ph_probs <- length(res$probs) > 0L
  ph_breaks <- length(res$breaks) > 0L
  ph_bandwidths <- length(res$bandwidths) > 0L

  blend_comps <- dist$get_components()
  blend <- dist$get_params()

  k <- length(blend_comps)
  n <- nrow(obs)

  if (ph_breaks) {
    breaks_approx <- weighted_quantile(
      x = obs$x,
      w = obs$w,
      probs = seq_len(k - 1) / (k + 1)
    )
    res$breaks <- as.list(breaks_approx)
    blend$breaks <- res$breaks
  }

  if (ph_bandwidths) {
    eps_default <- min(diff(c(min(obs$x), unlist(res$breaks), max(obs$x)))) /
      3.0
    res$bandwidths <- as.list(rep(eps_default, k - 1))
    blend$bandwidths <- res$bandwidths
  }

  comp_supp <- map_lgl_matrix(
    seq_len(k),
    function(i) {
      if (i == 1L) {
        obs$xmin <= blend$breaks[[1L]] + blend$bandwidths[[1L]]
      } else if (i == k) {
        obs$xmax >= blend$breaks[[k - 1L]] - blend$bandwidths[[k - 1L]]
      } else {
        obs$xmax >= blend$breaks[[i - 1L]] - blend$bandwidths[[i - 1L]] &
          obs$xmin <= blend$breaks[[i]] + blend$bandwidths[[i]]
      }
    },
    n
  )

  u_mat <- unlist(blend$breaks)
  eps_mat <- unlist(blend$bandwidths)
  obs_trans <- blended_transition(obs$x, u = u_mat, eps = eps_mat,
                                  .gradient = TRUE, .extend_na = TRUE)
  obs_trans_gradient <- attr(obs_trans, "gradient")
  obs_trans_min <- blended_transition(obs$xmin, u = u_mat, eps = eps_mat)
  obs_trans_max <- blended_transition(obs$xmax, u = u_mat, eps = eps_mat)
  tmin_comp <- blended_transition(obs$tmin, u_mat, eps_mat)
  tmax_comp <- blended_transition(obs$tmax, u_mat, eps_mat)

  # Precautions against (small) numerical error in blended_transition:
  obs_trans_min <- pmin(obs_trans_min, obs_trans, na.rm = TRUE)
  obs_trans_max <- pmax(obs_trans_max, obs_trans, na.rm = TRUE)
  tmin_comp <- pmin(tmin_comp, obs_trans_min)
  tmax_comp <- pmax(tmax_comp, obs_trans_max)

  res$dists[!ph_dists] <- rep_len(list(list()), sum(!ph_dists))
  if (any(ph_dists)) {
    if (!is.null(dists_start)) {
      for (comp in which(ph_dists)) {
        res$dists[[comp]] <- dists_start[[comp]]
      }
    } else {
      for (comp in which(ph_dists)) {
        # Transform blended parts and update truncation bounds
        if (comp > 1L) {
          comp_min <- blend$breaks[[comp - 1L]]
        } else {
          comp_min <- -Inf
        }

        if (comp < k) {
          comp_max <- blend$breaks[[comp]]
        } else {
          comp_max <- Inf
        }

        i_keep <- comp_supp[, comp]
        obs_comp <- trunc_obs(
          x = obs_trans[i_keep, comp],
          xmin = obs_trans_min[i_keep, comp],
          xmax = obs_trans_max[i_keep, comp],
          tmin = tmin_comp[i_keep, comp],
          tmax = tmax_comp[i_keep, comp]
        )
        obs_comp <- truncate_obs(obs_comp, comp_min, comp_max)

        res$dists[[comp]] <- fit_dist_start(
          dist = blend_comps[[comp]],
          obs = obs_comp,
          ...
        )
      }
    }
  }

  if (ph_probs) {
    # Compute partial component densities
    densmat <- map_dbl_matrix(
      seq_len(k),
      function(i) {
        c_obs <- !is.na(obs_trans[, i])
        c_cens <- !c_obs & comp_supp[, i]

        if (i == k) {
          pmax <- 1.0
        } else {
          pmax <- blend_comps[[i]]$probability(
            blend$breaks[[i]], with_params = res$dists[[i]]
          )
        }

        if (i == 1L) {
          pmin <- 0.0
        } else {
          pmin <- blend_comps[[i]]$probability(
            blend$breaks[[i - 1L]], with_params = res$dists[[i]]
          )
        }

        dens <- numeric(n)
        dens[c_obs] <- blend_comps[[i]]$density(
          obs_trans[c_obs, i], with_params = res$dists[[i]]
        ) * obs_trans_gradient[c_obs, i] / (pmax - pmin)
        dens[c_cens] <- (
          blend_comps[[i]]$probability(
            obs_trans_max[c_cens, i], with_params = res$dists[[i]]
          ) - blend_comps[[i]]$probability(
            obs_trans_min[c_cens, i], with_params = res$dists[[i]]
          )
        ) / (pmax - pmin)

        if (!blend_comps[[i]]$is_continuous()) {
          is_disc <- blend_comps[[i]]$is_discrete_at(
            obs_trans_min[c_cens, i], with_params = res$dists[[i]]
          )

          if (any(is_disc)) {
            dens[c_cens][is_disc] <- dens[c_cens][is_disc] + blend_comps[[i]]$density(
              obs_trans_min[c_cens, i][is_disc], with_params = res$dists[[i]]
            )
          }
        }

        dens
      },
      n
    )

    if (dist$get_type() == "mixed") {
      is_disc <- map_lgl_matrix(
        seq_len(k),
        function(i) {
          in_range <- !is.na(obs_trans[, i])

          disc <- logical(n)
          disc[in_range] <- blend_comps[[i]]$is_discrete_at(
            obs_trans[in_range, i], with_params = res$dists[[i]]
          )
          disc
        },
        n
      )

      has_disc <- matrixStats::rowAnys(is_disc)
      # Zero all non-discrete components if any component is discrete at that x
      densmat[!is_disc & has_disc] <- 0.0
    }

    res$probs <- local({
      p0 <- colSums(obs$w * (densmat / rowSums(densmat)))
      as.list(p0 / sum(p0))
    })
  } else {
    res$probs <- NULL
  }

  res
}

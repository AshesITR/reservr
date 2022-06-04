#' Mixture distribution
#'
#' Parameters of mixing components can be overridden with
#' `with_params = list(dists = list(..., ..., ...))`.
#' #' Mixing probabilites can be overridden with
#' `with_params = list(probs = list(..., ..., ...))`.
#' The **number of components** cannot be overridden.
#'
#' Does **not** support the `quantile()` capability!
#'
#' @param dists A list of mixing distributions.
#' May contain placeholders and duplicates.
#' @param probs A list of mixing probabilities with the same length as `dists`.
#' They are normalized to sum to one and `NULL` can be used as a placeholder
#' within probs.
#' To reduce the number of required parameters, probs should at least be partly
#' specified (`probs = list(NULL, NULL, ..., 1)` with k - 1 `NULL`s where k is
#' the number of mixing components).
#'
#' @return A `MixtureDistribution` object.
#' @export
#'
#' @examples
#'
#' # A complicated way to define a uniform distribution on \[0, 2\]
#' dist_mixture(
#'   dists = list(
#'     dist_uniform(min = 0, max = 1),
#'     dist_uniform(min = 1, max = 2)
#'   ),
#'   probs = list(0.5, 0.5)
#' )
#'
#' @family Distributions
dist_mixture <- function(dists = list(), probs = NULL) {
  if (is.null(probs)) probs <- vector("list", length(dists))
  MixtureDistribution$new(dists = dists, probs = probs)
}

MixtureDistribution <- distribution_class(
  name = "Mixture",
  params = list(
    dists = list(NULL),
    probs = list(I_UNIT_INTERVAL)
  ),
  sample = function(n, params) {
    slot <- runif(n)
    k <- length(params$dists)
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]
    slot <- rowSums(slot > probmat) + 1

    pick_idx <- function(p, idx) {
      if (is.list(p)) {
        lapply(p, pick_idx, idx = idx)
      } else if (is.null(p) || length(p) == 1) {
        p
      } else {
        p[idx]
      }
    }

    res <- numeric(n)
    for (i in seq_len(k)) {
      idx <- which(slot == i)
      res[idx] <- params$dists[[i]]$dist$sample(
        n = length(idx),
        with_params = pick_idx(params$dists[[i]]$params, idx)
      )
    }

    res
  },
  density = function(x, log = FALSE, params) {
    params$probs <- lapply(params$probs, rep_len, length(x))
    probmat <- do.call(cbind, params$probs)
    probmat <- probmat / rowSums(probmat)

    densmat <- map_dbl_matrix(
      params$dists,
      ~.$dist$density(
        x = x,
        with_params = .$params
      ),
      length(x)
    )

    comp_types <- purrr::map_chr(self$get_components(), ~.$get_type())

    if ("mixed" %in% comp_types) {
      warning(
        "Mixture densities for mixed distributions can't ",
        "currently be computed. Values might be wrong!"
      )
    }

    discrete_parts <- comp_types %in% "discrete"
    discrete_mass <- rowSums(densmat[, discrete_parts, drop = FALSE]) > 0.0
    densmat[discrete_mass, !discrete_parts] <- 0.0

    res <- rowSums(probmat * densmat)
    if (log) res <- log(res)
    res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    params$probs <- lapply(params$probs, rep_len, length(q))
    probmat <- do.call(cbind, params$probs)
    probmat <- probmat / rowSums(probmat)

    cdfmat <- map_dbl_matrix(
      params$dists,
      ~.$dist$probability(
        q = q,
        lower.tail = lower.tail,
        with_params = .$params
      ),
      length(q)
    )

    res <- rowSums(probmat * cdfmat)
    if (log.p) res <- log(res)
    res
  },
  support = function(x, params) {
    rowSums(map_lgl_matrix(
      params$dists,
      ~.$dist$is_in_support(x, .$params),
      length(x)
    )) > 0L
  },
  get_components = function() {
    private$.default_params$dists
  },
  get_param_constraints = function() {
    ph <- self$get_placeholders()
    comps <- self$get_components()
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
        comp_constr <- comps[[i]]$get_param_constraints()
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
  get_type = function() {
    if (all(purrr::map_lgl(self$get_components(), ~.$is_discrete()))) {
      "discrete"
    } else if (all(purrr::map_lgl(self$get_components(), ~.$is_continuous()))) {
      "continuous"
    } else {
      "mixed"
    }
  },
  is_discrete = function(x, params) {
    if (self$get_type() == "mixed") {
      k_types <- purrr::map_chr(params$dists, ~.$dist$get_type())

      are_discrete <- rowSums(map_lgl_matrix(
        params$dists[k_types == "discrete"],
        ~.$dist$is_in_support(x, .$params),
        length(x)
      )) > 0.0

      if (all(k_types != "mixed")) {
        are_discrete <- are_discrete | rowSums(map_lgl_matrix(
          params$dists[k_types == "mixed"],
          ~.$dist$is_discrete_at(x, .$params),
          length(x)
        )) > 0.0
      }

      are_discrete
    } else {
      super$is_discrete_at(x, params)
    }
  },
  get_dof = function() {
    sdof <- super$get_dof()
    if (length(self$get_placeholders()$probs)) {
      sdof <- sdof - 1L
    }
    sdof
  },
  has_capability = function(caps) {
    super$has_capability(caps) &
      rowSums(
        1 - map_lgl_matrix(
          self$get_components(),
          function(comp) comp$has_capability(caps),
          length(caps)
        )
      ) == 0
  },
  tf_is_discrete_at = function() {
    if (self$is_continuous()) return(super$tf_is_discrete_at())

    comp_discretes <- lapply(self$get_components(), function(comp) {
      comp$tf_is_discrete_at()
    })

    function(x, args) {
      are_discrete <- tf$stack(
        lapply(
          seq_along(comp_discretes),
          function(i) comp_discretes[[i]](x, args[["dists"]][[i]])
        )
      )
      tf$math$reduce_any(are_discrete, axis = 1L)
    }
  },
  tf_logdensity = function() {
    comps <- self$get_components()
    comps_discrete <- lapply(comps, function(comp) comp$tf_is_discrete_at())
    comps_logdens <- lapply(comps, function(comp) comp$tf_logdensity())
    k <- length(comps)

    function(x, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      dist_args <- args[["dists"]]

      are_discrete <- tf$stack(lapply(seq_along(comps_discrete), function(i) {
        comps_discrete[[i]](x, dist_args[[i]])
      }), axis = 1L)

      # Component densities
      logdenss <- tf$stack(lapply(seq_along(comps_logdens), function(i) {
        comps_logdens[[i]](x, dist_args[[i]])
      }), axis = 1L)

      # Remove continuous densities where there is discrete support
      logdens_real <- tf$where(
        are_discrete == tf$reduce_any(are_discrete, axis = 1L, keepdims = TRUE),
        logdenss,
        K$neg_inf
      )

      logdens_safe <- tf$where(logdens_real > K$neg_inf, logdens_real, K$zero)

      tf$reduce_logsumexp(tf$where(logdens_real > K$neg_inf, log(probs) + logdens_safe, K$neg_inf), axis = 1L)
    }
  },
  tf_logprobability = function() {
    comps_logprob <- lapply(
      self$get_components(),
      function(comp) comp$tf_logprobability()
    )
    k <- length(comps_logprob)

    function(qmin, qmax, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))
      dist_args <- args[["dists"]]

      all_logprobs <- lapply(seq_along(comps_logprob), function(i) {
        comps_logprob[[i]](qmin, qmax, dist_args[[i]])
      })
      logprobs <- tf$stack(all_logprobs, axis = 1L)

      logprobs_safe <- tf$where(logprobs > K$neg_inf, logprobs, K$zero)

      tf$reduce_logsumexp(tf$where(logprobs > K$neg_inf, log(probs) + logprobs_safe, K$neg_inf), axis = 1L)
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
          out$probs <- outputs[[1L]]
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
  compile_sample = function() {
    comps <- self$get_components()
    k <- length(comps)
    ph_probs <- length(self$get_placeholders()$probs) > 0L

    comp_sample <- lapply(comps, function(dist) dist$compile_sample())
    comp_param_counts <- vapply(comp_sample, function(fun) attr(fun, "n_params"), integer(1L))
    comp_param_ends <- cumsum(comp_param_counts)
    comp_param_starts <- comp_param_ends - comp_param_counts + 1L

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L

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
    for (i in seq_along(comp_sample)) {
      comp_param_expr <- if (comp_param_counts[i] > 0L) {
        bquote(param_matrix[slot == .(i), .(comp_param_starts[i]):.(comp_param_ends[i]), drop = FALSE])
      } else {
        NULL
      }

      sampling_code[[i + 3L]] <- bquote(
        res[slot == .(i)] <- comp_sample[[.(i)]](
          num_samples[.(i)],
          .(comp_param_expr)
        )
      )
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
    ph_probs <- length(self$get_placeholders()$probs) > 0L

    comp_density <- lapply(comps, function(dist) dist$compile_density())
    comp_param_counts <- vapply(comp_density, function(fun) attr(fun, "n_params"), integer(1L))
    comp_param_ends <- cumsum(comp_param_counts)
    comp_param_starts <- comp_param_ends - comp_param_counts + 1L

    comp_types <- vapply(comps, function(comp) comp$get_type(), character(1L))
    # TODO implement mixed type too
    stopifnot(all(comp_types %in% c("discrete", "continuous")))

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L

    component_code <- bquote({
      compmat <- matrix(nrow = length(x), ncol = .(k))
    })

    for (i in seq_len(k)) {
      comp_param_expr <- if (comp_param_counts[i] > 0L) {
        bquote(param_matrix[slot == .(i), .(comp_param_starts[i]):.(comp_param_ends[i]), drop = FALSE])
      } else {
        NULL
      }

      component_code[[i + 2L]] <- bquote(
        compmat[, .(i)] <- comp_density[[.(i)]](x, .(comp_param_expr))
      )
    }

    if (self$get_type() == "mixed") {
      component_code[[k + 3L]] <- bquote({
        is_discrete <- matrixStats::rowAnys(compmat[, .(which(comp_types == "discrete"))])
        compmat[is_discrete, .(which(comp_types == "continuous"))] <- 0.0
      })
    }

    mixture_code <- if (ph_probs) {
      bquote({
        probmat <- param_matrix[, .(n_params - k + 1L):.(n_params), drop = FALSE]
        probmat <- probmat / rowSums(probmat)
        res <- rowSums(compmat * probmat)
      })
    } else {
      probs <- as.numeric(self$get_params()$probs)
      probs <- probs / sum(probs)
      bquote(res <- drop(compmat %*% .(probs)))
    }

    as_compiled_distribution_function(
      eval(bquote(function(x, param_matrix, log = FALSE) {
        .(component_code)
        .(mixture_code)
        if (log) log(res) else res
      })),
      n_params
    )
  },
  compile_probability = function() {
    comps <- self$get_components()
    k <- length(comps)
    ph_probs <- length(self$get_placeholders()$probs) > 0L

    comp_probability <- lapply(comps, function(dist) dist$compile_probability())
    comp_param_counts <- vapply(comp_probability, function(fun) attr(fun, "n_params"), integer(1L))
    comp_param_ends <- cumsum(comp_param_counts)
    comp_param_starts <- comp_param_ends - comp_param_counts + 1L

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L

    component_code <- bquote({
      compmat <- matrix(nrow = length(x), ncol = .(k))
    })

    for (i in seq_len(k)) {
      comp_param_expr <- if (comp_param_counts[i] > 0L) {
        bquote(param_matrix[slot == .(i), .(comp_param_starts[i]):.(comp_param_ends[i]), drop = FALSE])
      } else {
        NULL
      }

      component_code[[i + 2L]] <- bquote(
        compmat[, .(i)] <- comp_probability[[.(i)]](q, .(comp_param_expr), lower.tail = lower.tail)
      )
    }

    mixture_code <- if (ph_probs) {
      bquote({
        probmat <- param_matrix[, .(n_params - k + 1L):.(n_params), drop = FALSE]
        probmat <- probmat / rowSums(probmat)
        res <- rowSums(compmat * probmat)
      })
    } else {
      probs <- as.numeric(self$get_params()$probs)
      probs <- probs / sum(probs)
      bquote(res <- drop(compmat %*% .(probs)))
    }

    as_compiled_distribution_function(
      eval(bquote(function(q, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        .(component_code)
        .(mixture_code)
        if (log.p) log(res) else res
      }))
    )
  },
  compile_probability_interval = function() {
    comps <- self$get_components()
    k <- length(comps)
    ph_probs <- length(self$get_placeholders()$probs) > 0L

    comp_probability_interval <- lapply(comps, function(dist) dist$compile_probability_interval())
    comp_param_counts <- vapply(comp_probability_interval, function(fun) attr(fun, "n_params"), integer(1L))
    comp_param_ends <- cumsum(comp_param_counts)
    comp_param_starts <- comp_param_ends - comp_param_counts + 1L

    n_params <- sum(comp_param_counts) + if (ph_probs) k else 0L

    component_code <- bquote({
      compmat <- matrix(nrow = length(x), ncol = .(k))
    })

    for (i in seq_len(k)) {
      comp_param_expr <- if (comp_param_counts[i] > 0L) {
        bquote(param_matrix[slot == .(i), .(comp_param_starts[i]):.(comp_param_ends[i]), drop = FALSE])
      } else {
        NULL
      }

      component_code[[i + 2L]] <- bquote(
        compmat[, .(i)] <- comp_probability_interval[[.(i)]](qmin, qmax, .(comp_param_expr))
      )
    }

    mixture_code <- if (ph_probs) {
      bquote({
        probmat <- param_matrix[, .(n_params - k + 1L):.(n_params), drop = FALSE]
        probmat <- probmat / rowSums(probmat)
        res <- rowSums(compmat * probmat)
      })
    } else {
      probs <- as.numeric(self$get_params()$probs)
      probs <- probs / sum(probs)
      bquote(res <- drop(compmat %*% .(probs)))
    }

    as_compiled_distribution_function(
      eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
        .(component_code)
        .(mixture_code)
        if (log.p) log(res) else res
      }))
    )
  }
)

#' @rdname fit_dist_start
#'
#' @param dists_start List of initial parameters for all component
#' distributions. If left empty, initialisation will be automatically performed
#' using [fit_dist_start()] with all observations in the support of each
#' respective component.
#' @export
fit_dist_start.MixtureDistribution <- function(dist, obs, dists_start = NULL,
                                               ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs)
  n <- nrow(obs)
  res <- dist$get_placeholders()
  ph_dists <- lengths(res$dists) > 0L
  ph_probs <- length(res$probs) > 0L

  comp_dists <- dist$get_components()
  comp_types <- purrr::map_chr(comp_dists, ~.$get_type())
  supp_mat <- map_lgl_matrix(
    comp_dists,
    ~.$is_in_support(obs$x),
    nrow(obs)
  )
  k_discrete <- which(comp_types %in% "discrete")
  k_continuous <- which(comp_types %in% "continuous")

  i_cens <- is.na(obs$x)
  i_obs <- !i_cens
  i_discrete <- rowSums(supp_mat[, k_discrete, drop = FALSE]) > 0L

  supp_mat[i_discrete, k_continuous] <- FALSE
  supp_mat[i_cens, ] <- TRUE

  res$dists[!ph_dists] <- rep_len(list(list()), sum(!ph_dists))
  if (any(ph_dists)) {
    if (!is.null(dists_start)) {
      for (comp in which(ph_dists)) {
        res$dists[[comp]] <- dists_start[[comp]]
      }
    } else {
      for (comp in which(ph_dists)) {
        res$dists[[comp]] <- fit_dist_start(
          dist = comp_dists[[comp]],
          obs = obs[supp_mat[, comp], ],
          ...
        )
      }
    }
  }

  if (ph_probs) {
    # Compute partial component densities
    densmat <- map_dbl_matrix(
      seq_along(comp_dists),
      function(i) {
        dens <- numeric(n)
        dens[i_obs] <- comp_dists[[i]]$density(
          x[i_obs], with_params = res$dists[[i]]
        )
        dens[i_cens] <- (
          comp_dists[[i]]$probability(
            obs$xmax[i_cens], with_params = res$dists[[i]]
          ) - comp_dists[[i]]$probability(
            obs$xmin[i_cens], with_params = res$dists[[i]]
          )
        )

        if (!comp_dists[[i]]$is_continuous()) {
          is_disc <- comp_dists[[i]]$is_discrete_at(
            obs$xmin[i_cens], with_params = res$dists[[i]]
          )

          if (any(is_disc)) {
            dens[i_cens][is_disc] <- dens[i_cens][is_disc] + comp_dists[[i]]$density(
              obs$xmin[i_cens][is_disc], with_params = res$dists[[i]]
            )
          }
        }

        dens
      },
      n
    )

    if (dist$get_type() == "mixed") {
      is_disc <- map_lgl_matrix(
        seq_along(comp_dists),
        function(i) {
          comp_dists[[i]]$is_discrete_at(
            x, with_params = res$dists[[i]]
          )
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

# TODO fit_mixture assumes probs to be tunable
# If probs is fixed, we should instead just fit all dists with appropriately
# weighted observations
#' @include fit_mixture.R
#' @export
fit_dist.MixtureDistribution <- fit_mixture

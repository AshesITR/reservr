#' Fit a generic mixture using an ECME-Algorithm
#'
#' @param dist A `MixtureDistribution` specifying the structure of the mixture.
#' Free parameters are to be optimised. The dominating measure for likelihoods
#' must be constant, so for example [dist_dirac()] may not have its `point`
#' parameter free.
#' @param min_iter Minimum number of EM-Iterations
#' @param max_iter Maximum number of EM-Iterations (weight updates)
#' @param skip_first_e Skip the first E-Step (update Probability weights)?
#' This can help if the initial values cause a mixture component to vanish in
#' the first E-Step before the starting values can be improved.
#' @param tolerance Numerical tolerance.
#' @param trace Include tracing information in output?
#' If `TRUE`, additional tracing information will be added to the result list.
#' @param ...  Passed to [fit_dist_start()] if `start` is missing.
#' @inheritParams fit_dist
#'
#' @return A list with elements
#'
#'   * `params` the fitted parameters in the same structure as `init`.
#'   * `params_hist` (if `trace` is TRUE) the history of parameters
#'     (after each e- and m- step)
#'   * `iter` the number of outer EM-iterations
#'   * `logLik` the final log-likelihood
#'
#' @examples
#' dist <- dist_mixture(
#'   list(
#'     dist_dirac(0.0),
#'     dist_exponential()
#'   )
#' )
#'
#' params <- list(
#'   probs = list(0.1, 0.9),
#'   dists = list(
#'     list(),
#'     list(rate = 1.0)
#'   )
#' )
#'
#' x <- dist$sample(100L, with_params = params)
#'
#' fit_mixture(dist, x)
#'
#' @export
#' @family distribution fitting functions
fit_mixture <- function(dist, obs, start,
                        min_iter = 0L, max_iter = 100L, skip_first_e = FALSE,
                        tolerance = 1e-5, trace = FALSE, ...) {
  stopifnot(
    "`dist` must be a MixtureDistribution." =
      inherits(dist, "MixtureDistribution")
  )

  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start, ...)

  mix_comps <- dist$get_components()
  comp_types <- purrr::map_chr(mix_comps, ~.$get_type())

  n <- nrow(obs)
  k <- length(mix_comps)

  # Determine dominating measure, i.e. discrete points.
  comp_discrete_supp <- map2_lgl_matrix(
    mix_comps[comp_types != "continuous"],
    start$dists[comp_types != "continuous"],
    ~.x$is_discrete_at(obs$x, with_params = .y),
    n
  )
  i_cens <- is.na(obs$x)
  i_discrete <- rowSums(comp_discrete_supp) > 0.0
  i_continuous <- !i_cens & !i_discrete

  n_cens <- sum(i_cens)
  n_discrete <- sum(i_discrete)
  n_continuous <- sum(i_continuous)

  k_discrete <- which(comp_types %in% "discrete")
  k_continuous <- which(comp_types %in% "continuous")
  k_tunable <- which(lengths(start$dists) > 0L)

  dens_matrix <- function(params) {
    res <- matrix(0.0, nrow = n, ncol = k)
    res[i_cens, ] <- map2_dbl_matrix(
      mix_comps,
      params$dists,
      ~.x$probability(obs$xmax[i_cens], with_params = .y) -
        .x$probability(obs$xmin[i_cens], with_params = .y),
      n_cens
    )
    res[i_discrete, k_discrete] <- map2_dbl_matrix(
      mix_comps[k_discrete],
      params$dists[k_discrete],
      ~.x$density(obs$x[i_discrete], with_params = .y),
      n_discrete
    )
    res[i_continuous, k_continuous] <- map2_dbl_matrix(
      mix_comps[k_continuous],
      params$dists[k_continuous],
      ~.x$density(obs$x[i_continuous], with_params = .y),
      n_continuous
    )
    res
  }

  trunc_matrix <- function(params) {
    map2_dbl_matrix(
      mix_comps,
      params$dists,
      function(dist, params) {
        res <- dist$probability(obs$tmax, with_params = params) -
          dist$probability(obs$tmin, with_params = params)

        if (!dist$is_continuous()) {
          disc <- dist$is_discrete_at(obs$tmin, with_params = params)
          res[disc] <- res[disc] +
            dist$density(obs$tmin[disc], with_params = params)
        }

        res
      },
      n
    )
  }

  neg_loglik <- function(probs, dens_mat, trunc_mat) {
    list(
      objective = -sum(
        obs$w * (
          log(dens_mat %*% probs) -
            log(trunc_mat %*% probs)
        )
      ),
      gradient = -colSums(
        obs$w * (
          dens_mat / drop(dens_mat %*% probs) -
            trunc_mat / drop(trunc_mat %*% probs)
        )
      )
    )
  }

  params <- start

  iter <- 0L
  dens_mat <- dens_matrix(params)
  trunc_mat <- trunc_matrix(params)

  if (is.null(params$probs)) {
    params$probs <- local({
      p0 <- colSums((dens_mat / rowSums(dens_mat)) / trunc_mat)
      p0 / sum(p0)
    })
  } else {
    params$probs <- as.numeric(params$probs)
  }

  llik_new <- -neg_loglik(params$probs, dens_mat, trunc_mat)$objective
  llik_improved <- TRUE

  if (trace) {
    # Initialize history of parameters

    params_hist <- list(
      init = params
    )
    params_hist[[1L]]$probs <- as.list(params$probs)
  }

  while (iter < max_iter && (llik_improved || iter < min_iter)) {
    # Backup current params if llik worsens so we can restore
    params_old <- params

    # E-Step ----
    # 1. Update probs
    if (iter > 0L || !skip_first_e) {
      params$probs <- local({
        p0 <- params$probs

        nloptr::slsqp(
          x0 = p0,
          fn = neg_loglik,
          dens_mat = dens_mat,
          trunc_mat = trunc_mat,
          lower = rep_len(0.0, k),
          upper = rep_len(1.0, k),
          heq = function(x) sum(x) - 1.0,
          heqjac = function(x) rep_len(1.0, k)
        )$par
      })

      if (trace) {
        curr_step <- list(params)
        curr_step[[1L]]$probs <- as.list(params$probs)
        names(curr_step) <- sprintf("e_%d", iter + 1L)

        params_hist <- c(
          params_hist,
          curr_step
        )
      }
    }

    # 2. Compute posterior probabilities
    z <- dens_mat %*% diag(params$probs) / trunc_mat
    z[is.na(z)] <- 0.0 # 0 / 0 where dens_mat == trunc_mat == 0
    z <- z / rowSums(z)

    # M-Step ----
    # 3. Update distribution parameters
    for (comp in k_tunable) {
      comp_dist <- mix_comps[[comp]]

      i_keep <- if (comp %in% k_discrete) {
        (i_discrete | i_cens) & z[, comp] > 0.0
      } else {
        (i_continuous | i_cens) & z[, comp] > 0.0
      }

      obs_comp <- obs[i_keep, ]
      obs_comp$w <- obs_comp$w * z[i_keep, comp]

      params$dists[[comp]] <- fit_dist(
        dist = comp_dist,
        obs = obs_comp,
        params$dists[[comp]],
        ...
      )$params

      if (trace) {
        curr_step <- list(params)
        curr_step[[1L]]$probs <- as.list(params$probs)
        names(curr_step) <- sprintf("m_%d_%d", iter + 1L, comp)

        params_hist <- c(
          params_hist,
          curr_step
        )
      }
    }

    # 4. Update loglik, density- and truncation-matrix
    iter <- iter + 1L
    dens_mat <- dens_matrix(params)
    trunc_mat <- trunc_matrix(params)
    llik_old <- llik_new
    llik_new <- -neg_loglik(params$probs, dens_mat, trunc_mat)$objective
    llik_improved <- (llik_new - llik_old) > tolerance

    if (llik_new < llik_old && iter > min_iter) {
      # log-Likelihood worse than previous iteration:
      # restore params of previous iteration
      params <- params_old
      llik_new <- llik_old
    }
  }

  params$probs <- as.list(params$probs)

  if (trace) {
    params_hist <- c(params_hist, list(final = params))
    list(
      params = params,
      params_hist = params_hist,
      iter = iter,
      logLik = structure(
        llik_new, class = "logLik",
        df = dist$get_dof(),
        nobs = n
      )
    )
  } else {
    list(
      params = params,
      iter = iter,
      logLik = structure(
        llik_new, class = "logLik",
        df = dist$get_dof(),
        nobs = n
      )
    )
  }
}

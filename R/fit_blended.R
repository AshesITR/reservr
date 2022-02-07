#' Fit a Blended mixture using an ECME-Algorithm
#'
#' @param dist A `BlendedDistribution`. It is assumed, that `breaks` and
#' `bandwidths` are not a placeholder and that `weights` are to be estimated.
#' @param ... Passed to [fit_dist_start()] if `start` is missing.
#' @inheritParams fit_mixture
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
#' dist <- dist_blended(
#'    list(
#'      dist_exponential(),
#'      dist_genpareto()
#'    )
#'  )
#'
#'  params <- list(
#'    probs = list(0.9, 0.1),
#'    dists = list(
#'      list(rate = 2.0),
#'      list(u = 1.5, xi = 0.2, sigmau = 1.0)
#'    ),
#'    breaks = list(1.5),
#'    bandwidths = list(0.3)
#'  )
#'
#'  x <- dist$sample(100L, with_params = params)
#'
#'  dist$default_params$breaks <- params$breaks
#'  dist$default_params$bandwidths <- params$bandwidths
#'  fit_blended(dist, x)
#'
#' @export
#' @family distribution fitting functions
fit_blended <- function(dist, obs, start,
                        min_iter = 0L, max_iter = 100L, skip_first_e = FALSE,
                        tolerance = 1e-5, trace = FALSE, ...) {
  stopifnot(
    "`dist` must be a BlendedDistribution." =
      inherits(dist, "BlendedDistribution")
  )

  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start, ...)

  blend_comps <- dist$get_components()
  blend <- dist$get_params()

  n <- nrow(obs)
  k <- length(blend_comps)

  k_tunable <- which(lengths(start$dists) > 0L)

  comp_supp <- map_lgl_matrix(
    seq_len(k),
    function(i) {
      if (i == 1L) {
        obs$xmin < blend$breaks[[1L]] + blend$bandwidths[[1L]]
      } else if (i == k) {
        obs$xmax > blend$breaks[[k - 1L]] - blend$bandwidths[[k - 1L]]
      } else {
        obs$xmax > blend$breaks[[i - 1L]] - blend$bandwidths[[i - 1L]] &
          obs$xmin < blend$breaks[[i]] + blend$bandwidths[[i]]
      }
    },
    n
  )

  u_mat <- unlist(blend$breaks)
  eps_mat <- unlist(blend$bandwidths)

  obs_trans <- blended_transition(obs$x, u = u_mat, eps = eps_mat, .gradient = TRUE, .extend_na = TRUE)
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

  # TODO ggfs. 1 / (pmax - pmin) aus dens_matrix und trunc_matrix kürzen
  # density matrix
  dens_matrix <- function(params) {
    map_dbl_matrix(
      seq_len(k),
      function(i) {
        c_obs <- !is.na(obs_trans[, i])
        c_cens <- !c_obs & comp_supp[, i]

        if (i == k) {
          pmax <- 1.0
        } else {
          pmax <- blend_comps[[i]]$probability(
            blend$breaks[[i]], with_params = params$dists[[i]]
          )
        }

        if (i == 1L) {
          pmin <- 0.0
        } else {
          pmin <- blend_comps[[i]]$probability(
            blend$breaks[[i - 1L]], with_params = params$dists[[i]]
          )
          if (blend_comps[[i]]$is_discrete_at(blend$breaks[[i - 1L]], with_params = params$dists[[i]])) {
            pmin <- pmin + blend_comps[[i]]$density(
              blend$breaks[[i - 1L]], with_params = params$dists[[i]]
            )
          }
        }

        dens <- numeric(n)
        dens[c_obs] <- blend_comps[[i]]$density(
          obs_trans[c_obs, i], with_params = params$dists[[i]]
        ) * obs_trans_gradient[c_obs, i] / (pmax - pmin)
        dens[c_cens] <- (
          blend_comps[[i]]$probability(
            obs_trans_max[c_cens, i], with_params = params$dists[[i]]
          ) - blend_comps[[i]]$probability(
            obs_trans_min[c_cens, i], with_params = params$dists[[i]]
          )
        )
        if (!blend_comps[[i]]$is_continuous()) {
          disc <- blend_comps[[i]]$is_discrete_at(
            obs_trans_min[, i], with_params = params$dists[[i]]
          )
          dens[disc & c_cens] <- dens[disc & c_cens] + blend_comps[[i]]$density(
            obs_trans_min[disc & c_cens, i], with_params = params$dists[[i]]
          )
        }

        dens[c_cens] <- dens[c_cens] / (pmax - pmin)
        dens
      },
      n
    )
  }

  # truncation probability matrix
  trunc_matrix <- function(params) {
    map_dbl_matrix(
      seq_len(k),
      function(i) {
        if (i == k) {
          pmax <- 1.0
        } else {
          pmax <- blend_comps[[i]]$probability(
            blend$breaks[[i]], with_params = params$dists[[i]]
          )
        }

        if (i == 1L) {
          pmin <- 0.0
        } else {
          pmin <- blend_comps[[i]]$probability(
            blend$breaks[[i - 1L]], with_params = params$dists[[i]]
          )
          if (blend_comps[[i]]$is_discrete_at(blend$breaks[[i - 1L]], with_params = params$dists[[i]])) {
            pmin <- pmin - blend_comps[[i]]$density(
              blend$breaks[[i - 1L]], with_params = params$dists[[i]]
            )
          }
        }

        prob_tmax <- blend_comps[[i]]$probability(
          tmax_comp[, i], with_params = params$dists[[i]]
        )
        prob_tmin <- blend_comps[[i]]$probability(
          tmin_comp[, i], with_params = params$dists[[i]]
        )

        if (!blend_comps[[i]]$is_continuous()) {
          disc <- blend_comps[[i]]$is_discrete_at(
            tmin_comp[, i], with_params = params$dists[[i]]
          )
          prob_tmin[disc] <- prob_tmin[disc] - blend_comps[[i]]$density(
            tmin_comp[disc, i], with_params = params$dists[[i]]
          )
        }

        (prob_tmax - prob_tmin) / (pmax - pmin)
      }, n
    )
  }

  neg_loglik <- function(probs, dens_mat, trunc_mat) {
    list(
      objective = -sum(
        obs$w * (log(dens_mat %*% probs) - log(trunc_mat %*% probs))
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
      p0 <- colSums(
        obs$w *
          (dens_mat / rowSums(dens_mat)) /
          (trunc_mat / rowSums(trunc_mat))
      )
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
        curr_step[[1L]]$weights <- as.list(params$weights)
        names(curr_step) <- sprintf("e_%d", iter + 1L)

        params_hist <- c(
          params_hist,
          curr_step
        )
      }
    }

    # 2. Compute posterior probabilities
    z <- dens_mat %*% diag(params$probs)
    z <- z / rowSums(z)

    # M-Step ----
    # 3. Update distribution parameters
    for (comp in k_tunable) {
      comp_dist <- blend_comps[[comp]]
      i_keep <- comp_supp[, comp]

      obs_comp <- trunc_obs(
        x = obs_trans[i_keep, comp],
        xmin = obs_trans_min[i_keep, comp],
        xmax = obs_trans_max[i_keep, comp],
        tmin = tmin_comp[i_keep, comp],
        tmax = tmax_comp[i_keep, comp],
        w = obs$w[i_keep] * z[i_keep, comp]
      )

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

      obs_comp <- truncate_obs(obs_comp, comp_min, comp_max)

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

    # 4. Update loglik, density matrix and kappa weights
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

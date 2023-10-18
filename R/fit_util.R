.mixture_neg_loglik <- function(obs, probs, dens_mat, trunc_mat) {
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

.mixture_e_step <- function(obs, probs, dens_mat, trunc_mat) {
  k <- length(probs)

  nloptr::slsqp(
    x0 = probs,
    fn = .mixture_neg_loglik,
    obs = obs,
    dens_mat = dens_mat,
    trunc_mat = trunc_mat,
    lower = rep_len(0.0, k),
    upper = rep_len(1.0, k),
    heq = function(x) sum(x) - 1.0,
    heqjac = function(x) rep_len(1.0, k)
  )$par
}

.get_init_x <- function(obs, .min = -Inf, .max = Inf) {
  res <- obs$x
  i_cens <- is.na(obs$x) & obs$xmin > .min & obs$xmax < .max
  r_cens <- is.na(obs$x) & obs$xmax >= .max & is.finite(obs$xmin)
  l_cens <- is.na(obs$x) & obs$xmin <= .min & is.finite(obs$xmax)
  res[i_cens] <- 0.5 * obs$xmin[i_cens] + 0.5 * obs$xmax[i_cens]
  res[r_cens] <- obs$xmin[r_cens]
  res[l_cens] <- obs$xmax[l_cens]
  res
}

.assert_set <- function(ph, disallow, distname) {
  if (any(disallow %in% ph)) {
    not_allowed <- intersect(disallow, ph)
    msg <- paste0(
      enumerate_strings(not_allowed),
      " cannot be a placeholder for ",
      distname,
      " distribution fitting."
    )
    stop(msg)
  }
}

.fit_dist_objective <- function(dist, obs, i_obs, i_cens, param_names,
                                param_bounds = NULL) {
  if (all(dist$has_capability(c("diff_density", "diff_probability")))) {
    function(par) {
      names(par) <- param_names
      wp <- inflate_params(par)

      if (anyNA(par) ||
            any(par < param_bounds$lower) ||
            any(par > param_bounds$upper)) {
        return(NaN)
      }

      log_f_obs <- dist$density(obs$x[i_obs], log = TRUE, with_params = wp)
      P_cens <- dist$probability(obs$xmax[i_cens], with_params = wp) -
        dist$probability(obs$xmin[i_cens], with_params = wp)
      disc_xmin <- dist$is_discrete_at(obs$xmin[i_cens], with_params = wp)
      P_cens[disc_xmin] <- P_cens[disc_xmin] + dist$density(obs$xmin[i_cens][disc_xmin], with_params = wp)

      P_trunc <- dist$probability(obs$tmax, with_params = wp) -
        dist$probability(obs$tmin, with_params = wp)
      disc_tmin <- dist$is_discrete_at(obs$tmin, with_params = wp)
      P_trunc[disc_tmin] <- P_trunc[disc_tmin] + dist$density(obs$tmin[disc_tmin], with_params = wp)

      log_df_obs <- flatten_params_matrix(
        dist$diff_density(obs$x[i_obs], log = TRUE, with_params = wp)
      )
      log_dP_cens <- (
        flatten_params_matrix(
          dist$diff_probability(obs$xmax[i_cens], with_params = wp)
        ) - flatten_params_matrix(
          dist$diff_probability(obs$xmin[i_cens], with_params = wp)
        )
      )
      log_dP_cens[disc_xmin, ] <- log_dP_cens[disc_xmin, ] + flatten_params_matrix(
        dist$diff_density(obs$xmin[i_cens][disc_xmin], with_params = wp)
      )
      log_dP_cens <- log_dP_cens / P_cens
      log_dP_trunc <- (
        flatten_params_matrix(
          dist$diff_probability(obs$tmax, with_params = wp)
        ) - flatten_params_matrix(
          dist$diff_probability(obs$tmin, with_params = wp)
        )
      )
      log_dP_trunc[disc_tmin, ] <- log_dP_trunc[disc_tmin, ] + flatten_params_matrix(
        dist$diff_density(obs$tmin[disc_tmin], with_params = wp)
      )
      log_dP_trunc <- log_dP_trunc / P_trunc

      list(
        objective = -sum(obs$w[i_obs] * log_f_obs) -
          sum(obs$w[i_cens] * log(P_cens)) +
          sum(obs$w * log(P_trunc)),
        gradient = -colSums(obs$w[i_obs] * log_df_obs) -
          colSums(obs$w[i_cens] * log_dP_cens) +
          colSums(obs$w * log_dP_trunc)
      )
    }
  } else {
    function(par) {
      names(par) <- param_names
      wp <- inflate_params(par)

      if (anyNA(par) ||
            any(par < param_bounds$lower) ||
            any(par > param_bounds$upper)) {
        return(NaN)
      }

      P_cens <- dist$probability(obs$xmax[i_cens], with_params = wp) -
        dist$probability(obs$xmin[i_cens], with_params = wp)
      disc_xmin <- dist$is_discrete_at(obs$xmin[i_cens], with_params = wp)
      P_cens[disc_xmin] <- P_cens[disc_xmin] + dist$density(obs$xmin[i_cens][disc_xmin], with_params = wp)

      P_trunc <- dist$probability(obs$tmax, with_params = wp) -
        dist$probability(obs$tmin, with_params = wp)
      disc_tmin <- dist$is_discrete_at(obs$tmin, with_params = wp)
      P_trunc[disc_tmin] <- P_trunc[disc_tmin] + dist$density(obs$tmin[disc_tmin], with_params = wp)

      -sum(obs$w[i_obs] * dist$density(obs$x[i_obs], log = TRUE, with_params = wp)) -
        sum(obs$w[i_cens] * log(P_cens)) +
        sum(obs$w * log(P_trunc))
    }
  }
}

.check_fit_dist_start <- function(dist, obs, start, ...) {
  if (missing(start)) {
    start <- fit_dist_start(dist, obs, ...)
  }

  ph_nms <- names(flatten_params(dist$get_placeholders()))
  init_nms <- names(flatten_params(start))

  if (!setequal(init_nms, ph_nms)) {
    extra_nms <- setdiff(init_nms, ph_nms)
    missing_nms <- setdiff(ph_nms, init_nms)

    if (length(extra_nms)) {
      if (length(extra_nms) > 5L) {
        extra_nms <- c(paste0("'", extra_nms[1L:4L], "'"), "...")
      } else {
        extra_nms <- paste0("'", extra_nms, "'")
      }
      extra_nms <- paste0(
        "extra names: ",
        paste(extra_nms, collapse = ", ")
      )
    }

    if (length(missing_nms)) {
      if (length(missing_nms) > 5L) {
        missing_nms <- c(paste0("'", missing_nms[1L:4L], "'"), "...")
      } else {
        missing_nms <- paste0("'", missing_nms, "'")
      }
      missing_nms <- paste0(
        "missing names: ",
        paste(missing_nms, collapse = ", ")
      )
    }

    detail_msg <- paste(c(missing_nms, extra_nms), collapse = "; ")

    stop("`start` is inconsistent with `dist$get_placeholders()`: ", detail_msg)
  }

  start
}

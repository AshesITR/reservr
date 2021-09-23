#' Fit an Erlang mixture using an ECME-Algorithm
#'
#' @param dist An `ErlangMixtureDistribution`. It is assumed, that both `probs`
#' and `scale` are to be estimated.
#' @param parallel Enable experimental parallel evaluation of expected log-likelihood?
#' @param ... Passed to [fit_dist_start()] if `start` is missing.
#' @inheritParams fit_mixture
#'
#' @return A list with elements
#'
#'   * `params` the fitted parameters in the same structure as `init`.
#'   * `params_hist` (if `trace` is TRUE) the history of parameters
#'      (after each e- and m- step). Otherwise an empty list.
#'   * `iter` the number of outer EM-iterations
#'   * `logLik` the final log-likelihood
#'
#' @examples
#' dist <- dist_erlangmix(list(NULL, NULL, NULL))
#' params <- list(
#'   shapes = list(1L, 4L, 12L),
#'   scale = 2.0,
#'   probs = list(0.5, 0.3, 0.2)
#' )
#' x <- dist$sample(100L, with_params = params)
#' fit_erlang_mixture(dist, x, init = "kmeans")
#'
#' @export
#' @family distribution fitting functions
fit_erlang_mixture <- function(dist, obs, start,
                               min_iter = 0L, max_iter = 100L,
                               skip_first_e = FALSE, tolerance = 1e-5,
                               trace = FALSE, parallel = FALSE, ...) {
  stopifnot(
    "`dist` must be an ErlangMixtureDistribution." =
      inherits(dist, "ErlangMixtureDistribution")
  )

  obs <- as_trunc_obs(obs)
  params <- .check_fit_dist_start(dist, obs, start, ...)

  if ("shapes" %in% names(params)) {
    shapes <- as.numeric(params$shapes)
  } else {
    shapes <- as.numeric(dist$get_params()$shapes)
  }

  probs <- as.numeric(params$probs)
  scale <- as.numeric(params$scale)

  if ("shapes" %in% names(params)) {
    .erlang_shape_search(
      obs = obs,
      probs = probs,
      shapes = shapes,
      scale = scale,
      min_iter = min_iter,
      max_iter = max_iter,
      skip_first_e = skip_first_e,
      tolerance = tolerance,
      trace = trace,
      parallel = parallel
    )
  } else {
    .erlang_em(
      obs = obs,
      probs = probs,
      shapes = shapes,
      scale = scale,
      min_iter = min_iter,
      max_iter = max_iter,
      skip_first_e = skip_first_e,
      tolerance = tolerance,
      trace = trace,
      free_vars = names(params),
      parallel = parallel
    )
  }
}

.erlang_shape_search <- function(obs, probs, shapes, scale,
                                 min_iter = 0L, max_iter = 100L,
                                 skip_first_e = FALSE, tolerance = 1e-5,
                                 trace = FALSE, parallel = FALSE) {
  # TODO implement shape search tracing
  ord <- order(shapes)
  shapes <- shapes[ord]
  probs <- probs[ord]

  i_rev <- rev(seq_along(shapes))
  i_fwd <- seq_along(shapes)
  k <- length(shapes)

  em_fit <- function(shapes, params) {
    .erlang_em(
      obs = obs,
      probs = as.numeric(params$probs),
      shapes = shapes,
      scale = params$scale,
      min_iter = min_iter,
      max_iter = max_iter,
      skip_first_e = skip_first_e,
      tolerance = tolerance,
      trace = trace,
      parallel = parallel
    )
  }

  # Helper to move shapes
  move_at <- function(shapes, i, by) {
    shapes[i] <- shapes[i] + by
    shapes
  }

  curr_fit <- em_fit(shapes, list(probs = as.list(probs), scale = scale))
  fit_improved <- TRUE

  while (fit_improved) {
    curr_shapes <- as.numeric(curr_fit$params$shapes)
    fit_improved <- FALSE

    # 1. Bump up shapes
    for (i in i_rev) {
      shape_changed <- TRUE

      while (shape_changed) {
        shape_changed <- FALSE
        if (i == k || curr_shapes[i] < curr_shapes[i + 1L] - 1L) {
          new_fit <- em_fit(move_at(curr_shapes, i, by = 1L), curr_fit$params)
          if (new_fit$logLik > curr_fit$logLik + tolerance) {
            curr_fit <- new_fit
            curr_shapes <- as.numeric(curr_fit$params$shapes)
            shape_changed <- TRUE
            fit_improved <- TRUE
          }
        }
      }
    }

    # 2. Bump down shapes
    for (i in i_fwd) {
      shape_changed <- TRUE

      while (shape_changed) {
        shape_changed <- FALSE
        if (i == 1L && curr_shapes[1L] > 1L ||
          i > 1L && curr_shapes[i] > curr_shapes[i - 1L] + 1L) {
          new_fit <- em_fit(move_at(curr_shapes, i, by = -1L), curr_fit$params)
          if (new_fit$logLik > curr_fit$logLik + tolerance) {
            curr_fit <- new_fit
            curr_shapes <- as.numeric(curr_fit$params$shapes)
            shape_changed <- TRUE
            fit_improved <- TRUE
          }
        }
      }
    }
  }

  curr_fit
}

.erlang_em <- function(obs, probs, shapes, scale,
                       min_iter = 0L, max_iter = 100L, skip_first_e = FALSE,
                       tolerance = 1e-5, trace = FALSE,
                       free_vars = c("probs", "shapes", "scale"),
                       parallel = FALSE) {
  i_cens <- is.na(obs$x)
  i_obs <- !i_cens

  iter <- 0L
  dens_mat <- .erlang_dens_mat(obs, shapes, scale, i_obs, i_cens)
  trunc_mat <- .erlang_trunc_mat(obs, shapes, scale)
  llik_new <- -.mixture_neg_loglik(obs, probs, dens_mat, trunc_mat)$objective
  llik_improved <- TRUE

  # Initialize history of parameters
  make_params <- function(free_vars) {
    params <- vector("list", length(free_vars))
    names(params) <- free_vars

    if ("scale" %in% free_vars) params$scale <- scale
    if ("probs" %in% free_vars) params$probs <- as.list(probs)
    if ("shapes" %in% free_vars) params$shapes <- as.list(shapes)

    params
  }

  params_hist <- list()
  if (trace) {

    param_snap <- function(...) {
      curr_step <- list(make_params(free_vars))
      names(curr_step) <- paste(..., sep = "_")
      params_hist <<- c(params_hist, curr_step)
      invisible()
    }

  } else {

    param_snap <- function(...) {
      invisible()
    }

  }

  param_snap("init")

  while (iter < max_iter && (llik_improved || iter < min_iter)) {
    # Backup current params if llik worsens so we can restore
    probs_old <- probs
    scale_old <- scale

    # E-Step ----
    # 1. Update probs
    if (iter > 0L || !skip_first_e) {
      probs <- .mixture_e_step(obs, probs, dens_mat, trunc_mat)
      param_snap("e", iter + 1L)
    }

    # 2. Compute posterior probabilities
    z <- dens_mat %*% diag(probs)
    z <- z / rowSums(z)

    # M-Step ----
    # 3. Update scale parameter
    scale <- nloptr::tnewton(
      x0 = scale,
      fn = .erlang_neg_ellik,
      obs = obs,
      shapes = shapes,
      zadj = z * obs$w,
      lower = 0.0,
      upper = Inf,
      parallel = parallel
    )$par

    param_snap("m", iter + 1L)

    # 4. Update loglik, density- and truncation-matrix
    iter <- iter + 1L
    dens_mat <- .erlang_dens_mat(obs, shapes, scale, i_obs, i_cens)
    trunc_mat <- .erlang_trunc_mat(obs, shapes, scale)
    llik_old <- llik_new
    llik_new <- -.mixture_neg_loglik(obs, probs, dens_mat, trunc_mat)$objective
    llik_improved <- (llik_new - llik_old) > tolerance

    if (llik_new < llik_old && iter > min_iter) {
      # log-Likelihood worse than previous iteration:
      # restore params of previous iteration
      probs <- probs_old
      scale <- scale_old
      llik_new <- llik_old
    }
  }

  param_snap("final")

  list(
    params = make_params(free_vars),
    params_hist = params_hist,
    iter = iter,
    logLik = structure(
      llik_new, class = "logLik",
      df = 2L * length(shapes),
      #> dof: k (shapes) + k (probs) - 1 (sum(probs) == 1.0) + 1 (scale)
      nobs = nrow(obs)
    )
  )
}

.erlang_dens_mat <- function(obs, shapes, scale, i_obs, i_cens) {
  # Prevent warnings from dgamma()
  # scale = 0 can happen if nloptr tries the boundary.
  if (is.na(scale) || scale == 0.0) return(matrix(NaN, nrow = nrow(obs), ncol = length(shapes)))

  res <- matrix(0.0, nrow = nrow(obs), ncol = length(shapes))
  res[i_obs, ] <- dgamma_matrix(obs$x[i_obs], shapes, scale)
  res[i_cens, ] <- pgamma_diff_matrix(obs$xmin[i_cens], obs$xmax[i_cens], shapes, scale)
  res
}

.erlang_trunc_mat <- function(obs, shapes, scale) {
  # Prevent warnings from pgamma()
  # scale = 0 can happen if nloptr tries the boundary.
  if (is.na(scale) || scale == 0.0) return(matrix(NaN, nrow = nrow(obs), ncol = length(shapes)))
  pgamma_diff_matrix(obs$tmin, obs$tmax, shapes, scale)
}

.erlang_neg_ellik <- function(scale, obs, shapes, zadj, parallel) {
  if (is.na(scale) || scale == 0.0) return(list(objective = NaN, gradient = c("scale" = NaN)))

  trunc_erlangmix_ellik(
    xmin = obs$xmin, xmax = obs$xmax, tmin = obs$tmin, tmax = obs$tmax, weight = obs$w,
    shapes = shapes, scale = scale, zadj = zadj, parallel = parallel
  )
}

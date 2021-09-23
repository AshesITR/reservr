#' Quantiles of Distributions
#'
#' Produces quantiles corresponding to the given probabilities with configurable distribution parameters.
#'
#' @param x A `Distribution`.
#' @param probs Quantiles to compute.
#' @param with_params Optional list of distribution parameters. Note that if `x$has_capability("quantile")` is false,
#' `with_params` is assumed to contain only one set of parameters.
#' @param ... ignored
#' @param .start Starting value if quantiles are computed numerically. Must be within the support of `x`.
#'
#' @details
#' If `x$has_capability("quantile")` is true, this returns the same as `x$quantile(probs, with_params = with_params)`.
#' In this case, `with_params` may contain separate sets of parameters for each quantile to be determined.
#'
#' Otherwise, a numerical estimation of the quantiles is done using the density and probability function.
#' This method assumes `with_params` to cantain only one set of parameters.
#' The strategy uses two steps:
#'
#' 1. Find the smallest and largest quantiles in `probs` using a newton method starting from `.start`.
#' 2. Find the remaining quantiles with bisection using [stats::uniroot()].
#'
#' @return The quantiles of `x` corresponding to `probs` with parameters `with_params`.
#'
#' @examples
#' # With quantiles available
#' dist <- dist_normal(sd = 1)
#' qqs <- quantile(dist, probs = rep(0.5, 3), with_params = list(mean = 1:3))
#' stopifnot(all.equal(qqs, 1:3))
#'
#' # Without quantiles available
#' dist <- dist_erlangmix(shapes = list(1, 2, 3), scale = 1.0)
#' my_probs <- c(0, 0.01, 0.25, 0.5, 0.75, 1)
#' qqs <- quantile(
#'   dist, probs = my_probs,
#'   with_params = list(probs = list(0.5, 0.3, 0.2)), .start = 2
#' )
#'
#' all.equal(dist$probability(qqs, with_params = list(probs = list(0.5, 0.3, 0.2))), my_probs)
#' # Careful: Numerical estimation of extreme quantiles can result in out-of-bounds values.
#' # The correct 0-quantile would be 0 in this case, but it was estimated < 0.
#' qqs[1L]
#'
#' @export
quantile.Distribution <- function(x, probs = seq(0, 1, 0.25),
                                  with_params = list(), ..., .start = 0.0) {
  assert_that(
    is.Distribution(x),
    msg = "`x` must be a Distribution."
  )
  assert_that(
    is.numeric(probs),
    all(probs >= 0),
    all(probs <= 1),
    msg = "`probs` must be a numeric vector with elements in [0, 1]."
  )

  if (x$has_capability("quantile")) {
    x$quantile(probs, with_params = with_params)
  } else {
    # Numeric quantile approximation.
    # Idea:
    # 1. Start at init, find q0_min and q0_max, the finite quantiles of the
    #    smallest and largest probability using nloptr::tnewton (probs 0 and
    #    1 might be infinite)
    # 2. All other quantiles are either in [q0_min, init] or in [init, q0_max].
    #    Find these remaining quantiles using uniroot

    x$require_capability(
      c("probability", "density"),
      fun_name = "quantile.Distribution() without capability 'quantile'"
    )
    assert_that(
      is_scalar_double(.start),
      x$is_in_support(.start, with_params = with_params),
      msg = "`start` must be a single numeric within the support of `x`."
    )

    p0 <- x$probability(.start, with_params = with_params)
    out <- rep_len(NA_real_, length(probs))

    out[probs == p0] <- .start

    if (1.0 %in% probs && x$is_in_support(.Machine$double.xmax)) {
      out[probs == 1.0] <- Inf
    }

    if (0.0 %in% probs && x$is_in_support(-.Machine$double.xmax)) {
      out[probs == 0.0] <- -Inf
    }

    d_squared <- function(q, prob_target) {
      delta_prob <- x$probability(q, with_params = with_params) - prob_target
      list(
        objective = delta_prob^2,
        gradient = 2 * delta_prob * x$density(q, with_params = with_params)
      )
    }

    diff_prob <- function(q, prob_target) {
      x$probability(q, with_params = with_params) - prob_target
    }

    if (any(is.na(out) & probs < p0)) {
      p0_min <- min(probs[is.na(out) & probs < p0])
      q0_min <- nloptr::tnewton(.start, fn = d_squared, prob_target = p0_min)$par
      if (diff_prob(q0_min, p0_min) > 0.0) {
        out[probs == p0_min] <- q0_min
      }

      for (prob in sort(unique(probs[is.na(out) & probs < p0]))) {
        q0_min <- uniroot(
          diff_prob, lower = q0_min, upper = .start, prob_target = prob
        )$root
        out[probs == prob] <- q0_min
      }
    }

    if (any(is.na(out) & probs > p0)) {
      p0_max <- max(probs[is.na(out) & probs > p0])
      q0_max <- nloptr::tnewton(.start, fn = d_squared, prob_target = p0_max)$par
      if (diff_prob(q0_max, p0_max) < 0.0) {
        out[probs == p0_max] <- q0_max
      }

      remaining_probs <- sort(
        unique(probs[is.na(out) & probs > p0]),
        decreasing = TRUE
      )
      for (prob in remaining_probs) {
        q0_max <- uniroot(
          diff_prob, lower = .start, upper = q0_max, prob_target = prob
        )$root
        out[probs == prob] <- q0_max
      }
    }

    out
  }
}

#' @rdname fit_dist
#' @param object same as parameter `dist`
#' @export
fit.Distribution <- function(object, obs, start, ...) {
  fit_dist(dist = object, obs = obs, start = start, ...)
}

#' @export
fit_dist.Distribution <- function(dist, obs, start, ...) {
  obs <- as_trunc_obs(obs)
  start <- .check_fit_dist_start(dist, obs, start, ...)

  n <- nrow(obs)
  init_flat <- flatten_params(start)
  param_names <- names(flatten_params(dist$get_placeholders()))

  # Ensure consistent ordering.
  init_flat <- init_flat[param_names]

  i_cens <- is.na(obs$x)
  i_obs <- !i_cens

  neg_loglik <- .fit_dist_objective(
    dist, obs, i_obs, i_cens, param_names
  )

  bounds <- flatten_bounds(dist$get_param_bounds())
  constraint <- dist$get_param_constraints()

  if (!is.null(constraint)) {
    constraint_impl <- constraint

    constraint <- function(par) {
      names(par) <- param_names
      wp <- inflate_params(par)
      constraint_impl(wp)
    }

  }

  if (length(init_flat) == 1L && is.null(constraint)) {
    xopt <- muffle_nans_produced(nloptr::tnewton(
      x0 = init_flat,
      fn = neg_loglik,
      lower = bounds$lower,
      upper = bounds$upper
    ))
  } else if (is.null(constraint)) {
    xopt <- muffle_nans_produced(nloptr::lbfgs(
      x0 = init_flat,
      fn = neg_loglik,
      lower = bounds$lower,
      upper = bounds$upper
    ))
  } else {
    xopt <- muffle_nans_produced(nloptr::slsqp(
      x0 = init_flat,
      fn = neg_loglik,
      lower = bounds$lower,
      upper = bounds$upper,
      heq = constraint
    ))
  }
  names(xopt$par) <- param_names

  list(
    params = inflate_params(xopt$par),
    opt = xopt,
    logLik = structure(
      -xopt$value, class = "logLik",
      df = dist$get_dof(),
      nobs = n
    )
  )
}

#' @export
format.Distribution <- function(x, short = FALSE, ...) {
  dists <- x$get_components()
  if (length(dists)) {
    dists <- vapply(dists, format, character(1L), short = TRUE)
    prefix <- gsub("Distribution$", "", class(x)[1L])
    dist <- paste0(prefix, "<", paste(dists, collapse = ", "), ">")
  } else {
    dist <- class(x)[1L]
  }
  if (short) return(dist)
  a_an <- if (substr(dist, 1L, 1L) %in% c("A", "E", "I", "O", "U")) "An" else "A"
  paste0(a_an, " ", dist, " with ", x$get_dof(), " dof")
}

#' @export
format.ErlangMixtureDistribution <- function(x, short = FALSE, ...) {
  free_shapes <- "shapes" %in% names(x$get_placeholders())
  spec <- if (free_shapes) {
    paste0("k = ", length(x$get_components()))
  } else {
    paste(unlist(x$default_params$shapes), collapse = ", ")
  }
  dist <- paste0("ErlangMixture<", spec, ">")
  if (short) return(dist)
  paste0("Am ", dist, " with ", x$get_dof(), " dof")
}

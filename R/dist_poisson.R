#' Poisson Distribution
#'
#' See [stats::Poisson]
#'
#' The parameter can be overridden with
#' `with_params = list(lambda = ...)`.
#'
#' @param lambda Scalar rate parameter, or `NULL` as a placeholder.
#'
#' @return A `PoissonDistribution` object.
#' @export
#'
#' @examples
#' d_pois <- dist_poisson(lambda = 5.0)
#' x <- d_pois$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_pois,
#'   estimated = d_pois,
#'   with_params = list(
#'     estimated = inflate_params(
#'       fitdistrplus::fitdist(x, distr = "pois")$estimate
#'     )
#'   ),
#'   .x = 0:max(x)
#' )
#'
#' @family Distributions
dist_poisson <- function(lambda = NULL) {
  PoissonDistribution$new(lambda = lambda)
}

PoissonDistribution <- distribution_class_simple(
  name = "Poisson",
  fun_name = "pois",
  type = "discrete",
  params = list(lambda = interval(0, Inf)),
  support = I_NATURALS,
  diff_density = function(x, vars, log, params) {
    res <- vars

    if ("lambda" %in% names(vars)) {
      res$lambda <- if (log) {
        x / params$lambda - 1.0
      } else {
        (x / params$lambda - 1.0) * dpois(x, params$lambda)
      }
      supp <- self$is_in_support(x)
      res$lambda[!supp] <- if (log) NaN else 0.0
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if ("lambda" %in% names(vars)) {
      # Compute numeric derivative
      res$lambda <- gradient(
        func = function(lambda) {
          sum(ppois(
            q = q, lambda = lambda, lower.tail = lower.tail, log.p = log.p
          ))
        },
        x = params$lambda
      )
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    lambda <- tf$broadcast_to(args[["lambda"]], x$shape)

    ok <- tf_is_integerish(x) & x >= K$zero
    x_safe <- tf$where(ok, x, K$zero)

    tf$where(
      ok,
      x_safe * log(lambda) - lambda - tf$math$lgamma(x_safe + K$one),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    lambda <- tf$broadcast_to(args[["lambda"]], qmin$shape)

    lambda_safe <- tf$where(lambda == K$zero, K$one, lambda)

    qmin0 <- qmin <= K$zero
    qmax_ok <- qmax >= K$zero & tf$math$is_finite(qmax)
    qmax_nok <- tf$where(qmax < 0, K$neg_inf, K$zero)

    qmin_safe <- tf$math$maximum(K$zero, tf$math$ceil(qmin)) - K$one
    qmax_safe <- tf$math$maximum(K$zero, tf$math$floor(qmax))

    tf$where(
      lambda == K$zero,
      tf$where(
        qmin0,
        tf$where(qmax_ok, K$zero, qmax_nok),
        K$neg_inf
      ),
      tf$where(
        qmin0,
        tf$where(
          qmax_ok,
          log(tf$math$igammac(K$one + qmax_safe, lambda_safe)),
          qmax_nok
        ),
        tf$where(
          qmax_ok,
          log(tf$math$igammac(K$one + qmax_safe, lambda_safe) - tf$math$igammac(K$one + qmin_safe, lambda_safe)),
          log(tf$math$igamma(K$one + qmin_safe, lambda_safe))
        )
      )
    )
  },
  tf_is_discrete_at = function() function(x, args) { # nolint: brace.
    lambda <- tf$broadcast_to(args[["lambda"]], x$shape)
    tf$where(
      lambda == K$zero,
      x == K$zero,
      tf_is_integerish(x) & x >= K$zero
    )
  }
)

#' @export
fit_dist_start.PoissonDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0L)
  res <- dist$get_placeholders()
  ph <- names(res)
  if ("lambda" %in% ph) {
    res$lambda <- weighted_moments(x, obs$w, n = 1L)
  }
  res
}

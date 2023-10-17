#' Exponential distribution
#'
#' See [stats::Exponential].
#'
#' The parameter can be overridden with `with_params = list(rate = ...)`.
#'
#' @param rate Scalar rate parameter, or `NULL` as a placeholder.
#'
#' @return An `ExponentialDistribution` object.
#' @export
#'
#' @examples
#' rate <- 1
#' d_exp <- dist_exponential()
#' x <- d_exp$sample(20, with_params = list(rate = rate))
#' d_emp <- dist_empirical(x, positive = TRUE)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_exp,
#'   estimated = d_exp,
#'   with_params = list(
#'     theoretical = list(rate = rate),
#'     estimated = list(rate = 1 / mean(x))
#'   ),
#'   .x = seq(1e-4, 5, length.out = 100)
#' )
#'
#' @family Distributions
dist_exponential <- function(rate = NULL) {
  ExponentialDistribution$new(rate = rate)
}

ExponentialDistribution <- distribution_class_simple(
  name = "Exponential",
  fun_name = "exp",
  params = list(rate = I_POSITIVE_REALS),
  support = I_POSITIVE_REALS,
  diff_density = function(x, vars, log, params) {
    res <- vars
    if ("rate" %in% names(vars)) {
      res$rate <- if (log) {
        1 / params$rate - x
      } else {
        (1 / params$rate - x) * dexp(x, params$rate)
      }
      res$rate[is.infinite(x)] <- 0.0
    }
    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars
    if ("rate" %in% names(vars)) {
      res$rate <- if (lower.tail) {
        if (log.p) {
          q * pexp(q, params$rate, lower.tail = FALSE) / pexp(q, params$rate)
        } else {
          q * pexp(q, params$rate, lower.tail = FALSE)
        }
      } else {
        if (log.p) {
          -q
        } else {
          -q * pexp(q, params$rate, lower.tail = FALSE)
        }
      }
      res$rate[is.infinite(q)] <- 0.0
    }
    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    tf <- tensorflow::tf
    rate <- tf$squeeze(args[["rate"]])

    x_ok <- tf$math$is_finite(x) & x >= 0
    x_safe <- tf$where(x_ok, x, K$zero)
    tf$where(x_ok, log(rate) - x_safe * rate, K$neg_inf)
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    tf <- tensorflow::tf
    rate <- tf$squeeze(args[["rate"]])

    qmax_safe <- tf$where(tf$math$is_finite(qmax), qmax, K$zero)
    logprob_low <- tf$where(qmin > 0, -rate * tf$where(qmin > 0, qmin, K$one), K$zero)
    logprob_high <- tf$where(qmax_safe > 0, -rate * tf$where(qmax_safe > 0, qmax_safe, K$one), K$zero)

    tf$where(
      qmax > 0 & qmin < K$inf,
      tf$where(
        tf$math$is_finite(qmax),
        log(exp(logprob_low) - exp(logprob_high)),
        logprob_low
      ),
      K$neg_inf
    )
  }
)

#' @export
fit_dist_start.ExponentialDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0.0)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 1L)
  if ("rate" %in% ph) {
    res$rate <- 1.0 / mom[1L]
  }
  res
}

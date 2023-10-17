#' Binomial Distribution
#'
#' See [stats::Binomial]
#'
#' Both parameters can be overridden with
#' `with_params = list(size = ..., prob = ...)`.
#'
#' @param size Number of trials parameter (integer), or `NULL` as a placeholder.
#' @param prob Success probability parameter, or `NULL` as a placeholder.
#'
#' @return A `BinomialDistribution` object.
#' @export
#'
#' @examples
#' d_binom <- dist_binomial(size = 10, prob = 0.5)
#' x <- d_binom$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_binom,
#'   estimated = d_binom,
#'   with_params = list(
#'     estimated = list(
#'       size = max(x),
#'       prob = mean(x) / max(x)
#'     )
#'   ),
#'   .x = 0:max(x)
#' )
#'
#' @family Distributions
dist_binomial <- function(size = NULL, prob = NULL) {
  BinomialDistribution$new(size = size, prob = prob)
}

BinomialDistribution <- distribution_class_simple(
  name = "Binomial",
  fun_name = "binom",
  type = "discrete",
  params = list(
    size = interval(1L, Inf, include_lowest = TRUE, integer = TRUE),
    prob = interval(0.0, 1.0, closed = TRUE)
  ),
  support = function(x, params) {
    x == trunc(x) & x >= 0L & x <= params$size
  },
  diff_density = function(x, vars, log, params) {
    res <- vars
    dens <- dbinom(x = x, size = params$size, prob = params$prob)

    if ("size" %in% names(vars)) {
      dens_n_plus_1 <- dbinom(
        x = x, size = params$size + 1L, prob = params$prob, log = log
      )

      res$size <- if (log) {
        dens_n_plus_1 / dens - 1.0
      } else {
        dens_n_plus_1 - dens
      }
    }

    if ("prob" %in% names(vars)) {
      log_diff <- x / params$prob - (params$size - x) / (1 - params$prob)

      res$prob <- if (log) {
        log_diff[dens == 0.0] <- NaN
        log_diff
      } else {
        log_diff * dens
      }
    }

    res
  },
  diff_probability = function(q, vars, lower.tail, log.p, params) {
    res <- vars

    if ("size" %in% names(vars)) {
      prob_n_plus_1 <- pbinom(q = q, size = params$size, prob = params$prob)
      prob <- pbinom(q = q, size = params$size, prob = params$prob)

      res$size <- if (log.p) {
        prob_n_plus_1 / prob - 1.0
      } else {
        prob_n_plus_1 - prob
      }

      if (!lower.tail) res$size <- -res$size
    }

    if ("prob" %in% names(vars)) {
      res$prob <- if (log.p) {
        -params$size *
          dbinom(x = q, size = params$size - 1L, prob = params$prob) /
          pbinom(q = q, size = params$size, prob = params$prob,
                 lower.tail = lower.tail)
      } else {
        -params$size * dbinom(
          x = q, size = params$size - 1L, prob = params$prob
        )
      }

      if (!lower.tail) res$prob <- -res$prob
    }

    res
  },
  tf_logdensity = function() function(x, args) { # nolint: brace.
    size <- args[["size"]]
    prob <- args[["prob"]]

    ok <- tf_is_integerish(x) & x >= K$zero & x <= size
    x_safe <- tf$where(ok, x, K$zero)

    beta_x <- tf$math$lbeta(tf$stack(list(K$one + x_safe, K$one + size - x_safe), axis = 1L))

    tf$where(
      ok,
      -beta_x - log1p(size) + x_safe * log(prob) + (size - x_safe) * log1p(-prob),
      K$neg_inf
    )
  },
  tf_logprobability = function() function(qmin, qmax, args) { # nolint: brace.
    size <- args[["size"]]
    prob <- args[["prob"]]

    qmin0 <- qmin <= K$zero
    qmax_ok <- qmax >= K$zero & qmax < size
    qmin_safe <- tf$math$maximum(K$zero, tf$math$minimum(size, tf$math$ceil(qmin))) - K$one
    qmax_safe <- tf$math$maximum(K$zero, tf$math$minimum(size, tf$math$floor(qmax)))

    qmax_nok <- tf$where(qmax < K$zero, K$neg_inf, K$zero)

    tf$where(
      qmin0,
      tf$where(
        qmax_ok,
        log(tf$math$betainc(size - qmax_safe, qmax_safe + K$one, K$one - prob)),
        qmax_nok
      ),
      tf$where(
        qmax_ok,
        log(tf$math$betainc(size - qmax_safe, qmax_safe + K$one, K$one - prob) -
              tf$math$betainc(size - qmin_safe, qmin_safe + K$one, K$one - prob)),
        log1p(-tf$math$betainc(size - qmin_safe, qmin_safe + K$one, K$one - prob))
      )
    )
  },
  tf_is_discrete_at = function() function(x, args) { # nolint: brace.
    size <- args[["size"]]

    tf_is_integerish(x) & x >= K$zero & x <= size
  }
)

#' @export
fit_dist_start.BinomialDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  x <- .get_init_x(obs, .min = 0L)
  res <- dist$get_placeholders()
  ph <- names(res)
  mom <- weighted_moments(x, obs$w, n = 1L)

  if ("size" %in% ph) {
    size <- max(obs$xmax)
    res$size <- size
  } else {
    size <- dist$get_params()$size
  }

  if ("prob" %in% ph) {
    res$prob <- mom[1L] / size
  }
  res
}

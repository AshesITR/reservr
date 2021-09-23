#' Dirac (degenerate point) Distribution
#'
#' A degenerate distribution with all mass at a single point.
#'
#' The parameter can be overridden with
#' `with_params = list(point = ...)`.
#'
#' @param point The point with probability mass 1.
#'
#' @return A `DiracDistribution` object.
#' @export
#'
#' @examples
#' d_dirac <- dist_dirac(1.5)
#' d_dirac$sample(2L)
#' d_dirac$sample(2L, list(point = 42.0))
#'
#' @family Distributions
dist_dirac <- function(point = NULL) {
  DiracDistribution$new(point = point)
}

DiracDistribution <- distribution_class(
  name = "Dirac",
  type = "discrete",
  params = list(point = I_REALS),
  sample = function(n, params) {
    rep_len(params$point, n)
  },
  density = function(x, log = FALSE, params) {
    res <- as.numeric(x == params$point)
    if (log) {
      res[res == 0.0] <- -Inf
      res[res == 1.0] <- 0.0
    }
    res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    res <- if (lower.tail) params$point <= q else params$point > q
    res <- as.numeric(res)
    if (log.p) {
      res[res == 0.0] <- -Inf
      res[res == 1.0] <- 0.0
    }
    res
  },
  quantile = function(p, lower.tail = TRUE, log.p = FALSE, params) {
    params$point
  },
  support = function(x, params) {
    x %in% params$point
  },
  tf_is_discrete_at = function() function(x, args) {
    point <- tensorflow::tf$squeeze(args[["point"]])
    tensorflow::tf$equal(x, point)
  },
  tf_logdensity = function() function(x, args) {
    point <- tensorflow::tf$squeeze(args[["point"]])
    tensorflow::tf$where(x == point, K$zero, K$neg_inf)
  },
  tf_logprobability = function() function(qmin, qmax, args) {
    point <- tensorflow::tf$squeeze(args[["point"]])
    tensorflow::tf$where(qmin > point | qmax < point, K$neg_inf, K$zero)
  }
)

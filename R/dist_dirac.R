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
  },
  compile_sample = function() {
    if ("point" %in% names(self$get_placeholders())) {
      as_compiled_distribution_function(function(n, param_matrix) {
        param_matrix[, 1L]
      }, 1L)
    } else {
      as_compiled_distribution_function(eval(bquote(function(n, param_matrix) {
        rep_len(.(self$default_params$point), n)
      })), 0L)
    }
  },
  compile_density = function() {
    ph <- "point" %in% names(self$get_placeholders())
    as_compiled_distribution_function(
      eval(substitute(
        function(x, param_matrix, log = FALSE) {
          res <- as.numeric(x == point_expr)
          if (log) {
            res[res == 0.0] <- -Inf
            res[res == 1.0] <- 0.0
          }
          res
        },
        list(point_expr = if (ph) quote(param_matrix[, 1L]) else self$default_params$point)
      )),
      n_params = ph
    )
  },
  compile_probability = function() {
    ph <- "point" %in% names(self$get_placeholders())
    as_compiled_distribution_function(
      eval(substitute(
        function(q, param_matrix, lower.tail = TRUE, log.p = FALSE) {
          res <- if (lower.tail) point_expr <= q else point_expr > q
          res <- as.numeric(res)
          if (log.p) {
            res[res == 0.0] <- -Inf
            res[res == 1.0] <- 0.0
          }
          res
        },
        list(point_expr = if (ph) quote(param_matrix[, 1L]) else self$default_params$point)
      )),
      n_params = ph
    )
  },
  compile_quantile = function() {
    if ("point" %in% names(self$get_placeholders())) {
      as_compiled_distribution_function(function(p, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        param_matrix[, 1L]
      }, 1L)
    } else {
      as_compiled_distribution_function(eval(bquote(function(p, param_matrix, lower.tail = TRUE, log.p = FALSE) {
        rep_len(.(self$default_params$point), length(p))
      })), 0L)
    }
  }
)

#' @rdname Pareto
#' @name Pareto
#' @title The Pareto Distribution
#'
#' @description
#' These functions provide information about the Pareto distribution.
#' `dpareto` gives the density, `ppareto` gives the distribution
#' function, `qpareto` gives the quantile function and `rpareto` generates random
#' deviates.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n integer number of observations.
#' @param shape shape parameter (must be positive).
#' @param scale scale parameter (must be positive).
#' @param log,log.p logical; if `TRUE`, probabilities/densities
#' `p` are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#' \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#'
#' @details
#' If `shape` or `scale` are not specified, they assume the default values of `1`.
#'
#' The Pareto distribution with scale \eqn{\theta} and shape \eqn{\xi} has density
#'
#' \deqn{f(x) = \xi \theta^\xi / (x + \theta)^(\xi + 1)}
#'
#' The support is \eqn{x \ge 0}.
#'
#' The Expected value exists if \eqn{\xi > 1} and is equal to
#'
#' \deqn{E(X) = \theta / (\xi - 1)}
#'
#' k-th moments exist in general for \eqn{k < \xi}.
#'
#' @examples
#'
#' x <- rpareto(1000, shape = 10, scale = 5)
#' xx <- seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dpareto(xx, shape = 10, scale = 5))
#'
#' plot(xx, dpareto(xx, shape = 10, scale = 5), type = "l")
#' lines(xx, dpareto(xx, shape = 3, scale = 5), col = "red", lwd = 2)
#'
#' plot(xx, dpareto(xx, shape = 10, scale = 10), type = "l")
#' lines(xx, dpareto(xx, shape = 10, scale = 5), col = "blue", lwd = 2)
#' lines(xx, dpareto(xx, shape = 10, scale = 20), col = "red", lwd = 2)
#'
#' @references
#' [https://en.wikipedia.org/wiki/Pareto_distribution]() - named Lomax therein.
NULL

#' @rdname Pareto
#' @export
#' @return
#' `rpareto` generates random deviates.
rpareto <- function(n = 1L, shape = 0.0, scale = 1.0) {
  if (n == 0L) return(numeric())
  check_lengths(shape, scale, .len = n, .msg = "n")
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)

  scale * (runif(n)^(-1.0 / shape) - 1.0)
}

#' @rdname Pareto
#' @export
#' @return
#' `dpareto` gives the density.
dpareto <- function(x, shape = 1.0, scale = 1.0, log = FALSE) {
  if (!length(x)) return(numeric())
  n <- check_lengths(x, shape, scale)
  x <- rep_len(x, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)

  ok <- x >= 0.0
  logdens <- numeric(n)
  logdens[!ok] <- -Inf
  logdens[ok] <- log(shape[ok]) - log(scale[ok]) - (shape[ok] + 1.0) * log1p(x[ok] / scale[ok])

  if (log) logdens else exp(logdens)
}

#' @rdname Pareto
#' @export
#' @return
#' `ppareto` gives the distribution function.
ppareto <- function(q, shape = 1.0, scale = 1.0, lower.tail = TRUE, log.p = FALSE) {
  if (!length(q)) return(numeric())
  n <- check_lengths(q, shape, scale)
  q <- rep_len(q, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)

  logfbar <- -shape * log1p(pmax(q, 0.0) / scale)

  if (lower.tail) {
    f <- -expm1(logfbar)
    if (log.p) log(f) else f
  } else {
    if (log.p) logfbar else exp(logfbar)
  }
}

#' @rdname Pareto
#' @export
#' @return
#' `qpareto` gives the quantile function.
qpareto <- function(p, shape = 1.0, scale = 1.0, lower.tail = TRUE, log.p = FALSE) {
  if (!length(p)) return(numeric())
  n <- check_lengths(p, shape, scale)
  p <- rep_len(p, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)

  if (lower.tail) {
    p <- if (log.p) exp(p) else p
    lgpbar <- log1p(-p)
  } else {
    lgpbar <- if (log.p) p else log(p)
  }

  scale * expm1(-lgpbar / shape)
}

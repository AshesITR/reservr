#' @rdname GenPareto
#' @name GenPareto
#' @title The Generalized Pareto Distribution (GPD)
#'
#' @description
#' These functions provide information about the generalized Pareto distribution
#' with threshold `u`. `dgpd` gives the density, `pgpd` gives the distribution
#' function, `qgpd` gives the quantile function and `rgpd` generates random
#' deviates.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n integer number of observations.
#' @param u threshold parameter (minimum value).
#' @param sigmau scale parameter (must be positive).
#' @param xi shape parameter
#' @param log,log.p logical; if `TRUE`, probabilities/densities
#' `p` are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#' \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#'
#' @details
#' If `u`, `sigmau` or `xi` are not specified, they assume the default values of
#' `0`, `1` and `0` respectively.
#'
#' The generalized Pareto distribution has density
#'
#' \deqn{f(x) = 1 / \sigma_u (1 + \xi z)^(- 1 / \xi - 1)}
#'
#' where \eqn{z = (x - u) / \sigma_u} and \eqn{f(x) = exp(-z)} if
#' \eqn{\xi} is 0.
#' The support is \eqn{x \ge u} for \eqn{\xi \ge 0} and
#' \eqn{u \le x \le u - \sigma_u / \xi} for \eqn{\xi < 0}.
#'
#' The Expected value exists if \eqn{\xi < 1} and is equal to
#'
#' \deqn{E(X) = u + \sigma_u / (1 - \xi)}
#'
#' k-th moments exist in general for \eqn{k\xi < 1}.
#'
#' @examples
#'
#' x <- rgpd(1000, u = 1, sigmau = 0.5, xi = 0.1)
#' xx <- seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgpd(xx, u = 1, sigmau = 0.5, xi = 0.1))
#'
#' plot(xx, dgpd(xx, u = 1, sigmau = 1, xi = 0), type = "l")
#' lines(xx, dgpd(xx, u = 0.5, sigmau = 1, xi = -0.3), col = "blue", lwd = 2)
#' lines(xx, dgpd(xx, u = 1.5, sigmau = 1, xi = 0.3), col = "red", lwd = 2)
#'
#' plot(xx, dgpd(xx, u = 1, sigmau = 1, xi = 0), type = "l")
#' lines(xx, dgpd(xx, u = 1, sigmau = 0.5, xi = 0), col = "blue", lwd = 2)
#' lines(xx, dgpd(xx, u = 1, sigmau = 2, xi = 0), col = "red", lwd = 2)
#'
#' @references
#' [https://en.wikipedia.org/wiki/Generalized_Pareto_distribution]()
NULL

#' @rdname GenPareto
#' @export
#' @return
#' `rgpd` generates random deviates.
rgpd <- function(n = 1L, u = 0.0, sigmau = 1.0, xi = 0.0) {
  if (n == 0L) return(numeric())
  check_lengths(u, sigmau, xi, .len = n, .msg = "n")
  u <- rep_len(u, n)
  sigmau <- rep_len(sigmau, n)
  xi <- rep_len(xi, n)

  xi0 <- xi == 0.0
  res <- numeric(n)
  res[xi0] <- rexp(sum(xi0))
  res[!xi0] <- (runif(sum(!xi0))^(-xi[!xi0]) - 1.0) / xi[!xi0]

  u + sigmau * res
}

#' @rdname GenPareto
#' @export
#' @return
#' `dgpd` gives the density.
dgpd <- function(x, u = 0.0, sigmau = 1.0, xi = 0.0, log = FALSE) {
  if (!length(x)) return(numeric())
  n <- check_lengths(x, u, sigmau, xi)
  x <- rep_len(x, n)
  u <- rep_len(u, n)
  sigmau <- rep_len(sigmau, n)
  xi <- rep_len(xi, n)

  xi0 <- xi == 0.0

  z <- (x - u) / sigmau
  zxip1 <- 1.0 + xi * z
  support <- (z >= 0.0) & (xi0 | zxip1 >= 0.0)
  support[is.na(z)] <- FALSE

  logdens <- numeric(length(x))
  logdens[!support] <- -Inf
  logdens[xi0 & support] <- -z[xi0 & support]
  logdens[!xi0 & support] <- -(1.0 / xi[!xi0 & support] + 1.0) * # nolint (lintr bug)
    log(zxip1[!xi0 & support])
  logdens[support] <- logdens[support] - log(sigmau[support])
  logdens[is.na(z)] <- z[is.na(z)]

  if (log) logdens else exp(logdens)
}

#' @rdname GenPareto
#' @export
#' @return
#' `pgpd` gives the distribution function.
pgpd <- function(q, u = 0.0, sigmau = 1.0, xi = 0.0,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!length(q)) return(numeric())
  n <- check_lengths(q, u, sigmau, xi)
  q <- rep_len(q, n)
  u <- rep_len(u, n)
  sigmau <- rep_len(sigmau, n)
  xi <- rep_len(xi, n)

  z <- pmax(0.0, (q - u) / sigmau)
  zxip1 <- pmax(0.0, 1.0 + xi * z)
  xi0 <- xi == 0.0

  res <- numeric(length(q))
  res[xi0] <- pexp(z[xi0], lower.tail = lower.tail, log.p = log.p)
  if (lower.tail) {
    f <- 1 - zxip1[!xi0]^(-1.0 / xi[!xi0])
    res[!xi0] <- if (log.p) log(f) else f
  } else {
    logfbar <- -1.0 / xi[!xi0] * log(zxip1[!xi0])
    res[!xi0] <- if (log.p) logfbar else exp(logfbar)
  }

  res
}

#' @rdname GenPareto
#' @export
#' @return
#' `qgpd` gives the quantile function.
qgpd <- function(p, u = 0.0, sigmau = 1.0, xi = 0.0,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!length(p)) return(numeric())
  n <- check_lengths(p, u, sigmau, xi)
  p <- rep_len(p, n)
  u <- rep_len(u, n)
  sigmau <- rep_len(sigmau, n)
  xi <- rep_len(xi, n)

  xi0 <- xi == 0.0

  res <- numeric(length(p))
  res[xi0] <- qexp(p[xi0], lower.tail = lower.tail, log.p = log.p)
  if (log.p) p[!xi0] <- exp(p[!xi0])
  if (lower.tail) p[!xi0] <- 1.0 - p[!xi0]
  res[!xi0] <- (p[!xi0]^(-xi[!xi0]) - 1.0) / xi[!xi0]

  u + sigmau * res
}

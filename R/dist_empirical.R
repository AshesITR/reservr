#' Empirical distribution
#'
#' Creates an empirical distribution object from a sample.
#' Assumes iid. samples. `with_params` should **not** be used with this
#' distribution because estimation of the relevant indicators happens during
#' construction.
#'
#'  * `sample()` samples iid. from `sample`. This approach is similar to
#'    bootstrapping.
#'  * `density()` evaluates a kernel density estimate, approximating with zero
#'    outside of the known support. This estimate is either obtained using
#'    [stats::density] or [logKDE::logdensity_fft], depending on `positive`.
#'  * `probability()` evaluates the empirical cumulative density function
#'    obtained by [stats::ecdf].
#'  * `quantile()` evaluates the empirical quantiles using [stats::quantile]
#'  * `hazard()` estimates the hazard rate using the density estimate and the
#'    empirical cumulative density function: `h(t) = df(t) / (1 - cdf(t))`.
#'
#' @param sample Sample to build the empirical distribution from
#' @param positive Is the underlying distribution known to be positive?
#' This will effect the density estimation procedure.
#' `positive = FALSE` uses a kernel density estimate produced by [density()],
#' `positive = TRUE` uses a log-kernel density estimate produced by
#' [logKDE::logdensity_fft()]. The latter can improve density estimation near
#' zero.
#' @param bw Bandwidth parameter for density estimation. Passed to the density
#' estimation function selected by `positive`.
#'
#' @return An `EmpiricalDistribution` object.
#'
#' @export
#'
#' @examples
#' x <- rexp(20, rate = 1)
#' dx <- dist_empirical(sample = x, positive = TRUE)
#'
#' y <- rnorm(20)
#' dy <- dist_empirical(sample = y)
#'
#' plot_distributions(
#'   exponential = dx,
#'   normal = dy,
#'   .x = seq(-3, 3, length.out = 100)
#' )
#'
#' @family Distributions
dist_empirical <- function(sample, positive = FALSE, bw = "nrd0") {
  if (positive) {
    check_installed("logKDE",
                    "for dist_empirical(..., positive = TRUE).")
    densest <- muffle_warning(
      logKDE::logdensity_fft(sample, bw = bw),
      "Auto-range choice cut-off at 0."
    )
  } else {
    densest <- stats::density(sample, bw = bw)
  }
  densfun <- approxfun(densest, yleft = 0, yright = 0)
  EmpiricalDistribution$new(
    x = as.list(sample),
    df = densfun,
    cdf = ecdf(sample)
  )
}

EmpiricalDistribution <- distribution_class(
  name = "empirical",
  type = "discrete",
  params = list(
    x = list(),
    df = list(),
    cdf = list()
  ),
  sample = function(n, params) {
    base::sample(unlist(params$x), n, replace = TRUE)
  },
  density = function(x, log = FALSE, params) {
    res <- params$df(x)
    if (log) log(res) else res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    res <- params$cdf(q)
    if (!lower.tail) res <- 1.0 - res
    if (log.p) log(res) else res
  },
  quantile = function(p, lower.tail = TRUE, log.p = FALSE, params) {
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1.0 - p
    quantile(unlist(params$x), p) # FIXME this only works if params$x is a list of scalars
  },
  is_in_support = function(x, params) {
    x %in% unlist(params$x) # FIXME this only works if params$x is a list of scalars
  }
)

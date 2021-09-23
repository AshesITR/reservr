#' Compute weighted moments
#'
#' @param x Observations
#' @param w Case weights (optional)
#' @param n Number of moments to calculate
#' @param center Calculate centralized moments (default) or noncentralized
#' moments, i.e. E((X - E(X))^k) or E(X^k).
#'
#' @return A vector of length `n` where the `k`th entry is the `k`th weighted
#' moment of `x` with weights `w`. If `center` is `TRUE` the moments are
#' centralized, i.e. E((X - E(X))^k). The first moment is never centralized.
#' The moments are scaled with 1 / sum(w), so they are not de-biased.
#'
#' e.g. the second central weighted moment
#' `weighted_moment(x, w)[2L]`
#' is equal to
#' `var(rep(x, w)) * (sum(w) - 1) / sum(w)`
#' for integer `w`
#' @export
#' @family weighted statistics
#'
#' @examples
#' weighted_moments(rexp(100))
#' weighted_moments(c(1, 2, 3), c(1, 2, 3))
#' c(mean(rep(1:3, 1:3)), var(rep(1:3, 1:3)) * 5 / 6)
weighted_moments <- function(x, w, n = 2L, center = TRUE) {
  if (missing(w)) {
    w <- rep_len(1.0, length(x))
  }
  mu <- weighted.mean(x, w)
  vapply(seq_len(n), function(k) {
    if (k == 1L) {
      mu
    } else if (center) {
      weighted.mean((x - mu)^k, w)
    } else {
      weighted.mean(x^k, w)
    }
  }, numeric(1L))
}

#' Compute weighted quantiles
#'
#' @param x Observations
#' @param w Case weights (optional)
#' @param probs Quantiles to calculate
#'
#' @return A vector the same length as `probs` with the corresponding weighted
#' quantiles of `x` with weight `w`. For integer weights, this is equivalent to
#' `quantile(rep(x, w), probs)`
#' @export
#' @family weighted statistics
#'
#' @examples
#' weighted_median(1:6)
#' weighted_median(1:3, c(1, 4, 9))
#' weighted_median(1:3, c(9, 4, 1))
#'
#' weighted_quantile(1:3, c(1, 4, 9), seq(0.0, 1.0, by = 0.25))
#' quantile(rep(1:3, c(1, 4, 9)), seq(0.0, 1.0, by = 0.25))
weighted_quantile <- function(x, w, probs) {
  if (missing(w)) return(unname(quantile(as.numeric(x), probs)))
  o <- order(x)
  wc <- cumsum(w[o])

  n <- length(x)

  targets <- 1.0 + (wc[n] - 1.0) * probs
  i <- findInterval(targets, wc)

  w_low <- wc[i]
  w_low[i == 0L] <- wc[1L]
  w_hi <- wc[i + 1L]
  w_hi[i == n] <- wc[n]

  h <- (w_hi - targets) / (w_hi - w_low)
  h[i %in% c(0L, n)] <- 0.0
  h[targets - w_low > 1.0] <- 0.0

  h * x[o[pmax(1L, i)]] + (1 - h) * x[o[pmax(1L, pmin(i + 1L, n))]]
}


#' @return The weighted median of `x` with weights `w`.
#' For integer weights, this is equivalent to `median(rep(x, w))`
#'
#' @export
#' @rdname weighted_quantile
weighted_median <- function(x, w) {
  weighted_quantile(x = x, w = w, probs = 0.5)
}

#' Compute weighted tabulations
#'
#' Computes the sum of `w` grouped by `bin`.
#' If `w` is missing the result is equivalent to
#' [`tabulate(bin, nbins)`][base::tabulate()]
#'
#' @param bin An integer vector with values from `1L` to `nbins`
#' @param w Weights per entry in `bin`.
#' @param nbins Number of bins
#'
#' @return A vector with length `nbins` where the `i`th result is equal to
#' `sum(w[bin == i])` or `sum(bin == i)` if `w` is missing.
#' For integer weights, this is equivalent to `tabulate(rep(bin, w), nbins)`.
#' @export
#' @family weighted statistics
#'
#' @examples
#' weighted_tabulate(c(1, 1, 2))
#' weighted_tabulate(c(1, 1, 2), nbins = 3L)
#' weighted_tabulate(c(1, 1, 2), w = c(0.5, 0.5, 1), nbins = 3L)
weighted_tabulate <- function(bin, w, nbins = max(1L, bin, na.rm = TRUE)) {
  if (missing(w)) return(tabulate(bin, nbins))
  res <- numeric(nbins)
  for (b in seq_len(nbins)) {
    res[b] <- sum(w[bin == b])
  }
  res
}

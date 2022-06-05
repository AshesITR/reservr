#' Transition functions for blended distributions
#'
#' @param x Points to evaluate at
#' @param u Sorted vector of blending thresholds, or rowwise sorted matrix of blending thresholds
#' @param eps Corresponding vector or matrix of blending bandwidths.
#' Must be positive and the same dimensions as `u`, or scalar.
#' No rowwise blending regions (u - eps, u + eps) may overlap.
#' @param .gradient Also evaluate the gradient with respect to `x`?
#' @param .extend_na Extend out-of range transitions by the last in-range value (i.e. the corresponding u) or by NA?
#'
#' @return `blended_transition` returns a matrix with `length(x)` rows and `length(u) + 1` columns containing the
#' transformed values for each of the blending components.
#' If `.gradient` is TRUE, an attribute `"gradient"` is attached with the same dimensions, containing the derivative
#' of the respective transition component with respect to `x`.
#'
#' @examples
#' library(ggplot2)
#' xx <- seq(from = 0, to = 20, length.out = 101)
#' blend_mat <- blended_transition(xx, u = 10, eps = 3, .gradient = TRUE)
#' ggplot(
#'   data.frame(
#'     x = rep(xx, 2L),
#'     fun = rep(c("p", "q"), each = length(xx)),
#'     y = as.numeric(blend_mat),
#'     relevant = c(xx <= 13, xx >= 7)
#'   ),
#'   aes(x = x, y = y, color = fun, linetype = relevant)
#' ) %+%
#'   geom_line() %+%
#'   theme_bw() %+%
#'   theme(
#'     legend.position = "bottom", legend.box = "horizontal"
#'   ) %+%
#'   guides(color = guide_legend(direction = "horizontal", title = ""), linetype = guide_none()) %+%
#'   scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 3))
#'
#' ggplot(
#'   data.frame(
#'     x = rep(xx, 2L),
#'     fun = rep(c("p'", "q'"), each = length(xx)),
#'     y = as.numeric(attr(blend_mat, "gradient")),
#'     relevant = c(xx <= 13, xx >= 7)
#'   ),
#'   aes(x = x, y = y, color = fun, linetype = relevant)
#' ) %+%
#'   geom_line() %+%
#'   theme_bw() %+%
#'   theme(
#'     legend.position = "bottom", legend.box = "horizontal"
#'   ) %+%
#'   guides(color = guide_legend(direction = "horizontal", title = ""), linetype = guide_none()) %+%
#'   scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 3))
#'
#' @export
blended_transition <- function(x, u, eps, .gradient = FALSE, .extend_na = FALSE) {
  n <- length(x)

  if (is.matrix(u)) {
    k <- ncol(u)
  } else if (is.vector(u)) {
    k <- length(u)
    u <- matrix(data = u, nrow = n, ncol = k, byrow = TRUE)
  }
  if (is.vector(eps)) {
    eps <- matrix(data = eps, nrow = n, ncol = k, byrow = TRUE)
  }

  assert_that(
    is_bool(.gradient),
    msg = "`.gradient` must be a bool."
  )
  assert_that(
    is_bool(.extend_na),
    msg = "`.extend_na` must be a bool."
  )
  u_diffs <- matrixStats::rowDiffs(u)
  assert_that(
    is.numeric(u),
    all(u_diffs >= 0),
    msg = "`u` must be a (rowwise) non-decreasing real vector or matrix."
  )
  assert_that(
    is.numeric(eps),
    all(eps >= 0),
    msg = "`eps` must be a non-negative real vector or matrix."
  )
  if (!all(u_diffs >= eps[, -k] + eps[, -1L])) {
    bad_i <- which(u_diffs >= eps[-k] + eps[-1L])[1L]
    bad_u1 <- u[bad_i]
    bad_u2 <- u[bad_i + n]
    bad_eps1 <- eps[bad_i]
    bad_eps2 <- eps[bad_i + n]
    stop(sprintf(
      "Blending regions must not overlap. (%g \u00b1 %g) overlaps (%g \u00b1 %g)",
      bad_u1, bad_eps1, bad_u2, bad_eps2
    ))
  }

  res <- matrix(data = x, nrow = n, ncol = k + 1L)
  if (.gradient) {
    dres <- matrix(data = 1.0, nrow = n, ncol = k + 1L)
  }

  for (i in seq_len(k)) {
    u_curr <- u[, i]
    e_curr <- eps[, i]
    b_curr <- u_curr - e_curr < x & x < u_curr + e_curr
    lo_curr <- x <= u_curr - e_curr
    hi_curr <- x >= u_curr + e_curr

    res[hi_curr, i] <- u_curr[hi_curr]
    res[lo_curr, i + 1L] <- u_curr[lo_curr]

    blend_curr <- e_curr[b_curr] / pi * cos(0.5 * pi * (x[b_curr] - u_curr[b_curr]) / e_curr[b_curr])
    res[b_curr, i] <- 0.5 * (x[b_curr] + u_curr[b_curr] - e_curr[b_curr]) + blend_curr
    res[b_curr, i + 1L] <- 0.5 * (x[b_curr] + u_curr[b_curr] + e_curr[b_curr]) - blend_curr

    if (.gradient) {
      dres[hi_curr, i] <- 0.0
      dres[lo_curr, i + 1L] <- 0.0

      dblend_curr <- 0.5 * sin(0.5 * pi * (x[b_curr] - u_curr[b_curr]) / e_curr[b_curr])
      dres[b_curr, i] <- 0.5 - dblend_curr
      dres[b_curr, i + 1L] <- 0.5 + dblend_curr
    }

    if (.extend_na) {
      na_lo <- x < u_curr - e_curr
      na_hi <- x > u_curr + e_curr

      res[na_hi, i] <- NA_real_
      res[na_lo, i + 1L] <- NA_real_

      if (.gradient) {
        dres[na_hi, i] <- NA_real_
        dres[na_lo, i + 1L] <- NA_real_
      }
    }
  }

  if (.gradient) {
    attr(res, "gradient") <- dres
  }

  res
}

#' @rdname blended_transition
#'
#' @param .component Component index (up to `length(u) + 1`) to invert.
#'
#' @return `blended_transition_inv` returns a vector with `length(x)` values containing the inverse of the transformed
#' values for the `.component`th blending component.
#' @export
blended_transition_inv <- function(x, u, eps, .component) {
  n <- length(x)

  if (is.matrix(u)) {
    k <- ncol(u)
  } else if (is.vector(u)) {
    k <- length(u)
    u <- matrix(data = u, nrow = n, ncol = k, byrow = TRUE)
  }
  if (is.vector(eps)) {
    eps <- matrix(data = eps, nrow = n, ncol = k, byrow = TRUE)
  }

  assert_that(
    is_scalar_integerish(.component),
    .component >= 1L,
    .component <= k + 1L,
    msg = sprintf("`.component` must be an index from 1 to %d.", k + 1L)
  )
  u_diffs <- matrixStats::rowDiffs(u)
  assert_that(
    is.numeric(u),
    all(u_diffs >= 0),
    msg = "`u` must be a (rowwise) non-decreasing real vector or matrix."
  )
  assert_that(
    is.numeric(eps),
    all(eps >= 0),
    msg = "`eps` must be a non-negative real vector or matrix."
  )
  if (!all(u_diffs >= eps[, -k] + eps[, -1L])) {
    bad_i <- which(u_diffs >= eps[-k] + eps[-1L])[1L]
    bad_u1 <- u[bad_i]
    bad_u2 <- u[bad_i + n]
    bad_eps1 <- eps[bad_i]
    bad_eps2 <- eps[bad_i + n]
    stop(sprintf(
      "Blending regions must not overlap. (%g \u00b1 %g) overlaps (%g \u00b1 %g)",
      bad_u1, bad_eps1, bad_u2, bad_eps2
    ))
  }

  res <- x

  if (.component > 1L) {
    u_lo <- u[, .component - 1L]
    e_lo <- eps[, .component - 1L]

    res[x < u_lo] <- NA_real_

    t_lo <- u_lo <= x & x < u_lo + e_lo
    z <- (x[t_lo] - u_lo[t_lo]) / e_lo[t_lo] - 0.5

    inv_tr <- xcx_inv(pi * z) * 2.0 / pi
    res[t_lo] <- u_lo[t_lo] + e_lo[t_lo] * inv_tr
  }

  if (.component <= k) {
    u_hi <- u[, .component]
    e_hi <- eps[, .component]

    res[x > u_hi] <- NA_real_

    t_hi <- u_hi - e_hi < x & x <= u_hi
    z <- (x[t_hi] - u_hi[t_hi]) / e_hi[t_hi] + 0.5

    inv_tr <- -xcx_inv(-pi * z) * 2.0 / pi
    res[t_hi] <- u_hi[t_hi] + e_hi[t_hi] * inv_tr
  }

  res
}

blended_transition_fst <- function(x, u_lo, u_hi, e_lo, e_hi, blend_left, blend_right) {
  xout <- x
  dout <- rep_len(1.0, length(x))

  if (blend_left) {
    i_low <- x < u_lo + e_lo
    if (length(u_lo) > 1L) u_lo <- u_lo[i_low]
    if (length(e_lo) > 1L) e_lo <- e_lo[i_low]
    blend_curr <- e_lo / pi * cospi(0.5 * (x[i_low] - u_lo) / e_lo)
    dblend_curr <- -0.5 * sinpi(0.5 * (x[i_low] - u_lo) / e_lo)
    xout[i_low] <- 0.5 * (x[i_low] + u_lo + e_lo) - blend_curr
    dout[i_low] <- 0.5 - dblend_curr
  }

  if (blend_right) {
    i_high <- x > u_hi - e_hi
    if (length(u_hi) > 1L) u_hi <- u_hi[i_high]
    if (length(e_hi) > 1L) e_hi <- e_hi[i_high]
    blend_curr <- e_hi / pi * cospi(0.5 * (x[i_high] - u_hi) / e_hi)
    dblend_curr <- -0.5 * sinpi(0.5 * (x[i_high] - u_hi) / e_hi)
    xout[i_high] <- 0.5 * (x[i_high] + u_hi - e_hi) + blend_curr
    dout[i_high] <- 0.5 + dblend_curr
  }

  list(xout, dout)
}

blended_transition_finv <- function(x, u_lo, u_hi, e_lo, e_hi, blend_left, blend_right) {
  res <- x

  if (blend_left) {
    res[x < u_lo] <- NA_real_

    t_lo <- u_lo <= x & x < u_lo + e_lo
    z <- (x[t_lo] - u_lo[t_lo]) / e_lo[t_lo] - 0.5

    inv_tr <- xcx_inv(pi * z) * 2.0 / pi
    res[t_lo] <- u_lo[t_lo] + e_lo[t_lo] * inv_tr
  }

  if (blend_right) {
    res[x > u_hi] <- NA_real_

    t_hi <- u_hi - e_hi < x & x <= u_hi
    z <- (x[t_hi] - u_hi[t_hi]) / e_hi[t_hi] + 0.5

    inv_tr <- -xcx_inv(-pi * z) * 2.0 / pi
    res[t_hi] <- u_hi[t_hi] + e_hi[t_hi] * inv_tr
  }

  res
}

# nodes and weights of the gauss-kronrod quadrature
GK_QUAD <- list(
  nodes = c(
    0.991455371120813, -0.991455371120813, 0.949107912342759,
    -0.949107912342759, 0.864864423359769, -0.864864423359769,
    0.741531185599394, -0.741531185599394, 0.586087235467691,
    -0.586087235467691, 0.405845151377397, -0.405845151377397,
    0.207784955007898, -0.207784955007898, 0.000000000000000
  ),
  weights = c(
    0.022935322010529, 0.022935322010529, 0.063092092629979,
    0.063092092629979, 0.104790010322250, 0.104790010322250,
    0.140653259715525, 0.140653259715525, 0.169004726639267,
    0.169004726639267, 0.190350578064785, 0.190350578064785,
    0.204432940075298, 0.204432940075298, 0.209482141084728
  ),
  weights7 = c(
    0.129484966168870, 0.129484966168870, 0.279705391489277,
    0.279705391489277, 0.381830050505119, 0.381830050505119,
    0.417959183673469
  ),
  idx7 = c(3, 4, 7, 8, 11, 12, 15)
)

#' Adaptive Gauss-Kronrod Quadrature for multiple limits
#'
#' Integrates fun over the bounds \[ lower, upper \] vectorized over `lower` and
#' `upper`. Vectorized list structures of parameters can also be passed.
#'
#' The integration error is estimated by the Gauss-Kronrod quadrature as the
#' absolute difference between the 7-point quadrature and the 15-point
#' quadrature. Integrals that did not converge will be bisected at the midpoint.
#' The `params` object will be recursively subsetted on all numeric vectors with
#' the same length as the number of observations.
#'
#' @param fun A function to integrate.
#' Must be vectorized and take one or two arguments, the first being points to
#' evaluate at and the second (optionally) being parameters to apply.
#' It must return a numeric vector the same length as its first input.
#'
#' Currently, infinite bounds are not supported.
#'
#' @param lower,upper Integration bounds. Must have the same length.
#' @param params Parameters to pass as a second argument to `fun`.
#' The actual parameters must have the same length as the number of integrals to
#' compute. Can be a possibly nested list structures containing numeric vectors.
#' Alternatively, can be a matrix with the same number of rows as the number of integrals to compute.
#'
#' @param .tolerance Absolute element-wise tolerance.
#' @param .max_iter Maximum number of iterations. The number of
#' integration intervals will be at most `length(lower) * .max_iter`. Therefor the maximum
#' number of function evaluations per integration interval will be
#' `15 * .max_iter`.
#'
#' @return A vector of integrals with the i-th entry containing an approximation
#' of
#'
#' int_{lower\[i\]}^{upper\[i\]} fun(t, pick_params_at(params, i)) dt
#'
#' @examples
#' # Argument recycling and parallel integration of two intervals
#' integrate_gk(sin, 0, c(pi, 2 * pi))
#'
#' dist <- dist_exponential()
#' integrate_gk(
#'   function(x, p) dist$density(x, with_params = p),
#'   lower = 0, upper = 1:10,
#'   params = list(rate = 1 / 1:10)
#' )
#' dist$probability(1:10, with_params = list(rate = 1 / 1:10))
#'
#' @export
integrate_gk <- function(fun, lower, upper, params = list(),
                         .tolerance = .Machine$double.eps^0.25,
                         .max_iter = 100L) {
  fun <- match.fun(fun)
  assert_that(is.function(fun), msg = "`fun` must be a function.")
  assert_that(is.numeric(lower), msg = "`lower` must be numeric.")
  assert_that(is.numeric(upper), msg = "`upper` must be numeric.")
  assert_that(!anyNA(lower), !anyNA(upper),
              msg = "Bounds may not contain NA / NaN.")
  assert_that(all(is.finite(lower)), all(is.finite(upper)),
              msg = "Infinite bounds are not supported yet.")
  n <- check_lengths(lower, upper)
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)
  n_args <- length(formals(args(fun)))
  assert_that(
    n_args %in% c(1L, 2L),
    msg = "`fun` must be a function of one or two arguments."
  )

  if (n_args == 1L && !missing(params) && !is.null(params)) {
    warning("`params` is specified but `fun` takes only one argument.",
            " Ignoring `params`.")
    params <- list()
  }

  if (is.null(params)) {
    params <- list()
  }

  if (n_args == 1L) {
    fun_noargs <- fun
    fun <- function(x, params) fun_noargs(x)
  }

  if (is.list(params)) {
    do_integrate_gk_lst(fun, lower, upper, params, .tolerance, .max_iter, FALSE)$value
  } else {
    do_integrate_gk_mat(fun, lower, upper, params, .tolerance, .max_iter, FALSE)$value
  }
}

#' Determine probability of reporting under a Poisson arrival Process
#'
#' Determines the probability that claims occuring under a Poisson process with
#' arrival intensity `expo` and reporting delay distribution `dist` during the
#' time between `t_min` and `t_max` are reported between `tau_min` and
#' `tau_max`.
#'
#' @param dist A reporting delay Distribution.
#' @param intervals A data frame with columns `xmin`, `xmax`, `tmin`, `tmax`.
#' Claims occur within `[xmin, xmax]` and be reported within `[tmin, tmax]`.
#' @param expo Poisson intensity. If given, must be a vectorised function that
#' yields the intensity of the claim arrival process at a specified time.
#' `expo = NULL` is equivalent to a constant intensity function. `expo` is only
#' relevant up to a multiplicative constant.
#' @param with_params Parameters of `dist` to use. Can be a parameter set with
#' different values for each interval.
#' @param .try_compile Try compiling the distributions probability function to speed up integration?
#' @inheritParams integrate_gk
#'
#' @details The reporting probability is given by
#'
#' P(x + d in \[tmin, tmax\] | x in \[xmin, xmax\])
#'    = E(P(x + d in \[tmin, tmax\] | x) | x in \[xmin, xmax\]) /
#'          P(x in \[xmin, xmax\])
#'    = int_\[xmin, xmax\] expo(x) P(x + d in \[tmin, tmax\]) dx
#'    = int_\[xmin, xmax\] expo(x) P(d in \[tmin - x, tmax - x\]) dx /
#'          int_\[xmin, xmax\] expo(x) dx
#'
#' `prob_report` uses [integrate_gk()] to compute the two integrals.
#'
#' @return A vector of reporting probabilities, with one entry per row of `intervals`.
#'
#' @examples
#' dist <- dist_exponential()
#' ints <- data.frame(
#'   xmin = 0,
#'   xmax = 1,
#'   tmin = seq_len(10) - 1.0,
#'   tmax = seq_len(10)
#' )
#' params <- list(rate = rep(c(1, 0.5), each = 5))
#'
#' prob_report(dist, ints, with_params = params)
#'
#' @export
prob_report <- function(dist, intervals, expo = NULL, with_params = list(),
                        .tolerance = .Machine$double.eps^0.5, .max_iter = 100L,
                        .try_compile = TRUE) {
  if (!is.null(expo)) {
    total_expo <- integrate_gk(
      function(x, params) {
        expo(x)
      },
      lower = intervals$xmin,
      upper = intervals$xmax,
      params = list(),
      .tolerance = .tolerance,
      .max_iter = .max_iter
    )
  } else {
    total_expo <- intervals$xmax - intervals$xmin
  }

  prob_int <- if (.try_compile) tryCatch(dist$compile_probability_interval(), error = function(e) NULL)
  if (!is.null(prob_int)) {
    wp_matrix <- flatten_params_matrix(with_params)
    if (nrow(wp_matrix) == 0L) wp_matrix <- matrix(nrow = nrow(intervals), ncol = 0L)
    exp_nclaims_reported <- integrate_gk(
      if (!is.null(expo)) {
        function(x, params) {
          p_report <- prob_int(
            qmin = params[, 1L] - x,
            qmax = params[, 2L] - x,
            param_matrix = params[, -(1L:2L)]
          )
          expo(x) * p_report
        }
      } else {
        function(x, params) {
          prob_int(
            qmin = params[, 1L] - x,
            qmax = params[, 2L] - x,
            param_matrix = params[, -(1L:2L)]
          )
        }
      },
      lower = intervals$xmin,
      upper = intervals$xmax,
      params = cbind(intervals$tmin, intervals$tmax, wp_matrix),
      .tolerance = .tolerance,
      .max_iter = .max_iter
    )
  } else {
    exp_nclaims_reported <- integrate_gk(
      function(x, params) {
        repdel_max <- params$tmax - x
        repdel_min <- params$tmin - x
        pars <- params$dist_params
        p_report <- dist$probability(repdel_max, with_params = pars) -
          dist$probability(repdel_min, with_params = pars)

        if (!dist$is_continuous()) {
          disc <- dist$is_discrete_at(repdel_min, with_params = pars)
          p_report[disc] <- p_report[disc] + dist$density(
            repdel_min[disc], with_params = pick_params_at(pars, disc)
          )
        }

        if (!is.null(expo)) expo(x) * p_report else p_report
      },
      lower = intervals$xmin,
      upper = intervals$xmax,
      params = list(
        tmax = intervals$tmax,
        tmin = intervals$tmin,
        dist_params = with_params
      ),
      .tolerance = .tolerance,
      .max_iter = .max_iter
    )
  }


  exp_nclaims_reported / total_expo
}

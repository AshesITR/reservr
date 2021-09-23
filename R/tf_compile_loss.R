#' Compile log-likelihood functions for use with TensorFlow
#'
#' The TensorFlow functions take two inputs
#'  * `obs`, a Tensor consisting of six columns in order
#'    x, xmin, xmax, tmin, tmax, w. (see [trunc_obs()])
#'  * `inputs`, a list of parameter input tensors.
#'
#' Truncation support uses the columns `tmin` and `tmax` for the truncation
#' intervals.
#' Censoring support uses the columns `xmin` and `xmax` for the censoring
#' intervals.
#'
#' @param dist A distribution object with placeholders.
#' Non-placeholder parameters are assumed to be constant for the purpose of
#' training.
#'
#' @return `tf_compile_loss` returns a Tensorflow function computing
#' the log-Likelihood of an uncensored sample without truncation.
#'
#' @noRd
tf_compile_loss <- function(dist) {
  constants <- dist$tf_make_constants()
  logdens <- dist$tf_logdensity()

  function(obs, inputs) {
    args <- tf_merge_constants(inputs, constants)
    x <- obs[, 1L]
    log_dens <- logdens(x, args)
    -log_dens * obs[, 6L] # weight
  }
}

#' @noRd
#'
#' @return `tf_compile_loss_trunc_cens` returns a Tensorflow function computing
#' the log-Likelihood of a censored sample with truncation.
tf_compile_loss_trunc_cens <- function(dist) {
  tf <- tensorflow::tf

  constants <- dist$tf_make_constants()
  logdens <- dist$tf_logdensity()
  logprob <- dist$tf_logprobability()

  function(obs, inputs) {
    args <- tf_merge_constants(inputs, constants)

    x <- obs[, 1L]
    xmin <- obs[, 2L]
    xmax <- obs[, 3L]
    tmin <- obs[, 4L]
    tmax <- obs[, 5L]

    x_safe <- tf$where(tf$math$is_nan(x), K$one_half * xmin + K$one_half * xmax, x)
    xmin_safe <- tf$where(tf$math$is_nan(x), xmin, K$neg_inf)
    xmax_safe <- tf$where(tf$math$is_nan(x), xmax, K$inf)

    log_dens <- tf$where(
      tf$math$is_nan(x),
      logprob(xmin_safe, xmax_safe, args),
      logdens(x_safe, args)
    )

    trunc <- logprob(tmin, tmax, args)

    (trunc - log_dens) * obs[, 6L] # weight
  }
}

#' @noRd
#'
#' @return `tf_compile_loss_trunc` returns a Tensorflow function computing
#' the log-Likelihood of an uncensored sample with truncation.
tf_compile_loss_trunc <- function(dist) {
  constants <- dist$tf_make_constants()
  logdens <- dist$tf_logdensity()
  logprob <- dist$tf_logprobability()

  function(obs, inputs) {
    args <- tf_merge_constants(inputs, constants)

    x <- obs[, 1L]
    tmin <- obs[, 4L]
    tmax <- obs[, 5L]

    log_dens <- logdens(x, args)
    trunc <- logprob(tmin, tmax, args)

    (trunc - log_dens) * obs[, 6L] # weight
  }
}

#' @noRd
#'
#' @return `tf_compile_loss_cens` returns a Tensorflow function computing
#' the log-Likelihood of a censored sample without truncation.
tf_compile_loss_cens <- function(dist) {
  tf <- tensorflow::tf

  constants <- dist$tf_make_constants()
  logdens <- dist$tf_logdensity()
  logprob <- dist$tf_logprobability()

  function(obs, inputs) {
    args <- tf_merge_constants(inputs, constants)

    x <- obs[, 1L]
    xmin <- obs[, 2L]
    xmax <- obs[, 3L]

    x_safe <- tf$where(tf$math$is_nan(x), K$one_half * xmin + K$one_half * xmax, x)
    xmin_safe <- tf$where(tf$math$is_nan(x), xmin, K$neg_inf)
    xmax_safe <- tf$where(tf$math$is_nan(x), xmax, K$inf)

    log_dens <- tf$where(
      tf$math$is_nan(x),
      logprob(xmin_safe, xmax_safe, args),
      logdens(x_safe, args)
    )

    -log_dens * obs[, 6L] # weight
  }
}

#' Predict individual distribution parameters
#'
#' @param object A compiled and trained `reservr_keras_model`.
#' @param data Input data compatible with the model.
#' @param ... ignored
#'
#' @return A parameter list suitable for the `with_params` argument of the distribution family used for the model.
#' Contains one set of parameters per row in `data`.
#'
#' @examples
#' if (keras::is_keras_available()) {
#'   dist <- dist_exponential()
#'   params <- list(rate = 1.0)
#'   N <- 100L
#'   rand_input <- runif(N)
#'   x <- dist$sample(N, with_params = params)
#'
#'   tf_in <- keras::layer_input(1L)
#'   mod <- tf_compile_model(
#'     inputs = list(tf_in),
#'     intermediate_output = tf_in,
#'     dist = dist,
#'     optimizer = keras::optimizer_adam(),
#'     censoring = FALSE,
#'     truncation = FALSE
#'   )
#'
#'   tf_fit <- fit(
#'     object = mod,
#'     x = k_matrix(rand_input),
#'     y = x,
#'     epochs = 10L,
#'     callbacks = list(
#'       callback_debug_dist_gradients(mod, k_matrix(rand_input), x)
#'     )
#'   )
#'
#'   tf_preds <- predict(mod, data = k_matrix(rand_input))
#' }
#'
#' @export
predict.reservr_keras_model <- function(object, data, ...) {
  keras_preds <- object$model(data)
  keras_preds <- object$output_splitter(keras_preds)
  keras_preds <- object$output_inflater(keras_preds)
  as_params(keras_preds)
}

#' @importFrom generics fit
#' @export
generics::fit

#' Fit a neural network based distribution model to data
#'
#' This function delegates most work to [keras::fit.keras.engine.training.Model()] and performs additional consistency
#' checks to make sure [tf_compile_model()] was called with the appropriate options to support fitting the observations
#' `y` as well as automatically converting `y` to a n x 6 matrix needed by the compiled loss function.
#'
#' Additionally, the default `batch_size` is `min(nrow(y), 10000)` instead of keras default of `32` because the latter
#' is a very bad choice for fitting most distributions since the involved loss is much less stable than typical losses
#' used in machine learning, leading to divergence for small batch sizes.
#'
#' @param object A compiled `reservr_keras_model` as obtained by [tf_compile_model()].
#' @param x A list of input tensors (predictors)
#' @param y A `trunc_obs` tibble of observed outcomes, or something convertible via [as_trunc_obs()].
#' @inheritParams keras::fit.keras.engine.training.Model
#' @param ... Unused. If old arguments are supplied, an error message will be raised informing how to fix the issue.
#'
#' @return A `history` object that contains all information collected during training.
#' The model object will be updated in-place as a side-effect.
#'
#' @seealso predict.reservr_keras_model tf_compile_model keras::fit.keras.engine.training.Model
#'
#' @examples
#' dist <- dist_exponential()
#' params <- list(rate = 1.0)
#' N <- 100L
#' rand_input <- runif(N)
#' x <- dist$sample(N, with_params = params)
#'
#' if (keras::is_keras_available()) {
#'   tf_in <- keras::layer_input(1L)
#'   mod <- tf_compile_model(
#'     inputs = list(tf_in),
#'     intermediate_output = tf_in,
#'     dist = dist,
#'     optimizer = keras::optimizer_adam(),
#'     censoring = FALSE,
#'     truncation = FALSE
#'   )
#'
#'   tf_fit <- fit(
#'     object = mod,
#'     x = k_matrix(rand_input),
#'     y = x,
#'     epochs = 10L,
#'     callbacks = list(
#'       callback_debug_dist_gradients(mod, k_matrix(rand_input), x, keep_grads = TRUE)
#'     )
#'   )
#' }
#'
#' @export
fit.reservr_keras_model <- function(object, x, y, batch_size = NULL, epochs = 10,
                                    verbose = getOption("keras.fit_verbose", default = 1), callbacks = NULL,
                                    view_metrics = getOption("keras.view_metrics", default = "auto"),
                                    validation_split = 0, validation_data = NULL, shuffle = TRUE,
                                    class_weight = NULL, sample_weight = NULL, initial_epoch = 0,
                                    steps_per_epoch = NULL, validation_steps = NULL, ...) {
  check_installed(c("tensorflow", "keras"))
  old_args <- c(
    "data", "obs", "n_epochs", "trace", ".debug_gradients", ".lr_decay", ".lr_patience", ".lr_min", ".lr_delta_rel",
    ".lr_delta_abs"
  )
  extra_args <- list(...)
  if (length(extra_args) && any(old_args %in% names(extra_args))) {
    warns <- unique(vapply(
      intersect(old_args, names(extra_args)),
      function(old_arg) {
        switch(
          old_arg,
          data = "`data` has been renamed to `x`.",
          obs = "`obs` has been renamed to `y`.",
          n_epochs = "`n_epochs` has been renamed to `epochs`.",
          trace = "Write a custom callback if `trace` is TRUE or > 1.",
          .debug_gradients =
            "Use `callback_debug_dist_gradients(object, x, y)` in `callbacks`. WARNING: Degrades performance!",
          .lr_decay = ,
          .lr_patience = ,
          .lr_min = ,
          .lr_delta_rel = ,
          .lr_delta_abs = "Use `callback_adaptive_lr(\"loss\", ...)` in `callbacks`."
        )
      },
      character(1L)
    ))
    stop("Old function arguments detected:\n", paste0(" - ", warns, collapse = "\n"))
  } else if (length(extra_args)) {
    warning("Unused arguments: ", paste0("'", names(extra_args), "'", collapse = ", "))
  }

  mod <- object$model

  y <- as_trunc_obs(y)
  keras_y <- k_matrix(y)

  handles_trunc <- object$loss_trunc
  handles_cens <- object$loss_cens
  has_trunc <- any(y$tmin > -Inf | y$tmax < Inf)
  has_cens <- anyNA(y$x)

  if (has_trunc && !handles_trunc) {
    # We can continue here because disregarding the truncation is possible.
    # The fit will just overweight observations with low truncation probability.
    # It may also be correct e.g. if dist has support on [0, Inf) and
    # tmin is always <= 0 and tmax is always Inf.
    warning(
      "`y` seems to contain truncated observations but the model wasn't ",
      "compiled for truncated data.\nResults might be wrong."
    )
  }

  if (has_cens && !handles_cens) {
    bad_obs <- which(is.na(y$x))
    if (length(bad_obs) > 5) {
      bad_obs <- c(bad_obs[1:4], "...")
    }
    # We can stop here because non-censored loss functions will get NaN in their
    # log-likelihood for censored observations, i.e. the gradient would NaN at
    # the first iteration.
    stop(
      "`y` contains censored observations but the model wasn't compiled for ",
      "censored data.\nCensored observations: ", paste(bad_obs, collapse = ", ")
    )
  }

  if (is.null(batch_size) && is.null(steps_per_epoch) && !is_tensorflow_dataset(x)) {
    # This is a better default for reservr_keras_models than keras default of 32.
    batch_size <- min(nrow(y), 10000L)
  }

  fit(
    object = mod,
    x = x,
    y = keras_y,
    batch_size = batch_size,
    epochs = epochs,
    verbose = verbose,
    callbacks = callbacks,
    view_metrics = view_metrics,
    validation_split = validation_split,
    validation_data = validation_data,
    shuffle = shuffle,
    class_weight = class_weight,
    sample_weight = sample_weight,
    initial_epoch = initial_epoch,
    steps_per_epoch = steps_per_epoch,
    validation_steps = validation_steps
  )
}

#' @export
format.reservr_keras_model <- function(x, ...) {
  dist_fmt <- format(x$dist)
  mod_fmt <- reticulate::py_str(x$model)
  cens <- if (x$loss_cens) "enabled" else "disabled"
  trunc <- if (x$loss_trunc) "enabled" else "disabled"
  paste0(
    "A reservr_keras_model (censoring: ", cens, ", truncation: ", trunc, ").\n",
    "Distribution:\n",
    dist_fmt, "\n",
    mod_fmt
  )
}

#' Callback to monitor likelihood gradient components
#'
#' Provides a keras callback to monitor the individual components of the
#' censored and truncated likelihood.
#' Useful for debugging TensorFlow implementations of Distributions.
#'
#' @param object A `reservr_keras_model` created by [tf_compile_model()].
#' @param data Input data for the model.
#' @param obs Observations associated to `data`.
#' @param keep_grads Log actual gradients? (memory hungry!)
#' @param stop_on_na Stop if any likelihood component as NaN in its gradients?
#' @param verbose Print a message if training is halted?
#' The Message will contain information about which likelihood components have
#' NaN in their gradients.
#'
#' @return A `KerasCallback` suitable for passing to [keras3::fit()].
#'
#' @examples
#' dist <- dist_exponential()
#' group <- sample(c(0, 1), size = 100, replace = TRUE)
#' x <- dist$sample(100, with_params = list(rate = group + 1))
#' global_fit <- fit(dist, x)
#'
#' if (interactive() && keras3::is_keras_available()) {
#'   library(keras3)
#'   l_in <- layer_input(shape = 1L)
#'   mod <- tf_compile_model(
#'     inputs = list(l_in),
#'     intermediate_output = l_in,
#'     dist = dist,
#'     optimizer = optimizer_adam(),
#'     censoring = FALSE,
#'     truncation = FALSE
#'   )
#'   tf_initialise_model(mod, global_fit$params)
#'   gradient_tracker <- callback_debug_dist_gradients(mod, k_constant(group), x, keep_grads = TRUE)
#'   fit_history <- fit(
#'     mod,
#'     x = k_constant(group),
#'     y = x,
#'     epochs = 20L,
#'     callbacks = list(
#'       callback_adaptive_lr("loss", factor = 0.5, patience = 2L, verbose = 1L, min_lr = 1.0e-4),
#'       gradient_tracker,
#'       callback_reduce_lr_on_plateau("loss", min_lr = 1.0) # to track lr
#'     )
#'   )
#'   gradient_tracker$gradient_logs[[20]]$dens
#'
#'   plot(fit_history)
#'
#'   predicted_means <- predict(mod, data = k_constant(c(0, 1)))
#' }
#'
#' @export
callback_debug_dist_gradients <- function(object, data, obs,
                                          keep_grads = FALSE,
                                          stop_on_na = TRUE,
                                          verbose = TRUE) {
  assert_that(inherits(object, "reservr_keras_model"),
              msg = "`object` must ba a reservr_keras_model.")
  obs <- as_trunc_obs(obs)
  assert_that(is_bool(keep_grads),
              msg = "`keep_grads` must be a bool.")
  assert_that(is_bool(stop_on_na),
              msg = "`stop_on_na` must be a bool.")
  DebugDistGradientsCallback(
    object = object, data = data, obs = obs,
    keep_grads = keep_grads, stop_on_na = stop_on_na, verbose = verbose
  )
}

DebugDistGradientsCallback <- keras3::Callback(
  "DebugDistGradientsCallback",
  public = list(
    initialize = function(object, data, obs, keep_grads, stop_on_na, verbose) {
      private$.object <- object
      private$.data <- data
      private$.keep_grads <- keep_grads
      private$.stop_on_na <- stop_on_na
      private$.verbose <- verbose
      private$.logd <- object$dist$tf_logdensity()
      if (object$loss_cens || object$loss_trunc) {
        private$.logp <- object$dist$tf_logprobability()
      }
      private$.const <- object$dist$tf_make_constants()
      nobs <- nrow(obs)
      if (!all(is.na(obs$x))) {
        private$.xd <- keras3::as_tensor(ifelse(is.na(obs$x), Inf, obs$x), keras3::config_floatx(), shape = list(nobs))
      }
      if (object$loss_cens && anyNA(obs$x)) {
        private$.xc_lower <- keras3::as_tensor(
          ifelse(is.na(obs$x), obs$xmin, -Inf),
          keras3::config_floatx(),
          shape = list(nobs)
        )
        private$.xc_upper <- keras3::as_tensor(
          ifelse(is.na(obs$x), obs$xmax, Inf),
          keras3::config_floatx(),
          shape = list(nobs)
        )
      }
      if (object$loss_trunc && any(is.finite(obs$tmin) | is.finite(obs$tmax))) {
        private$.xt_lower <- keras3::as_tensor(obs$tmin, keras3::config_floatx(), shape = list(nobs))
        private$.xt_upper <- keras3::as_tensor(obs$tmax, keras3::config_floatx(), shape = list(nobs))
      }
      private$reset()
    },
    on_train_begin = function(logs) {
      private$reset()
    },
    on_epoch_end = function(epoch, logs) {
      private$log_epoch(epoch)
    }
  ),
  private = list(
    .gradient_logs = list(),
    .object = NULL,
    .data = NULL,
    .logd = NULL,
    .logp = NULL,
    .const = NULL,
    .xd = NULL,
    .xc_lower = NULL,
    .xc_upper = NULL,
    .xt_lower = NULL,
    .xt_upper = NULL,
    .keep_grads = FALSE,
    .stop_on_na = TRUE,
    .verbose = TRUE,
    reset = function() {
      private$.gradient_logs <- list()
    },
    log_epoch = function(epoch) {
      `%as%` <- tensorflow::`%as%`
      with(tensorflow::tf$GradientTape(persistent = TRUE) %as% tape, {
        curr_args <- keras3::op_cast(self$model(private$.data), keras3::config_floatx())
        curr_args <- private$.object$output_splitter(curr_args)
        curr_args <- private$.object$output_inflater(curr_args)
        curr_args <- tf_merge_constants(curr_args, private$.const)
        loss_dens <- if (!is.null(private$.xd)) private$.logd(private$.xd, curr_args)
        loss_cens <- if (!is.null(private$.xc_lower))
          private$.logp(private$.xc_lower, private$.xc_upper, curr_args)
        loss_trunc <- if (!is.null(private$.xt_lower))
          private$.logp(private$.xt_lower, private$.xt_upper, curr_args)
      })

      grad_dens <- if (!is.null(loss_dens))
        tape$gradient(loss_dens, self$model$trainable_variables)
      grad_cens <- if (!is.null(loss_cens))
        tape$gradient(loss_cens, self$model$trainable_variables)
      grad_trunc <- if (!is.null(loss_trunc))
        tape$gradient(loss_trunc, self$model$trainable_variables)

      dens_ok <- if (!is.null(grad_dens)) !any(vapply(grad_dens, anyNA, logical(1L))) else TRUE
      cens_ok <- if (!is.null(grad_dens)) !any(vapply(grad_cens, anyNA, logical(1L))) else TRUE
      trunc_ok <- if (!is.null(grad_trunc)) !any(vapply(grad_trunc, anyNA, logical(1L))) else TRUE

      if (!all(dens_ok, cens_ok, trunc_ok) && private$.stop_on_na) {
        self$model$stop_training <- TRUE
        if (private$.verbose) {
          message(sprintf(
            paste(
              "\nEpoch %05d: NaN in Gradients: %s.",
              "Stopped trainig."
            ),
            epoch,
            private$debug_info(grad_dens, grad_cens, grad_trunc)
          ))
        }
      }

      if (!private$.keep_grads) {
        grad_dens <- NULL
        grad_cens <- NULL
        grad_trunc <- NULL
        curr_args <- NULL
      }

      private$.gradient_logs <- c(private$.gradient_logs, list(list(
        epoch = epoch,
        args = curr_args,
        dens = grad_dens,
        dens_ok = dens_ok,
        cens = grad_cens,
        cens_ok = cens_ok,
        trunc = grad_trunc,
        trunc_ok = trunc_ok
      )))
    },
    debug_info = function(grad_dens, grad_cens, grad_trunc) {
      nms <- vapply(
        self$model$trainable_variables, function(var) var$name, character(1L)
      )
      all_ok <- rep_len(FALSE, length(nms))
      bad_dens <- if (!is.null(grad_dens)) vapply(grad_dens, anyNA, logical(1L)) else all_ok
      bad_cens <- if (!is.null(grad_cens)) vapply(grad_cens, anyNA, logical(1L)) else all_ok
      bad_trunc <- if (!is.null(grad_trunc)) vapply(grad_trunc, anyNA, logical(1L)) else all_ok
      any_bad <- bad_dens | bad_cens | bad_trunc
      bad_idx <- which(any_bad)
      fmt <- paste0("%1$s[", paste(
        c(if (any(bad_dens)) "dens: %2$s", if (any(bad_cens)) "cens: %3$s", if (any(bad_trunc)) "trunc: %4$s"),
        collapse = ", "
      ), "]")
      infos <- suppressWarnings(sprintf( # "arguments not used by format" is intentional
        fmt, nms,
        ifelse(bad_dens, "bad", "ok"),
        ifelse(bad_cens, "bad", "ok"),
        ifelse(bad_trunc, "bad", "ok")
      ))[bad_idx]
      paste(infos, collapse = "; ")
    }
  )
)

#' Keras Callback for adaptive learning rate with weight restoration
#'
#' Provides a keras callback similar to [keras::callback_reduce_lr_on_plateau()] but which also restores the weights
#' to the best seen so far whenever a learning rate reduction occurs, and with slightly more restrictive improvement
#' detection.
#'
#' @param monitor quantity to be monitored.
#' @param factor factor by which the learning rate will be reduced. `new_lr = old_lr * factor`.
#' @param patience number of epochs with no significant improvement after which the learning rate will be reduced.
#' @param verbose integer. Set to 1 to receive update messages.
#' @param mode Optimisation mode. "auto" detects the mode from the name of `monitor`. "min" monitors for decreasing
#' metrics. "max" monitors for increasing metrics.
#' @param delta_abs Minimum absolute metric improvement per epoch. The learning rate will be reduced if the average
#' improvement is less than `delta_abs` per epoch for `patience` epochs.
#' @param delta_rel Minimum relative metric improvement per epoch. The learning rate will be reduced if the average
#' improvement is less than `|metric| * delta_rel` per epoch for `patience` epochs.
#' @param cooldown number of epochs to wait before resuming normal operation after learning rate has been reduced.
#' The minimum number of epochs between two learning rate reductions is `patience + cooldown`.
#' @param min_lr lower bound for the learning rate. If a learning rate reduction would lower the learning rate below
#' `min_lr`, it will be clipped at `min_lr` instead and no further reductions will be performed.
#' @param restore_weights Bool. If TRUE, the best weights will be restored at each learning rate reduction.
#' This is very useful if the metric oscillates.
#'
#' @details Note that while [callback_reduce_lr_on_plateau()] automatically logs the learning rate as a metric 'lr',
#' this is currently impossible from R.
#' Thus, if you want to also log the learning rate, you should add [callback_reduce_lr_on_plateau()] with a high
#' `min_lr` to effectively disable the callback but still monitor the learning rate.
#'
#' @return A `KerasCallback` suitable for passing to [keras::fit()].
#'
#' @examples
#' dist <- dist_exponential()
#' group <- sample(c(0, 1), size = 100, replace = TRUE)
#' x <- dist$sample(100, with_params = list(rate = group + 1))
#' global_fit <- fit(dist, x)
#'
#' if (interactive() && keras::is_keras_available()) {
#'   library(keras)
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
#'   fit_history <- fit(
#'     mod,
#'     x = k_constant(group),
#'     y = as_trunc_obs(x),
#'     epochs = 20L,
#'     callbacks = list(
#'       callback_adaptive_lr("loss", factor = 0.5, patience = 2L, verbose = 1L, min_lr = 1.0e-4),
#'       callback_reduce_lr_on_plateau("loss", min_lr = 1.0) # to track lr
#'     )
#'   )
#'
#'   plot(fit_history)
#'
#'   predicted_means <- predict(mod, data = k_constant(c(0, 1)))
#' }
#'
#' @export
callback_adaptive_lr <- function(monitor = "val_loss", factor = 0.1, patience = 10L, verbose = 0L,
                                 mode = c("auto", "min", "max"), delta_abs = 1.0e-4, delta_rel = 0.0, cooldown = 0L,
                                 min_lr = 0, restore_weights = TRUE) {
  mode <- match.arg(mode)
  AdaptiveLRCallback$new(
    monitor = monitor, factor = factor, patience = patience, verbose = verbose, mode = mode,
    delta_abs = delta_abs, delta_rel = delta_rel, cooldown = cooldown, min_lr = min_lr,
    restore_weights = restore_weights
  )
}

AdaptiveLRCallback <- R6Class(
  "AdaptiveLRCallback",
  inherit = keras::KerasCallback,
  public = list(
    initialize = function(monitor = "val_loss", factor = 0.1, patience = 10L, verbose = 0L, mode = "auto",
                          delta_abs = 1.0e-4, delta_rel = 0.0, cooldown = 0L, min_lr = 0,
                          restore_weights = TRUE, ...) {
      self$monitor <- monitor
      self$factor <- factor
      self$patience <- patience
      self$verbose <- verbose
      self$mode <- mode
      self$delta_abs <- delta_abs
      self$delta_rel <- delta_rel
      self$cooldown <- cooldown
      self$min_lr <- min_lr
      self$restore_weights <- restore_weights
      private$reset()
    },
    on_train_begin = function(logs = NULL) {
      private$reset()
      if (private$.restore_weights) {
        private$.best_weights <- self$model$get_weights()
      }
    },
    on_epoch_end = function(epoch, logs = NULL) {
      if (!private$.active) return()

      new_metric <- logs[[self$monitor]]
      if (is.null(new_metric)) {
        warning("Metric '", self$monitor, "' is not available for adaptive learning rate controlling. ",
                "Available metrics are: ", paste(names(logs), collapse = ","))
        return()
      }

      if (private$.cooldown_counter > 0L) {
        private$.cooldown_counter <- private$.cooldown_counter - 1L
        private$.patience_counter <- 0L
      }

      if (private$has_improved(epoch, new_metric)) {
        private$.patience_counter <- 0L
      } else {
        private$.patience_counter <- private$.patience_counter + 1L
      }

      if (private$has_improved_slightly(new_metric)) {
        private$.best_epoch <- epoch
        private$.best_metric <- new_metric
        if (private$.restore_weights) {
          curr_weights <- self$model$get_weights()
          if (any(vapply(curr_weights, anyNA, logical(1L)))) {
            warning(sprintf(
              "\nEpoch %05d: AdaptiveLR got NaN weights for new loss %g, keeping old weights.",
              epoch, new_metric
            ))
          } else {
            private$.best_weights <- curr_weights
          }
        }
      }

      if (private$.patience_counter >= private$.patience || is.na(new_metric)) {
        old_lr <- as.numeric(self$model$optimizer$lr)
        if (old_lr > private$.min_lr) {
          new_lr <- pmax(private$.min_lr, old_lr * private$.factor)
          self$model$optimizer$lr <- new_lr
          private$.patience_counter <- 0L
          private$.cooldown_counter <- private$.cooldown
          if (private$.restore_weights && private$.best_epoch < epoch) {
            self$model$set_weights(private$.best_weights)
            if (private$.verbose) {
              message(sprintf("\nEpoch %05d: AdaptiveLR reducing learning rate to %g and restoring weights from Epoch
               %05d.", epoch, new_lr, private$.best_epoch))
            }
          } else {
            if (private$.verbose) {
              message(sprintf("\nEpoch %05d: AdaptiveLR reducing learning rate to %g.", epoch, new_lr))
            }
          }

          if (new_lr == private$.min_lr) {
            private$.active <- FALSE
            if (private$.verbose) {
              message(sprintf("\nEpoch %05d: AdaptiveLR deactivated. min_lr reached.", epoch))
            }
          }
        }
      }

      # Track learning rate (requires keras >= 2.4.0.9000 to actually work)
      logs[["lr"]] <- self$model$optimizer$lr

      invisible(NULL)
    }
  ),
  private = list(
    .active = TRUE,
    .best_weights = list(),
    .best_metric = NA_real_,
    .best_epoch = NA_integer_,
    .verbose = NA_integer_,
    .patience = NA_integer_,
    .patience_counter = 0L,
    .cooldown = NA_integer_,
    .cooldown_counter = 0L,
    .mode = NA_character_,
    .monitor = NA_character_,
    .restore_weights = NA,
    .factor = NA_real_,
    .min_lr = NA_real_,
    .delta_rel = NA_real_,
    .delta_abs = NA_real_,
    has_improved = NULL,
    has_improved_slightly = NULL,
    reset = function() {
      private$.active <- TRUE
      if (private$.mode == "min") {

        private$has_improved <- function(epoch, new_metric) {
          isTRUE(
            is.infinite(private$.best_metric) ||
              private$.best_metric - new_metric >
                (epoch - private$.best_epoch) * pmax(abs(private$.best_metric) * private$.delta_rel, private$.delta_abs)
          )
        }

        private$has_improved_slightly <- function(new_metric) {
          isTRUE(new_metric < private$.best_metric)
        }

        private$.best_metric <- Inf
      } else {

        private$has_improved <- function(epoch, new_metric) {
          isTRUE(
            is.infinite(private$.best_metric) ||
              new_metric - private$.best_metric >
                (epoch - private$.best_epoch) * pmax(abs(private$.best_metric) * private$.delta_rel, private$.delta_abs)
          )
        }

        private$has_improved_slightly <- function(new_metric) {
          isTRUE(new_metric > private$.best_metric)
        }

        private$.best_metric <- -Inf
      }
      if (private$.restore_weights) {
        private$.best_weights <- list()
      }
      private$.best_epoch <- 0L
      private$.cooldown_counter <- 0L
      private$.patience_counter <- 0L
    }
  ),
  active = list(
    best_weights = function(value) {
      if (!missing(value)) stop("`best_weights` is read-only.")
      private$.best_weights
    },
    best_metric = function(value) {
      if (!missing(value)) stop("`best_metric` is read-only.")
      private$.best_metric
    },
    best_epoch = function(value) {
      if (!missing(value)) stop("`best_epoch` is read-only.")
      private$.best_epoch
    },
    patience = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_integerish(value, finite = TRUE), value > 0L,
                    msg = "`patience` must be a positive integer.")
        private$.patience <- as.integer(value)
      }
      private$.patience
    },
    mode = function(value) {
      if (!missing(value)) {
        value <- match.arg(value, c("min", "max", "auto"))
        if (value == "auto") {
          value <- if (grepl("acc", private$.monitor)) "max" else "min"
        }
        private$.mode <- value
      }
      private$.mode
    },
    monitor = function(value) {
      if (!missing(value)) {
        assert_that(is_string(value), msg = "`monitor` must be a string.")
        private$.monitor <- value
      }
      private$.monitor
    },
    restore_weights = function(value) {
      if (!missing(value)) {
        assert_that(is_bool(value), msg = "`restore_weights` must be a bool.")
        private$.restore_weights <- value
      }
      private$.restore_weights
    },
    factor = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_double(value), value >= 0, value < 1,
                    msg = "`factor` must be a number in [0, 1).")
        private$.factor <- value
      }
      private$.factor
    },
    min_lr = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_double(value), value >= 0,
                    msg = "`min_lr` must be a non-negative number.")
        private$.min_lr <- value
      }
      private$.min_lr
    },
    delta_rel = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_double(value), value >= 0,
                    msg = "`delta_rel` must be a non-negative number.")
        private$.delta_rel <- value
      }
      private$.delta_rel
    },
    delta_abs = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_double(value), value >= 0,
                    msg = "`delta_abs` must be a non-negative number.")
        private$.delta_abs <- value
      }
      private$.delta_abs
    },
    verbose = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_integerish(value, finite = TRUE))
        private$.verbose <- value
      }
      private$.verbose
    },
    cooldown = function(value) {
      if (!missing(value)) {
        assert_that(is_scalar_integerish(value), value >= 0,
                    msg = "`cooldown` must be a non-negative integer.")
        private$.cooldown <- value
      }
      private$.cooldown
    }
  )
)

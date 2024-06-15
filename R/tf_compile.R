#' Compile a Keras model for truncated data under dist
#'
#' @param inputs List of keras input layers
#' @param intermediate_output Intermediate model layer to be used as input to
#' distribution parameters
#' @param dist A `Distribution` to use for compiling the loss and parameter
#' outputs
#' @param censoring A flag, whether the compiled model should support censored
#' observations. Set to `FALSE` for higher efficiency. `fit(...)` will error if
#' the resulting model is used to fit censored observations.
#' @param truncation A flag, whether the compiled model should support truncated
#' observations. Set to `FALSE` for higher efficiency. `fit(...)` will warn if
#' the resuting model is used to fit truncated observations.
#' @inheritParams keras3::compile.keras.src.models.model.Model
#'
#' @return A `reservr_keras_model` that can be used to train truncated
#' and censored observations from `dist` based on input data from `inputs`.
#'
#' @examples
#' dist <- dist_exponential()
#' params <- list(rate = 1.0)
#' N <- 100L
#' rand_input <- runif(N)
#' x <- dist$sample(N, with_params = params)
#'
#' if (interactive()) {
#'   tf_in <- keras3::layer_input(1L)
#'   mod <- tf_compile_model(
#'     inputs = list(tf_in),
#'     intermediate_output = tf_in,
#'     dist = dist,
#'     optimizer = keras3::optimizer_adam(),
#'     censoring = FALSE,
#'     truncation = FALSE
#'   )
#' }
#'
#' @export
tf_compile_model <- function(inputs, intermediate_output, dist, optimizer,
                             censoring = TRUE, truncation = TRUE,
                             metrics = NULL, weighted_metrics = NULL) {
  check_installed(c("tensorflow", "keras3"))

  loss <- if (censoring && truncation) {
    tf_compile_loss_trunc_cens(dist)
  } else if (truncation) {
    tf_compile_loss_trunc(dist)
  } else if (censoring) {
    tf_compile_loss_cens(dist)
  } else {
    tf_compile_loss(dist)
  }

  params <- dist$tf_compile_params(input = intermediate_output)
  params_compressed <- tf_compress_output(params$outputs)

  loss_compressed <- function(obs, output_compressed) {
    outputs_compressed <- params_compressed$output_splitter(output_compressed)
    outputs <- params$output_inflater(outputs_compressed)
    tensorflow::tf$reduce_mean(loss(obs, outputs))
  }

  out <- list(
    dist = dist,
    model = keras3::keras_model(
      inputs = inputs,
      outputs = params_compressed$output
    ),
    output_inflater = params$output_inflater,
    output_splitter = params_compressed$output_splitter,
    loss = loss_compressed,
    loss_cens = censoring,
    loss_trunc = truncation
  )
  keras3::compile(
    out$model,
    optimizer = optimizer,
    loss = loss_compressed,
    metrics = metrics,
    weighted_metrics = weighted_metrics
  )
  class(out) <- "reservr_keras_model"
  out
}

tf_compress_output <- function(outputs) {
  output_names <- names(outputs)
  output_dims <- purrr::map_int(outputs, ~.$shape[[2L]])

  if (length(outputs) > 1L) {
    list(
      output = keras3::layer_concatenate(unname(outputs)),
      output_splitter = eval(bquote(function(output_condensed) {
        output_names <- .(output_names)
        output_dims <- .(output_dims)
        offset <- c(0, cumsum(output_dims))

        out <- vector("list", length = length(output_names))
        names(out) <- output_names

        for (i in seq_along(out)) {
          curr_i <- offset[i] + seq_len(output_dims[i])
          out[[i]] <- output_condensed[, curr_i]
        }

        out
      }))
    )
  } else {
    list(
      output = outputs[[1L]],
      output_splitter = eval(bquote(function(output_condensed) {
        output_name <- .(output_names)
        out <- list(output_condensed)
        names(out) <- output_name
        out
      }))
    )
  }
}

tf_merge_constants <- function(inputs, constants) {
  out <- inputs

  input_keys <- names(inputs) %||% seq_along(inputs)
  const_keys <- names(constants) %||% seq_along(constants)

  missing_vals <- setdiff(const_keys, input_keys)
  out[missing_vals] <- constants[missing_vals]

  merging_vals <- intersect(const_keys, input_keys)
  merging_ok <- vapply(inputs[merging_vals], is.list, logical(1L))
  if (!all(merging_ok)) {
    stop("Can't merge key '", merging_vals[which.min(merging_ok)], "'.")
  }
  out[merging_vals] <- lapply(merging_vals, function(key) {
    tf_merge_constants(inputs[[key]], constants[[key]])
  })
  out
}

maybe_tf_function <- function(f) {
  if (inherits(f, "tensorflow.python.eager.def_function.Function")) {
    f
  } else {
    tensorflow::tf_function(f)
  }
}

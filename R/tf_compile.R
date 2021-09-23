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
#' @inheritParams keras::compile.keras.engine.training.Model
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
#' }
#'
#' @export
tf_compile_model <- function(inputs, intermediate_output, dist, optimizer,
                             censoring = TRUE, truncation = TRUE,
                             metrics = NULL, sample_weight_mode = NULL,
                             weighted_metrics = NULL, target_tensors = NULL) {
  check_installed(c("tensorflow", "keras"))

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
    model = keras::keras_model(
      inputs = inputs,
      outputs = params_compressed$output
    ),
    output_inflater = params$output_inflater,
    output_splitter = params_compressed$output_splitter,
    loss = loss_compressed,
    loss_cens = censoring,
    loss_trunc = truncation
  )
  keras::compile(
    out$model,
    optimizer = optimizer,
    loss = loss_compressed,
    metrics = metrics,
    sample_weight_mode = sample_weight_mode,
    weighted_metrics = weighted_metrics,
    target_tensors = target_tensors
  )
  class(out) <- "reservr_keras_model"
  out
}

tf_compress_output <- function(outputs) {
  output_names <- names(outputs)
  output_dims <- purrr::map_int(outputs, ~.$shape[[2L]])

  if (length(outputs) > 1L) {
    list(
      output = keras::layer_concatenate(unname(outputs)),
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

#' Temporarily disable TensorFlow graph function compilation
#'
#' This can be very useful for debugging tensorflow functions.
#' Errors raised during TensorFlow functions code evaluation do not yield usable
#' line numbers for the R implementation. With graph compilation disabled, the
#' resulting functions will be run as regular R functions instead.
#'
#' @details
#' `tf_with_disable_graph` accomplishes its function by temporarily overwriting [tensorflow::tf_function()] in the
#' package tensorflow. After evaluation, [tensorflow::tf_function()] is automatically restored.
#'
#' @param expr expression to evaluate with disabled graph compilation
#'
#' @return `tf_with_disable_graph` evaluates `expr` with all calls to [tensorflow::tf_function()] transparently ignored.
#'
#' @examples
#' dist <- dist_dirac()
#' # TensorFlow functions are functions of class
#' # c("tensorflow.python.eager.def_function.Function", "python.builtin.object")
#' if (keras::is_keras_available()) {
#'   stopifnot(class(dist$tf_logdensity()) != "function")
#'   stopifnot(class(tf_with_disable_graph(dist$tf_logdensity())) == "function")
#' }
#'
#' @export
#'
#' @importFrom utils assignInNamespace
tf_with_disable_graph <- function(expr) {
  old_tf_function <- tensorflow::tf_function
  assignInNamespace("tf_function", function(f, ...) f, ns = "tensorflow")
  on.exit(assignInNamespace("tf_function", old_tf_function, ns = "tensorflow"))
  eval.parent(expr)
}

#' @rdname tf_with_disable_graph
#'
#' @return `tf_is_graph_disabled` checks if graph compilation is currenty disabled
#'
#' @examples
#' stopifnot(!tf_is_graph_disabled())
#' stopifnot(tf_with_disable_graph(tf_is_graph_disabled()))
#'
#' @export
tf_is_graph_disabled <- function() {
  length(formals(tensorflow::tf_function)) != 4L
}

maybe_tf_function <- function(f) {
  if (inherits(f, "tensorflow.python.eager.def_function.Function")) {
    f
  } else {
    tensorflow::tf_function(f)
  }
}

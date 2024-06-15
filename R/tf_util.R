#' Convert TensorFlow tensors to distribution parameters recursively
#'
#' @param x possibly nested list structure of `tensorflow.tensor`s
#'
#' @return A nested list of vectors suitable as distribution parameters
#'
#' @examples
#' if (interactive()) {
#'   tf_params <- list(
#'     probs = k_matrix(t(c(0.5, 0.3, 0.2))),
#'     shapes = k_matrix(t(c(1L, 2L, 3L)), dtype = "int32"),
#'     scale = keras3::as_tensor(1.0, keras3::config_floatx())
#'   )
#'   params <- as_params(tf_params)
#'   dist <- dist_erlangmix(vector("list", 3L))
#'   dist$sample(10L, with_params = params)
#' }
#'
#' @export
as_params <- function(x) {
  if (is.list(x)) {
    lapply(x, as_params)
  } else if (inherits(x, "tensorflow.tensor")) {
    x <- x$numpy()
    if (length(dim(x)) > 2L) {
      stop(
        "Don't know how to handle tensor with ", length(dim(x)), " dimensions."
      )
    }
    if (length(dim(x)) == 2L) {
      n_cols <- ncol(x)
      if (n_cols == 1L) {
        drop(x)
      } else {
        lapply(
          seq_len(n_cols),
          function(i) x[, i]
        )
      }
    } else {
      x
    }
  } else {
    stop("Don't know how to handle ", class(x)[1L], ".")
  }
}

#' Cast to a TensorFlow matrix
#'
#' @param x Numeric object to be converted to a matrix Tensor.
#' @param dtype Type of the elements of the resulting tensor. Defaults to [keras3::config_floatx()].
#'
#' @return A two-dimensional `tf.Tensor` with values from `x`.
#' The shape will be `(nrow(x), ncol(x))` where `x` is first converted to an R matrix via [as.matrix()].
#'
#' @examples
#' if (interactive()) {
#'   k_matrix(diag(1:3))
#'   k_matrix(diag(1:3), dtype = "int32")
#'   # Vectors are converted to columns:
#'   k_matrix(1:3)
#' }
#'
#' @export
k_matrix <- function(x, dtype = NULL) {
  if (is.null(dtype)) {
    dtype <- keras3::config_floatx()
  }
  keras3::as_tensor(as.matrix(x), dtype = dtype)
}

#' @export
anyNA.tensorflow.tensor <- function(x, recursive = FALSE) {
  check_installed("tensorflow")
  as.logical(tensorflow::tf$math$reduce_any(tensorflow::tf$math$is_nan(x)))
}

#' @export
anyNA.tensorflow.python.framework.indexed_slices.IndexedSlices <- function(x, recursive = FALSE) {
  check_installed("tensorflow")
  as.logical(tensorflow::tf$math$reduce_any(tensorflow::tf$math$is_nan(x)))
}

tf_is_integerish <- function(x) {
  tensorflow::tf$math$mod(x, K$one) == K$zero
}

# Copied from keras to avoid :::
is_tensorflow_dataset <- function(x) {
  inherits(x, "tensorflow.python.data.ops.dataset_ops.DatasetV2") ||
    inherits(x, "tensorflow.python.data.ops.dataset_ops.Dataset")
}

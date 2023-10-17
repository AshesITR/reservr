#' Initialise model weights to a global parameter fit
#'
#' Initialises a compiled `reservr_keras_model` weights such that the predictions are equal to, or close to, the
#' distribution parameters given by `params`.
#'
#' @param model A `reservr_compiled_model` obtained by [tf_compile_model()].
#' @param params A list of distribution parameters compatible with `model`.
#' @param mode An initialisation mode
#' \describe{
#'   \item{scale}{Initialise the biases according to `params` and the kernels uniform on \[-0.1, 0.1] *
#'   bias scale.}
#'   \item{perturb}{Initialise the biases according to `params` and leave the kernels as is.}
#'   \item{zero}{Initialise the biases according to `params` and set the kernel to zero.}
#'   \item{none}{Don't modify the weights.}
#' }
#'
#' @return Invisibly `model` with changed weights
#'
#' @examples
#' dist <- dist_exponential()
#' group <- sample(c(0, 1), size = 100, replace = TRUE)
#' x <- dist$sample(100, with_params = list(rate = group + 1))
#' global_fit <- fit(dist, x)
#'
#' if (keras::is_keras_available()) {
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
#'     x = group,
#'     y = x,
#'     epochs = 200L
#'   )
#'
#'   predicted_means <- predict(mod, data = k_constant(c(0, 1)))
#' }
#'
#' @export
tf_initialise_model <- function(model, params, mode = c("scale", "perturb", "zero", "none")) {
  mode <- match.arg(mode)
  if (mode == "none") {
    return(invisible(model))
  }

  tf_params <- model$dist$tf_make_constants(params)

  init_list <- function(x, prefix = "") {
    if (is.list(x) && length(x) > 0L) {
      nms <- names(x)
      if (is.null(nms)) {
        nms <- as.character(seq_along(x))
      }
      if (nzchar(prefix)) {
        nms <- paste(prefix, nms, sep = "_")
      }
      mapply(init_list, x = x, prefix = nms)
    } else if (inherits(x, "tensorflow.tensor")) {
      layer <- tryCatch(model$model$get_layer(prefix), error = function(e) NULL)
      if (is.null(layer)) return()

      linkfun <- layer$activation[["__name__"]]
      bias <- inverse_linkfun(x, linkfun)
      bias <- tensorflow::tf$reshape(bias, layer$bias$shape)

      layer$bias$assign(bias)

      switch(
        mode,
        scale = {
          new_weights <- tensorflow::tf$random$uniform(
            layer$kernel$shape,
            minval = -0.1, maxval = 0.1, dtype = keras::k_floatx()
          ) * bias[tensorflow::tf$newaxis, ]
          layer$kernel$assign(new_weights)
        },
        perturb = {
          # Leave as-is
        },
        zero = {
          layer$kernel$assign(tensorflow::tf$fill(layer$kernel$shape, K$zero))
        }
      )
    }
  }

  init_list(tf_params)

  invisible(model)
}

inverse_linkfun <- function(tensor, linkfun) {
  switch(
    linkfun,
    softplus = {
      tensor_safe <- tensorflow::tf$where(tensor > 50, 50, tensor)
      tensorflow::tf$where(tensor > 50, tensor, log(exp(tensor_safe) - 1.0))
    },
    softmax = {
      tensor_soft <- tensorflow::tf$where(tensor == 0, 1.0e-7, tensor)
      log(tensor_soft) - log(tensorflow::tf$math$reduce_max(tensor_soft))
    },
    sigmoid = {
      log(tensor) - log(1.0 - tensor)
    },
    linear = {
      tensor
    },
    stop("Unsupported link function '", linkfun, "'.")
  )
}

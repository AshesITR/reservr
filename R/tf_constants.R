#' Keras Constants for frequently used numerical constants
#'
#' Provides active bindings for frequently used numerical constants.
#' Will cache the resulting TensorFlow Tensor objects to avoid re-creation of
#' the same constant objects.
#'
#' @noRd
Constants <- R6Class(
  "Constants",
  active = local({
    consts <- c(
      "neg_inf" = -Inf,
      "inf" = Inf,
      "zero" = 0,
      "neg_one" = -1,
      "one" = 1,
      "one_half" = 0.5,
      "pi" = base::pi,
      "two" = 2,
      "log_sqrt_2pi" = 0.5 * log(2.0 * base::pi),
      "log_2" = log(2.0),
      "sqrt_2" = sqrt(2.0)
    )

    lapply(consts, function(r_value) {
      eval(substitute(function(value) {
        check_installed("keras")
        assert_that(missing(value), msg = "constants are read-only.")
        keras::k_constant(r_value)
      }, list(r_value = r_value)))
    })
  })
)

# Initialized .onLoad with Constants$new()
K <- NULL

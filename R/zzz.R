#' @importFrom assertthat assert_that
#' @importFrom rlang is_string is_null is_integerish eval_tidy enquo
#' is_installed is_scalar_double is_function is_bool is_scalar_integerish
#' is_bare_numeric check_installed %||%
#' @importFrom utils head tail modifyList hasName assignInMyNamespace
#' @importFrom R6 R6Class
#' @importFrom RcppParallel RcppParallelLibs
#' @import stats
NULL

utils::globalVariables(c(
  "self", "super", "private", # R6 Classes
  "curr_lr", "data", "last_loss" # Used by make_tracer from tf_fit
))

#' Textually enumerate strings
#'
#' @param strs A character vector.
#' @param quote Quotes
#'
#' @return A natural enumeration of `strs`, each element being surrounded by
#' `quote`s
#'
#' @examples
#' # Empty string
#' enumerate_strings(character())
#' # Quoted string
#' enumerate_strings("A", quote = "'")
#' # Multiple strings
#' enumerate_strings(LETTERS[1:10], quote = "'")
#'
#' @keywords internal
#' @noRd
enumerate_strings <- function(strs, quote = "") {
  n <- length(strs)
  if (n > 1) {
    paste0(
      paste0(quote, strs[-n], quote, collapse = ", "),
      " and ",
      quote, strs[n], quote
    )
  } else if (n == 1) {
    paste0(quote, strs, quote)
  } else { #> n == 0
    ""
  }
}

map_dbl_matrix <- function(.x, .f, .n, ...) {
  .f <- purrr::as_mapper(.f, ...)
  res <- vapply(.x, .f, numeric(.n), ...)
  dim(res) <- c(.n, length(.x))
  res
}

map2_dbl_matrix <- function(.x, .y, .f, .n, ...) {
  .f <- purrr::as_mapper(.f, ...)
  map_dbl_matrix(
    seq_along(.x),
    function(i, ...) .f(.x[[i]], .y[[i]], ...),
    .n,
    ...
  )
}

map_lgl_matrix <- function(.x, .f, .n, ...) {
  .f <- purrr::as_mapper(.f, ...)
  res <- vapply(.x, .f, logical(.n), ...)
  dim(res) <- c(.n, length(.x))
  res
}

map2_lgl_matrix <- function(.x, .y, .f, .n, ...) {
  .f <- purrr::as_mapper(.f, ...)
  map_lgl_matrix(
    seq_along(.x),
    function(i, ...) .f(.x[[i]], .y[[i]], ...),
    .n,
    ...
  )
}

#' Muffle all NaNs produced warnings
#'
#' @param expr Expression to be evaluated
#'
#' @return The value of `expr` with all Warnings with exact text "NaNs produced"
#' suppressed
#'
#' @keywords internal
#' @noRd
muffle_nans_produced <- function(expr, .envir = parent.frame()) {
  old_lang <- Sys.getenv("LANGUAGE", unset = NA_character_)
  Sys.setenv(LANGUGAGE = "en")
  on.exit({
    if (!is.na(old_lang)) {
      Sys.setenv(LANGUAGE = old_lang)
    } else {
      Sys.unsetenv("LANGUAGE")
    }
  }, add = TRUE)
  muffle_warning(expr, "NaNs produced", .envir = .envir)
}

muffle_warning <- function(expr, text, .envir = parent.frame()) {
  withCallingHandlers(
    eval(expr, .envir),
    warning = function(w) {
      if (identical(w$message, text)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# nocov start
# covr can't get coverage information for .onLoad
.onLoad <- function(libname, pkgname) {
  assignInMyNamespace("K", Constants$new())
}
# nocov end

deep_multiply <- function(x, factor) {
  if (is.list(x)) {
    lapply(x, deep_multiply, factor = factor)
  } else {
    x * factor
  }
}

xcx_inv <- function(yt) {
  vapply(yt, function(y) {
    uniroot(
      function(x) x - cos(x) - y,
      interval = c(-pi / 2, pi / 2)
    )$root
  }, numeric(1L))
}

gradient <- function(func, x, side = NULL) {
  as.numeric(numDeriv::grad(func, x, side = side))
}

jacobian <- function(func, x) {
  as.numeric(numDeriv::jacobian(func, x))
}

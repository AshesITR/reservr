#' @importFrom assertthat assert_that
#' @importFrom rlang is_string is_null is_integerish eval_tidy enquo
#' is_installed is_scalar_double is_function is_bool is_scalar_integerish
#' is_double check_installed %||%
#' @importFrom utils head tail modifyList hasName
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
  K <<- Constants$new()
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

make_interval_union <- function(lowest, highest, include_lowest, include_highest, integer) {
  n <- check_lengths(lowest, highest, include_lowest, include_highest, integer)
  res <- matrix(
    data = c(
      rep_len(lowest, n),
      rep_len(highest, n),
      rep_len(include_lowest, n),
      rep_len(include_highest, n),
      rep_len(integer, n)
    ),
    ncol = 5L,
    dimnames = list(NULL, c("lowest", "highest", "include_lowest", "include_highest", "integer"))
  )
  class(res) <- "interval_union"
  res
}

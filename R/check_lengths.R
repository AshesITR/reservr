#' Check compatibility of lengths
#'
#' @param ... Arguments (unnamed, should be syntactically nice names)
#' @param .len Expected length (defaults to length of longest input in `...`)
#' @param .msg Description of expected length (such as `"nrow(data)"`)
#'
#' @return Invisibly `.len`.
#'
#' @noRd
check_lengths <- function(...,
                          .len = max(lengths(list(...))),
                          .msg = "longest input") {
  nms <- vapply(rlang::exprs(...), purrr::possibly(rlang::expr_name, NA_character_), character(1))
  nms[is.na(nms)] <- sprintf("<argument %d>", seq_along(nms)[is.na(nms)])
  lens <- lengths(list(...))
  len_ok <- lens %in% c(1L, .len)
  if (!all(len_ok)) {
    bad_nms <- nms[!len_ok]
    bad_lens <- lens[!len_ok]

    if (length(bad_nms) == 1L) {
      # One bad length
      stop(sprintf(
        "`%s` must have 1 or %d elements (%s), but has %d.",
        bad_nms, .len, .msg, bad_lens
      ))
    } else {
      # Multiple problematic lengths
      stop(sprintf(
        paste0(
          "Some arguments have problematic lengths: %s, and %s ",
          "must have 1 or %d elements (%s)."
        ),
        paste0(
          "`", head(bad_nms, -1L), "` (", head(bad_lens, -1L), ")",
          collapse = ", "
        ),
        paste0(
          "`", tail(bad_nms,  1L), "` (", tail(bad_lens,  1L), ")"
        ),
        .len,
        .msg
      ))
    }
  }

  invisible(.len)
}

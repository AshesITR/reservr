#' Compute a union of intervals
#'
#' Computes a union of a (possibly nested) list of intervals and points, simplifying as much as possible.
#'
#' @details An `interval_union` is a matrix with columns "lowest", "highest", "include_lowest", "include_highest" and
#' "integer", efficiently representing a list of intervals with ranges (lowest, highest).
#' Use `as.list(union)` to construct a list of `Interval`s representing the same union if needed.
#'
#' @param ... Possibly list-nested intervals, interval_unions or numerics
#'
#' @return An `interval_union`
#' @examples
#' interval_union(interval(0, 1), interval(1, 3))
#' interval_union(interval(0, 1, closed = TRUE), interval(1, 3))
#' interval_union(interval(0, 1, closed = TRUE), interval(2, 3, closed = TRUE))
#' interval_union(interval(0, 1, closed = TRUE, integer = TRUE), interval(2, 3, closed = TRUE, integer = TRUE))
#' interval_union(interval(0, 1, closed = TRUE), interval(0.5, 2, include_lowest = TRUE))
#' interval_union(
#'   interval(0, 1, closed = TRUE, integer = TRUE),
#'   interval(0.5, 2, include_lowest = TRUE, integer = TRUE),
#'   2.5, 3.5, -1.0
#' )
#'
#' # Compute a list of atomic intervals necessary to represent the interval union
#' as.list(
#'   interval_union(interval(0, 1, closed = TRUE, integer = TRUE), interval(0.5, 2, closed = TRUE, integer = TRUE))
#' )
interval_union <- function(...) {
  # Unnest all lists given to ..., but make sure integer-only arguments stay as a list.
  args <- as.list(unlist(list(...)))

  is_interval <- vapply(args, is.Interval, logical(1L))
  is_numeric <- vapply(args, is.numeric, logical(1L))
  is_interval_union <- vapply(args, inherits, logical(1L), what = "interval_union")

  if (!all(is_interval | is_numeric | is_interval_union)) {
    stop("`interval_union()` only accepts Intervals, interval unions and numerics.")
  }

  n <- length(args) - sum(is_interval_union)
  interval_data <- matrix(
    data = NA_real_, nrow = n, ncol = 5,
    dimnames = list(NULL, c("lowest", "highest", "include_lowest", "include_highest", "integer"))
  )

  # Constants are represented as [a, a] (int), i.e. c(a, a, 1, 1, 1)
  interval_data[is_numeric, ] <- c(
    rep(unlist(args[is_numeric]), 2L),
    rep(1, 3L * sum(is_numeric))
  )
  interval_data[is_interval, ] <- t(vapply(args[is_interval], function(int) {
    c(int$range, int$include_lowest, int$include_highest, int$integer)
  }, numeric(5L)))

  if (any(is_interval_union)) {
    interval_data <- do.call(rbind, c(list(interval_data), args[is_interval_union]))
  }

  interval_data <- simplify_interval_data(interval_data)
  row.names(interval_data) <- NULL
  class(interval_data) <- "interval_union"

  interval_data
}

simplify_interval_data <- function(interval_data) {
  # 1. deduplicate
  interval_data <- unique(interval_data)

  # 2. sort
  interval_data <- interval_data[order(
    interval_data[, "integer"],
    interval_data[, "lowest"],
    interval_data[, "highest"],
    interval_data[, "include_lowest"],
    interval_data[, "include_highest"]
  ), , drop = FALSE]

  # 3. collapse contained intervals
  n_discrete <- sum(interval_data[, "integer"])
  n_continuous <- nrow(interval_data) - n_discrete

  find_contained <- function(x, interval_list, side = "inner", integer = FALSE) {
    low_ok <- interval_list[, "lowest"] < x |
      (interval_list[, "include_lowest"] | side == "left") & interval_list[, "lowest"] == x
    high_ok <- interval_list[, "highest"] > x |
      (interval_list[, "include_highest"] | side == "right") & interval_list[, "highest"] == x
    if (integer) {
      int_ok <- x - interval_list[, "lowest"] == trunc(x - interval_list[, "lowest"])
      which(low_ok & high_ok & int_ok)
    } else {
      which(low_ok & high_ok)
    }
  }
  sequential_merge <- function(interval_list, integer) {
    keep <- interval_list[1L, , drop = FALSE]
    for (i in seq_len(nrow(interval_list) - 1L)) {
      curr <- interval_list[i + 1L, ]
      low_candidates <- find_contained(
        curr[["lowest"]], keep,
        side = if (curr[["include_lowest"]]) "inner" else "left", integer = integer
      )
      high_candidates <- find_contained(
        curr[["highest"]], keep,
        side = if (curr[["include_highest"]]) "inner" else "right", integer = integer
      )
      if (length(high_candidates) == 0L && length(low_candidates) == 0L) {
        # disjoint and non-overlapping
        keep <- rbind(keep, interval_data[i + 1L, , drop = FALSE])
      } else if (length(intersect(low_candidates, high_candidates)) > 0L) {
        # not needed, fully contained in another interval
      } else if (length(setdiff(high_candidates, low_candidates)) > 0L) {
        # overlapping to the right
        i_extend <- min(setdiff(high_candidates, low_candidates))
        keep[i_extend, c("lowest", "include_lowest")] <- curr[c("lowest", "include_lowest")]
      } else if (length(setdiff(low_candidates, high_candidates)) > 0L) {
        # overlapping to the left
        i_extend <- min(setdiff(low_candidates, high_candidates))
        keep[i_extend, c("highest", "include_highest")] <- curr[c("highest", "include_highest")]
      }
    }

    if (integer) {
      # 1. Transform all integer intervals to [a, b) for simple merging
      keep[keep[, "include_lowest"] == 0.0, "lowest"] <- keep[keep[, "include_lowest"] == 0.0, "lowest"] + 1.0
      keep[, "include_lowest"] <- 1.0

      # make sure highest is actually an integerish endpoint (e.g. [
      interval_widths <- trunc(keep[, "highest", drop = TRUE] - keep[, "lowest", drop = TRUE])
      new_highest <- keep[, "lowest", drop = TRUE] + interval_widths
      # new_highest <= keep[, "highest"]. If <, include_highest must be set to 1.0
      # e.g. [0, 1.5) -> [0, 1]
      if (any(new_highest < keep[, "highest"])) {
        keep[new_highest < keep[, "highest"], "indclude_highest"] <- 1.0
        keep[, "highest"] <- new_highest
      }

      # [0, 1] -> [0, 2)
      keep[keep[, "include_highest"] == 1.0, "highest"] <- keep[keep[, "include_highest"] == 1.0, "highest"] + 1.0
      keep[, "include_highest"] <- 0.0


      # Merge [a, b) [b, c) -> [a, c); use equality comparison
      for (i in rev(seq_len(nrow(keep)))) {
        curr_lowest <- keep[i, "lowest"]
        if (curr_lowest %in% keep[, "highest"]) {
          i_merge <- min(which(curr_lowest == keep[, "highest"]))
          keep[i_merge, "highest"] <- keep[i, "highest"]
          keep <- keep[-i, , drop = FALSE]
        }
      }

      # 2. Transform all integer intervals to [a, b - 1] for easier reading
      keep[, "highest"] <- keep[, "highest"] - 1.0
      keep[, "include_highest"] <- 1.0
    }

    keep
  }

  if (n_continuous > 1L) {
    keep_continuous <- sequential_merge(head(interval_data, n_continuous), integer = FALSE)
    interval_data <- rbind(keep_continuous, tail(interval_data, n_discrete))
    n_continuous <- nrow(keep_continuous)
  }

  if (n_discrete > 1L) {
    keep_discrete <- sequential_merge(tail(interval_data, n_discrete), integer = TRUE)
    interval_data <- rbind(head(interval_data, n_continuous), keep_discrete)
  }

  interval_data
}

#' @export
as.list.interval_union <- function(interval_union, ...) {
  if (nrow(interval_union) == 0L) return(list())
  apply(interval_union, 1L, function(row) {
    interval(
      range = c(row[["lowest"]], row[["highest"]]),
      include_lowest = as.logical(row[["include_lowest"]]),
      include_highest = as.logical(row[["include_highest"]]),
      integer = as.logical(row[["integer"]])
    )
  }, simplify = FALSE)
}

#' @export
format.interval_union <- function(interval_union, ...) {
  ilst <- as.list(interval_union)
  fmts <- vapply(ilst, format, character(1L))
  paste0("{", paste(fmts, collapse = ", "), "}")
}

#' @export
print.interval_union <- function(interval_union, ...) {
  cat(format(interval_union), "\n", sep = "")
}

interval_union <- function(...) {
  # Unnest all lists given to ..., but make sure integer-only arguments stay as a list.
  args <- as.list(unlist(list(...)))

  is_interval <- vapply(args, is.Interval, logical(1L))
  is_numeric <- vapply(args, is.numeric, logical(1L))

  if (!all(is_interval | is_numeric)) {
    stop("`interval_union()` only accepts Intervals and numerics.")
  }

  n <- length(args)
  interval_data <- matrix(
    data = NA_real_, nrow = n, ncol = 5,
    dimnames = list(NULL, c("lowest", "highest", "include_lowest", "include_highest", "integer"))
  )

  interval_data[is_numeric, ] <- c(
    rep(unlist(args[is_numeric]), 2L),
    rep(c(1, 1, 1), each = sum(is_numeric))
  )
  interval_data[is_interval, ] <- t(vapply(args[is_interval], function(int) {
    c(int$range, int$include_lowest, int$include_highest, int$integer)
  }, numeric(5L)))

  interval_data <- simplify_interval_data(interval_data)
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
  ), ]

  # 3. collapse contained intervals
  n_continuous <- nrow(interval_data) - sum(interval_data[, "integer"])
  find_contained <- function(x, interval_list, side = "inner") {
    low_ok <- interval_list[, "lowest"] < x |
      (interval_list[, "include_lowest"] | side == "left") & interval_list[, "lowest"] == x
    high_ok <- interval_list[, "highest"] > x |
      (interval_list[, "include_highest"] | side == "right") & interval_list[, "highest"] == x
    which(low_ok & high_ok)
  }
  if (n_continuous > 1L) {
    keep_continuous <- interval_data[1L, , drop = FALSE]
    for (i in seq_len(n_continuous - 1L)) {
      curr <- interval_data[i + 1L, ]
      low_candidates <- find_contained(
        curr[["lowest"]], keep_continuous, side = if (curr[["include_lowest"]]) "inner" else "left"
      )
      high_candidates <- find_contained(
        curr[["highest"]], keep_continuous, side = if (curr[["include_highest"]]) "inner" else "right"
      )
      if (length(high_candidates) == 0L && length(low_candidates) == 0L) {
        # disjoint and non-overlapping
        keep_continuous <- rbind(keep_continuous, interval_data[i + 1L, , drop = FALSE])
      } else if (length(low_candidates) == 0L && length(high_candidates) > 0L) {
        # overlapping to the right
        i_extend <- min(high_candidates)
        keep_continuous[i_extend, c("lowest", "include_lowest")] <- curr[c("lowest", "include_lowest")]
      } else if (length(high_candidates) == 0L && length(low_candidates) > 0L) {
        # overlapping to the left
        i_extend <- min(low_candidates)
        keep_continuous[i_extend, c("highest", "include_highest")] <- curr[c("highest", "include_highest")]
      } else if (length(intersect(low_candidates, high_candidates)) == 0L) {
        # not needed (contained completely within another)
      }
    }
  }

  interval_data
}

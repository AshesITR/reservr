Interval <- R6Class(
  classname = "Interval",
  public = list(
    initialize = function(range, include_lowest, include_highest, integer,
                          read_only) {
      self$range <- range
      self$include_lowest <- include_lowest
      self$include_highest <- include_highest
      self$integer <- integer
      if (read_only) private$.read_only <- TRUE
    },
    format = function() {
      if (self$range[1] < self$range[2]) {
        paste0(
          if (self$include_lowest) "[" else "(",
          self$range[1],
          ", ",
          self$range[2],
          if (self$include_highest) "]" else ")",
          if (self$integer) " (int)" else NULL
        )
      } else if (self$is_empty()) {
        "{}"
      } else {
        paste0("{", self$range[1], "}")
      }
    },
    print = function() {
      cat(self$format(), "\n", sep = "")
    },
    contains = function(x) {
      if (is.Interval(x)) {
        int_ok <- x$integer | !self$integer
        lo_ok <- x$range[1] > self$range[1] ||
          (x$range[1] == self$range[1] &&
            (!x$include_lowest || self$include_lowest))
        hi_ok <- x$range[2] < self$range[2] ||
          (x$range[2] == self$range[2] &&
            (!x$include_highest || self$include_highest))
        int_ok && lo_ok && hi_ok
      } else if (is.numeric(x)) {
        res <- logical(length(x))
        res[x > self$range[1] & x < self$range[2]] <- TRUE
        if (self$include_lowest) {
          res[x == self$range[1]] <- TRUE
        }
        if (self$include_highest) {
          res[x == self$range[2]] <- TRUE
        }
        if (self$integer) {
          res[(x - self$range[1]) != trunc(x - self$range[1])] <- FALSE
        }
        res
      } else {
        stop("`x` must either be an Interval or a numeric vector.")
      }
    },
    is_empty = function() {
      self$range[1] == self$range[2] &&
        !(self$include_lowest && self$include_highest)
    },
    equals = function(b) {
      if (!is.Interval(b)) return(FALSE)
      identical(self$range, b$range) &&
        identical(self$include_lowest, b$include_lowest) &&
        identical(self$include_highest, b$include_highest) &&
        identical(self$integer, b$integer)
    },
    tf_make_layer = function(input, name = NULL, size = 1L) {
      check_installed("keras")

      # TODO support closed ranges as well

      if (self$integer) {
        stop(
          "Unsupported range ", self$format(), " for parameter ",
          name, ". Integer ranges are not supported."
        )
      } else if (I_REALS$equals(self)) { # (-Inf, Inf)
        int_case <- "r"
      } else if (I_POSITIVE_REALS$equals(self)) { # (0, Inf)
        int_case <- "r_plus"
      } else if (I_UNIT_INTERVAL$equals(self)) { # (0, 1)
        int_case <- "unit"
      } else if (all(is.finite(self$range))) { # (a, b)
        int_case <- "interval"
      } else { # (-Inf, a) or (a, Inf)
        int_case <- "half_line"
      }

      if (int_case %in% c("r", "r_plus", "unit")) {
        activation <- switch(
          int_case,
          r = "linear",
          r_plus = "softplus",
          unit = "sigmoid"
        )

        keras::layer_dense(
          object = input,
          units = size,
          activation = activation,
          name = name
        )
      } else if (int_case == "interval") {
        # (0, 1) -> (a, b)
        multip <- diag(self$range[2L] - self$range[1L], size, size)
        const <- rep_len(self$range[1L], size)

        keras::layer_dense(
          object = input,
          units = size,
          activation = "sigmoid",
          name = name
        ) %>%
          keras::layer_dense(
            units = size, activation = "linear",
            weights = list(
              keras::k_constant(multip, shape = c(size, size)),
              keras::k_constant(const, shape = size)
            ), trainable = FALSE
          )
      } else { # half_line
        if (self$range[1L] == -Inf) {
          # (0, Inf) -> (-Inf, a)
          multip <- diag(-1.0, size, size)
          const <- rep_len(self$range[2L], size)
        } else {
          # (0, Inf) -> (a, Inf)
          multip <- diag(1.0, size, size)
          const <- rep_len(self$range[1L], size)
        }

        keras::layer_dense(
          object = input,
          units = size,
          activation = "softplus",
          name = name
        ) %>%
          keras::layer_dense(
            units = size, activation = "linear",
            weights = list(
              keras::k_constant(multip, shape = c(size, size)),
              keras::k_constant(const, shape = size)
            ),
            trainable = FALSE
          )
      }
    }
  ),
  private = list(
    .range = c(-Inf, Inf),
    .include_lowest = FALSE,
    .include_highest = FALSE,
    .integer = FALSE,
    .read_only = FALSE,
    check_write = function() {
      assert_that(
        !private$.read_only,
        msg = "This interval is read-only and cannot be changed."
      )
    }
  ),
  active = list(
    range = function(value) {
      if (missing(value)) {
        private$.range
      } else {
        private$check_write()
        assert_that(
          is_double(value, n = 2L),
          value[1] <= value[2],
          msg = "`range` must be a sorted vector of two numbers."
        )
        private$.range <- value
      }
    },
    include_highest = function(value) {
      if (missing(value)) {
        private$.include_highest
      } else {
        private$check_write()
        assert_that(is_bool(value), msg = "`include_highest` must be a bool.")
        private$.include_highest <- value
      }
    },
    include_lowest = function(value) {
      if (missing(value)) {
        private$.include_lowest
      } else {
        private$check_write()
        assert_that(is_bool(value), msg = "`include_lowest` must be a bool.")
        private$.include_lowest <- value
      }
    },
    integer = function(value) {
      if (missing(value)) {
        private$.integer
      } else {
        private$check_write()
        assert_that(is_bool(value), msg = "`integer` must be a bool.")
        if (value) {
          rng <- private$.range
          rng[2L] <- rng[1L] + trunc(rng[2L] - rng[1L])
          if (rng[2L] != private$.range[2L]) private$.include_highest <- TRUE
          private$.range <- rng
        }
        private$.integer <- value
      }
    }
  )
)

#' Intervals
#'
#' @param range The interval boundaries as a sorted two-element numeric vector.
#' @param ... First argument is used as the endpoint if `range` has length 1.
#' Additional arguments, or any if `range` has length 2, cause a warning and
#' will be ignored.
#' @param include_lowest Is the lower boundary part of the interval?
#' @param include_highest Is the upper boundary part of the interval?
#' @param closed Is the interval closed?
#' @param integer Is the interval only over the integers?
#' @param read_only Make the interval object read-only?
#'
#' @return `interval` returns an `Interval`.
#' `is.Interval` returns `TRUE` if `x` is an `Interval`, `FALSE` otherwise.
#'
#' @seealso interval-operations
#' @export
#'
#' @examples
#' # The real line
#' interval()
#'
#' # Closed unit interval
#' interval(c(0, 1), closed = TRUE)
#' # Alternative form
#' interval(0, 1, closed = TRUE)
#'
#' # Non-negative real line
#' interval(c(0, Inf), include_lowest = TRUE)
interval <- function(
  range = c(-Inf, Inf), ..., include_lowest = closed, include_highest = closed,
  closed = FALSE, integer = FALSE, read_only = FALSE
) {
  stopifnot("`range` must be numeric." = is.numeric(range))

  if (length(range) == 2L) {
    if (...length() > 0L) {
      warning(...length(), " dot arguments provided but range is length 2.")
    }
  } else if (length(range) == 1L &&
    ...length() > 0L &&
    length(..1) == 1L) {
    if (...length() > 1L) {
      warning(...length(), " dot arguments provided, only using the first.")
    }

    stopifnot("upper bound (..1) must be numeric." = is.numeric(..1))
    range <- c(range, ..1)
  } else {
    stop(
      "Invalid arguments. Provide either a two-element numeric vector or",
      " two scalar numerics."
    )
  }

  Interval$new(
    range = range,
    include_lowest = include_lowest,
    include_highest = include_highest,
    integer = integer,
    read_only = read_only
  )
}

#' Convex union and intersection of intervals
#'
#' @param intervals A list of `Interval`s.
#' @param ... appened to `intervals` if present.
#'
#' @return
#' `interval_convex_union` returns the convex union of all intervals in `intervals`.
#' This is the smallest interval completely containing all intervals.
#'
#' @seealso [interval_union] for constructing non-convex interval unions.
#'
#' @export
#'
#' @examples
#' interval_convex_union(
#'   interval(c(0, 1), closed = TRUE),
#'   interval(c(1, 2))
#' )
#'
#' interval_convex_union(
#'   interval(c(0, 5)),
#'   interval(c(1, 4), closed = TRUE)
#' )
#'
#' # Convex union is not equal to set union:
#' interval_convex_union(
#'   interval(c(0, 1)),
#'   interval(c(2, 3))
#' )
#'
#' # The empty union is {}
#' interval_convex_union()
#' @name interval-operations
#' @seealso interval
interval_convex_union <- function(..., intervals = list()) {
  # FIXME consider integer intervals
  if (rlang::dots_n(...) > 0)
    intervals <- c(intervals, list(...))
  assert_that(
    is.list(intervals),
    all(vapply(intervals, is.Interval, logical(1))),
    msg = "`intervals` must be a list of intervals."
  )
  if (length(intervals) == 0L) return(interval(c(0, 0)))

  dat <- vapply(
    intervals,
    function(i) c(i$range, i$include_lowest, i$include_highest),
    numeric(4)
  )

  range <- c(min(dat[1, ]), max(dat[2, ]))
  include_lowest <- sum(dat[3, dat[1, ] == range[1]]) > 0
  include_highest <- sum(dat[3, dat[2, ] == range[2]]) > 0

  interval(
    range = range,
    include_lowest = include_lowest,
    include_highest = include_highest
  )
}

#' @rdname interval-operations
#' @return
#' `interval_intersection` returns the set intersection of all intervals in
#' `intervals`. The empty set is represented by the open interval (0, 0).
#'
#' @export
#'
#' @examples
#'
#' interval_intersection(
#'   interval(c(0, 1)),
#'   interval(c(0.5, 2))
#' )
#'
#' interval_intersection(
#'   interval(c(0, Inf)),
#'   interval(c(-Inf, 0))
#' )
#'
#' interval_intersection(
#'   interval(c(0, Inf), include_lowest = TRUE),
#'   interval(c(-Inf, 0), include_highest = TRUE)
#' )
#'
#' interval_intersection(
#'   interval(c(0, 5)),
#'   interval(c(1, 6), closed = TRUE)
#' )
#'
#' # The empty intersection is (-Inf, Inf)
#' interval_intersection()
interval_intersection <- function(..., intervals = list()) {
  # FIXME consider integer intervals
  if (rlang::dots_n(...) > 0)
    intervals <- c(intervals, list(...))
  assert_that(
    is.list(intervals),
    all(vapply(intervals, is.Interval, logical(1))),
    msg = "`intervals` must be a list of intervals."
  )
  if (length(intervals) == 0L) return(interval())

  dat <- vapply(
    intervals,
    function(i) c(i$range, i$include_lowest, i$include_highest),
    numeric(4)
  )

  range <- c(max(dat[1, ]), min(dat[2, ]))
  include_lowest <- sum(1 - dat[3, dat[1, ] == range[1]]) == 0
  include_highest <- sum(1 - dat[4, dat[2, ] == range[2]]) == 0

  if (range[1] > range[2]) {
    range <- c(0, 0)
    include_lowest <- FALSE
    include_highest <- FALSE
  }

  interval(
    range = range,
    include_lowest = include_lowest,
    include_highest = include_highest
  )
}

#' @rdname interval
#' @param x An object.
#' @export
is.Interval <- function(x) {
  inherits(x, "Interval")
}

I_REALS <- interval(read_only = TRUE)
I_POSITIVE_REALS <- interval(0, Inf, read_only = TRUE)
I_UNIT_INTERVAL <- interval(0, 1, closed = TRUE, read_only = TRUE)
I_NATURALS <- interval(0, Inf,
                       include_lowest = TRUE, integer = TRUE, read_only = TRUE)

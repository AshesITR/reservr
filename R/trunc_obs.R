#' Define a set of truncated observations
#'
#' If `x` is missing, both `xmin` and `xmax` must be specified.
#'
#' Uncensored observations must satisfy `tmin <= xmin = x = xmax <= tmax`.
#' Censored observations must satisfy `tmin <= xmin < xmax <= tmax` and `x = NA`.
#'
#' @param x Observations
#' @param xmin,xmax Censoring bounds. If `xmin != xmax`, `x` must be `NA`.
#' @param tmin,tmax Truncation bounds. May vary per observation.
#' @param w Case weights
#'
#' @return **trunc_obs**: A `trunc_obs` tibble with columns `x`, `xmin`, `xmax`,
#' `tmin` and `tmax` describing possibly interval-censored observations with
#' truncation
#' @export
#'
#' @examples
#' N <- 100
#' x <- rexp(N, 0.5)
#'
#' # Random, observation dependent truncation intervals
#' tmin <- runif(N, 0, 1)
#' tmax <- tmin + runif(N, 1, 2)
#'
#' oob <- x < tmin | x > tmax
#' x <- x[!oob]
#' tmin <- tmin[!oob]
#' tmax <- tmax[!oob]
#'
#' # Number of observations after truncation
#' N <- length(x)
#'
#' # Randomly interval censor 30% of observations
#' cens <- rbinom(N, 1, 0.3) == 1L
#' xmin <- x
#' xmax <- x
#' xmin[cens] <- pmax(tmin[cens], floor(x[cens]))
#' xmax[cens] <- pmin(tmax[cens], ceiling(x[cens]))
#' x[cens] <- NA
#'
#' trunc_obs(x, xmin, xmax, tmin, tmax)
#'
trunc_obs <- function(x, xmin = x, xmax = x, tmin = -Inf, tmax = Inf, w = 1.0) {
  if (missing(x) && !missing(xmin) && !missing(xmax)) {
    n <- max(length(xmin), length(xmax))
  } else {
    n <- length(x)
  }
  check_lengths(xmin, max, tmin, tmax, w,
                .len = n, .msg = "number of observations")
  xmin <- rep_len(xmin, n)
  xmax <- rep_len(xmax, n)
  tmin <- rep_len(tmin, n)
  tmax <- rep_len(tmax, n)
  w <- rep_len(w, n)
  if (missing(x)) {
    x <- ifelse(xmin == xmax, xmin, NA_real_)
  }

  check_trunc(x, xmin, xmax, tmin, tmax)
  res <- data.frame(x = x, xmin = xmin, xmax = xmax,
                    tmin = tmin, tmax = tmax, w = w)
  class(res) <- c("trunc_obs", "tbl_df", "tbl", "data.frame")

  res
}

#' @rdname trunc_obs
#'
#' @param .data A data frame or numeric vector.
#'
#' @return `as_trunc_obs` returns a `trunc_obs` tibble.
#' @export
#'
#' @examples
#' as_trunc_obs(c(1, 2, 3))
#' as_trunc_obs(data.frame(x = 1:3, tmin = 0, tmax = 10))
#' as_trunc_obs(data.frame(x = c(1, NA), xmin = c(1, 2), xmax = c(1, 3)))
as_trunc_obs <- function(.data) {
  if (inherits(.data, "trunc_obs")) {
    return(.data)
  }

  if (is.numeric(.data)) {
    return(trunc_obs(x = .data))
  }

  assert_that(
    is.data.frame(.data),
    msg = "`.data` must either be numeric or a data frame."
  )

  has_obs <- "x" %in% colnames(.data)

  has_cens <- "xmin" %in% colnames(.data)
  assert_that(
    "xmin" %in% colnames(.data) == "xmax" %in% colnames(.data),
    msg = "`xmin` and `xmax` must occur together."
  )

  assert_that(
    has_obs | has_cens,
    msg = "`.data` must have at least `x` or (`xmin`, `xmax`) as columns."
  )


  if (has_obs) {
    x <- .data$x
    xmin <- if (has_cens) .data$xmin else x
    xmax <- if (has_cens) .data$xmax else x
  } else {
    # Auto-comupte x from xmin and xmax
    xmin <- .data$xmin
    xmax <- .data$xmax
    x <- xmin
    x[xmin != xmax] <- NA
  }

  tmin <- if ("tmin" %in% colnames(.data)) .data$tmin else -Inf
  tmax <- if ("tmax" %in% colnames(.data)) .data$tmax else Inf

  w <- if ("w" %in% colnames(.data)) .data$w else 1.0

  trunc_obs(x, xmin, xmax, tmin, tmax, w)
}

#' @rdname trunc_obs
#'
#' @param .data A data frame or numeric vector.
#' @param tmin_new New truncation minimum
#' @param tmax_new New truncation maximum
#' @param .partial Enable partial truncation of censored observations?
#' This could potentially create inconsistent data if the actual observation
#' lies outside of the truncation bounds but the censoring interval overlaps.
#'
#' @return `truncate_obs` returns a `trunc_obs` tibble with possibly fewer
#' observations than `.data` and updated truncation bounds.
#'
#' @export
#' @examples
#' truncate_obs(1:10, tmin_new = 2.0, tmax_new = 8.0)
truncate_obs <- function(.data, tmin_new = -Inf, tmax_new = Inf,
                         .partial = FALSE) {
  obs <- as_trunc_obs(.data)

  check_lengths(
    tmin_new, tmax_new,
    .len = nrow(obs), .msg = "number of observations"
  )

  tmin_new <- rep_len(tmin_new, nrow(obs))
  tmax_new <- rep_len(tmax_new, nrow(obs))

  # Interval-censored observations are truncated if they intersect with the
  # out-of-bounds range.
  oob <- if (.partial) {
    obs$xmax < tmin_new | obs$xmin > tmax_new
  } else {
    obs$xmin < tmin_new | obs$xmax > tmax_new
  }
  res <- obs[!oob, ]

  res$tmin <- pmax(res$tmin, tmin_new[!oob])
  res$tmax <- pmin(res$tmax, tmax_new[!oob])
  if (.partial) {
    res$xmin <- pmax(res$xmin, res$tmin)
    res$xmax <- pmin(res$xmax, res$tmax)
  }

  res
}

#' @rdname trunc_obs
#'
#' @param accident accident time (unquoted, evaluated in `.data`)
#' @param delay reporting delay (unquoted, evaluated in `.data`)
#' @param time evaluation time (unquoted, evaluated in `.data`)
#' @param .truncate Should claims reported after `time` be silently discarded?
#' If there are claims reported after `time` and `.truncate` is `FALSE`,
#' an error will be raised.
#'
#' @return `repdel_obs` returns a `trunc_obs` tibble corresponding to the reporting
#' delay observations of each claim. If `.truncate` is `FALSE`, the result is
#' guaranteed to have the same number of rows as `.data`.
#'
#' @export
repdel_obs <- function(.data, accident, delay, time, .truncate = FALSE) {
  accident <- eval_tidy(enquo(accident), data = .data)
  delay <- eval_tidy(enquo(delay), data = .data)
  time <- eval_tidy(enquo(time), data = .data)

  assert_that(
    is.numeric(accident),
    is.numeric(delay),
    is.numeric(time),
    msg = "`accident`, `delay` and `time` must be numeric."
  )
  assert_that(
    all(delay >= 0.0),
    msg = "`delay` must be non-negative."
  )

  n <- check_lengths(accident, delay, time, .len = nrow(.data),
                     .msg = "nrow(.data)")

  accident <- rep_len(accident, n)
  delay <- rep_len(delay, n)
  time <- rep_len(time, n)

  obs <- accident + delay <= time
  if (!all(obs)) {
    if (!.truncate) {
      stop("`.truncate` must be TRUE if there are observations with ",
           "`accident` + `delay` > `time`.")
    }

    accident <- accident[obs]
    delay <- delay[obs]
    time <- time[obs]
    n <- sum(obs)
  }

  res <- data.frame(x = delay, xmin = delay, xmax = delay,
                    tmin = 0.0, tmax = time - accident, w = rep_len(1.0, n))
  class(res) <- c("trunc_obs", "tbl_df", "tbl", "data.frame")

  res
}

check_trunc <- function(x, xmin, xmax, tmin, tmax, .min = -Inf, .max = Inf) {
  cens <- is.na(x)
  in_range <- function(v, lower = .min, upper = .max,
                       include = FALSE, na_rm = FALSE) {
    if (include) {
      res <- (lower <= v & v <= upper) | (na_rm & is.na(v))
    } else {
      res <- (lower < v & v < upper) | (na_rm & is.na(v))
    }
    isTRUE(all(res))
  }

  assert_that(
    in_range(x, na_rm = TRUE),
    in_range(xmin, include = TRUE),
    in_range(xmax, include = TRUE),
    in_range(tmin, include = TRUE),
    in_range(tmax, include = TRUE),
    msg = glue::glue("`x` must be in ({.min}, {.max}); ",
                     "`xmin`, `xmax`, `tmin` and `tmax` ",
                     "must be in [{.min}, {.max}].")
  )

  assert_that(
    all(tmin < tmax),
    msg = "`tmin` < `tmax` must always hold."
  )

  assert_that(
    in_range(x, tmin, tmax, na_rm = TRUE, include = TRUE),
    in_range(xmin, tmin, tmax, include = TRUE),
    in_range(xmax, tmin, tmax, include = TRUE),
    msg = "`x`, `xmin` and `xmax` must be in [`tmin`, `tmax`]."
  )

  assert_that(
    all(xmin[cens] < xmax[cens]),
    msg = "Censored values must have `x` = NA and `xmin` < `xmax`."
  )

  assert_that(
    all(xmin[!cens] == x[!cens]),
    all(xmax[!cens] == x[!cens]),
    msg = "Non-censored values must have `x` = `xmin` = `xmax`."
  )

  invisible(TRUE)
}

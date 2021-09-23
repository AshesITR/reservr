#' Truncate claims data subject to reporting delay
#'
#' @param data Full claims data including IBNR
#' @param accident Accident times. May be an unquoted column name from `data`.
#' @param delay Reporting delays. May be an unquoted column name from `data`.
#' @param time Observation time (scalar number or one per claim).
#' Claims with `accident + delay > time` will be truncated.
#' Set `time = Inf` to only compute reporting times and perform no truncation.
#' @param .report_col `NULL` or a column name to store the reporting time
#' `report = accident + delay`.
#'
#' @return Truncated `data`.
#' The reporting time is stored in a colnumn named by `.report_col` unless
#' `.report_col` is `NULL`.
#' If both `.report_col` is `NULL` and `time` contains only `Inf`s,
#' a warning will be issued since `data` will be
#' returned unchanged and no work will be done.
#'
#' @examples
#' claims_full <- data.frame(
#'   acc = runif(100),
#'   repdel = rexp(100)
#' )
#' tau <- 2.0
#' truncate_claims(claims_full, acc, repdel, tau)
#'
#' @export
truncate_claims <- function(data, accident, delay, time,
                            .report_col = "report") {
  assert_that(is_string(.report_col) || is_null(.report_col),
              msg = "`.report_col` must be NULL or a string.")

  accident <- eval_tidy(enquo(accident), data = data)
  delay <- eval_tidy(enquo(delay), data = data)
  time <- eval_tidy(enquo(time), data = data)

  check_lengths(accident, delay, time, .len = nrow(data),
                .msg = "nrow(data)")

  # Sanity check
  if (is.null(.report_col) && all(is.infinite(time))) {
    warning("`time` is Inf and `.report_col` is NULL. ",
            "No work will be done.")
    return(data)
  }

  obs <- accident + delay <= time
  if (!is.null(.report_col)) {
    data[[.report_col]] <- accident + delay
  }

  data[obs, , drop = FALSE]
}

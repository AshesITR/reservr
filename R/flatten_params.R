#' Flatten / Inflate parameter lists / vectors
#'
#' @param params A named list of parameters to be flattened.
#' Should be in a form to be passed as the `with_params` argument to most
#' distribution functions.
#'
#' @return `flatten_params` returns a 'flattened' vector of parameters.
#' It is intended as an adapter for multi-dimensional optimisation functions
#' to distribution objects.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' mm <- dist_mixture(list(
#'   dist_exponential(NULL),
#'   dist_lognormal(0.5, NULL)
#' ), list(NULL, 1))
#'
#' ph <- mm$get_placeholders()
#' ph_flat <- flatten_params(ph)
#' ph_reinflated <- inflate_params(ph_flat)
#' ph_flat[] <- c(1, 1, 6)
#' ph_sample <- inflate_params(ph_flat)
#'
#' x <- mm$sample(
#'   100,
#'   with_params = ph_sample
#' )
#'
#' emp_cdf <- ecdf(x)
#'
#' ggplot(data.frame(t = seq(from = min(x), to = max(x), length.out = 100))) %+%
#'   geom_point(aes(x = t, y = emp_cdf(t))) %+%
#'   geom_line(aes(x = t, y = mm$probability(t, with_params = ph_sample)),
#'             linetype = 2)
flatten_params <- function(params) {
  drop(flatten_params_matrix(params))
}

#' @rdname flatten_params
#'
#' @return `flatten_params_matrix` returns a 'flattened' matrix of parameters.
#' It is intended as an adapter for multi-dimensional optimisation functions
#' to distribution objects. Each column corresponds to one input element.
#'
#' @export
flatten_params_matrix <- function(params) {
  # nolint start
  #
  # Input      | Output
  #     --     |     --
  # named list | c(.nm = ..., .nm2 = ...)
  # list       | c([1] = ..., [2] = ...)
  # list()     | numeric()
  # NULL       | NA_real_
  #
  # nolint end

  flatten_element_matrix <- function(elem) {
    if (is.null(elem)) {
      matrix(nrow = 1L, ncol = 1L)
    } else if (length(elem) == 0) {
      matrix(nrow = 0L, ncol = 1L)
    } else if (is.list(elem) &&
      length(elem) == 2L &&
      hasName(elem, "dist") &&
      hasName(elem, "params") &&
      is.Distribution(elem$dist)) {
      flatten_element_matrix(elem$params)
    } else if (is.list(elem)) {
      if (!is.null(names(elem))) {
        nms <- sprintf(".%s", names(elem))
      } else {
        nms <- sprintf("[%d]", seq_along(elem))
      }
      ress <- lapply(elem, flatten_element_matrix)
      nrow_ress <- vapply(ress, nrow, integer(1L))
      if (any(nrow_ress > 0L)) {
        ress <- ress[nrow_ress > 0L]
        nms <- nms[nrow_ress > 0L]
      }
      nmss <- mapply(
        function(nm, res) {
          paste0(nm, colnames(res))
        },
        nm = nms,
        res = ress,
        SIMPLIFY = FALSE
      )
      res <- do.call(cbind, ress)
      colnames(res) <- do.call(c, nmss)
      res
    } else if (is.numeric(elem)) {
      as.matrix(elem)
    }
  }

  res <- flatten_element_matrix(params)
  colnames(res) <- substring(colnames(res), first = 2)
  res
}

#' @rdname flatten_params
#'
#' @param bounds List of parameter bounds as returned by
#' `dist$get_param_bounds()`
#'
#' @return `flatten_bounds` returns a named list of vectors with names `lower`
#' and `upper`. Containing the upper and lower bounds of each parameter.
#' @export
flatten_bounds <- function(bounds) {
  flatten_interval <- function(elem) {
    if (is.Interval(elem)) {
      elem$range
    } else if (is.list(elem)) {
      lapply(elem, flatten_interval)
    } else {
      NULL
    }
  }
  bounds <- lapply(bounds, flatten_interval)
  bounds <- flatten_params_matrix(bounds)
  list(
    lower = bounds[1, ],
    upper = bounds[2, ]
  )
}

#' @rdname flatten_params
#'
#' @param flat_params A named numeric vector of parameters
#'
#' @return `inflate_params` returns an 'inflated' list of parameters.
#' This can be passed as the `with_params` argument to most distribution
#' functions.
#'
#' @export
inflate_params <- function(flat_params) {

  if (length(flat_params) == 0L) return(list())

  inflate_element <- function(name_parts, values) {
    # Handle raw numerics
    if (any(lengths(name_parts) == 0)) {
      # If there is a raw number in the values, there cannot be a named element.
      stopifnot(all(lengths(name_parts)) == 0)
      return(unname(values))
    }

    # Split the values according to the first element of name_parts
    # Note that empty names (length(name_parts[[i]]) == 0) will be omited in the
    # split and handled later.
    elements <- split(
      seq_along(values),
      vapply(name_parts, `[`, character(1), i = 1)
    )

    res <- lapply(
      elements,
      function(idx) {
        # Inflate each child with its name_parts, removing the head of name
        # parts
        inflate_element(
          name_parts = lapply(name_parts[idx], `[`, -1),
          values = values[idx]
        )
      }
    )

    # An unnamed list is flattened to names [1], ..., [imax]
    # Its name parts will be "1]", ..., "imax]", because strsplit removes the [
    # we used to split on.
    if (length(res) && endsWith(names(res)[1], "]")) {
      idxes <- as.numeric(substr(
        x = names(res),
        start = 1,
        stop = nchar(names(res)) - 1
      ))
      len <- max(idxes)
      old_res <- res
      res <- vector(mode = "list", length = len)
      res[idxes] <- unname(old_res)
    }

    res
  }

  # Create a list with name parts per entry.
  # We will "shave" the head of these qualifiers in turn.
  name_parts <- strsplit(names(flat_params), "[\\[\\.]")

  inflate_element(name_parts, flat_params)
}

#' Subset parameters
#'
#' @param params Params
#' @param idx A logical vector for subsetting.
#' Only numeric vectors within params and its sub-lists will be subsetted and
#' only if those have length equal to `idx`.
#'
#' @return Params with all vectors subsetted by `idx`
#'
#' @noRd
pick_params_at <- function(params, idx) {
  if (is.numeric(params) && length(params) == length(idx)) {
    params[idx]
  } else if (is.list(params)) {
    lapply(params, pick_params_at, idx = idx)
  } else {
    params
  }
}

#' Subset parameters by index
#'
#' @param params Params
#' @param idx A zero-based index vector for subsetting.
#' Only numeric vectors within params and its sub-lists will be subsetted and
#' only if those have length greater than 1.
#'
#' @return Params with all vectors subsetted by `idx`
#'
#' @noRd
pick_params_at_idx <- function(params, idx) {
  if (is.numeric(params) && length(params) > 1L) {
    params[idx + 1]
  } else if (is.list(params)) {
    lapply(params, pick_params_at_idx, idx = idx)
  } else {
    params
  }
}

empty_derivative <- function(ph) {
  empty_and_compact <- function(lst) {
    nulls <- vapply(lst, is.null, logical(1L))
    lsts <- vapply(lst, is.list, logical(1L))
    empty_lsts <- lsts & (lengths(lst) == 0L)
    nonempty_lsts <- lsts & !empty_lsts

    lst[nulls] <- setNames(
      rep(list(numeric()), sum(nulls)),
      names(lst)[nulls]
    )
    lst[nonempty_lsts] <- lapply(lst[nonempty_lsts], empty_and_compact)

    lst[!empty_lsts]
  }

  empty_and_compact(ph)
}

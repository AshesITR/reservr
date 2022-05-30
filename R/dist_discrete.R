#' Discrete Distribution
#'
#' A full-flexibility discrete distribution with values from 1 to `size`.
#'
#' Both parameters can be overridden with
#' `with_params = list(size = ..., probs = ...)`.
#'
#' @param size Number of classes parameter (integer), or `NULL` as a
#'   placeholder.
#' @param probs Vector of probabilties parameter, or `NULL` as a placeholder.
#'
#' @return A `DiscreteDistribution` object.
#' @export
#'
#' @examples
#' d_discrete <- dist_discrete(size = 4, probs = list(0.5, 0.25, 0.15, 0.1))
#' x <- d_discrete$sample(100)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_discrete,
#'   estimated = d_discrete,
#'   with_params = list(
#'     estimated = list(
#'       size = max(x),
#'       probs = as.list(unname(table(x)) / 100)
#'     )
#'   ),
#'   .x = 0:max(x)
#' )
#'
#' @family Distributions
dist_discrete <- function(size = NULL, probs = NULL) {
  if (is.null(probs)) probs <- vector("list", as.integer(size))
  DiscreteDistribution$new(size = size, probs = probs)
}

DiscreteDistribution <- distribution_class(
  name = "Discrete",
  type = "discrete",
  params = list(
    size = I_POSITIVE_INTEGERS,
    probs = list(I_UNIT_INTERVAL)
  ),
  sample = function(n, params) {
    k <- params$size[[1L]]
    slot <- runif(n)
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]
    rowSums(slot > probmat) + 1L
  },
  density = function(x, log = FALSE, params) {
    params$probs <- lapply(params$probs, rep_len, length(x))
    densmat <- do.call(cbind, params$probs)
    densmat <- densmat / rowSums(densmat)

    x_ok <- as.integer(x)
    x_ok[x < 1L] <- NA_integer_
    x_ok[x > params$size[[1L]]] <- NA_integer_
    x_ok[!is_integerish(x)] <- NA_integer_

    res <- densmat[seq_along(x) + (x_ok - 1L) * nrow(densmat)]
    res[is.na(x_ok)] <- 0.0
    if (log) log(res) else res
  },
  probability = function(q, lower.tail = TRUE, log.p = FALSE, params) {
    k <- params$size[[1L]]
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]

    q <- floor(q)
    q[q > k] <- k
    q[q < 0] <- 0

    res <- numeric(length(q))
    res[q > 0] <- probmat[
      which(q > 0.0) + (q[q > 0.0] - 1.0) * nrow(probmat)
    ]
    res[q == 0] <- 0.0

    if (!lower.tail) res <- 1.0 - res
    if (log.p) res <- log(res)

    res
  },
  quantile = function(p, lower.tail = TRUE, log.p = FALSE, params) {
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1.0 - p

    k <- params$size[[1L]]
    probmat <- matrixStats::rowCumsums(do.call(cbind, params$probs))
    probmat <- probmat / probmat[, k]
    rowSums(p > probmat) + 1L
  },
  get_dof = function() {
    sdof <- super$get_dof()
    if (length(self$get_placeholders()$probs)) {
      sdof <- sdof - 1L
    }
    sdof
  },
  support = function(x, params) {
    is_integerish(x) & x > 0.0 & x <= params$size[[1L]]
  },
  tf_logdensity = function() {
    k <- as.integer(self$default_params$size)
    if (is.null(k)) stop("`size` must be fixed for tf_logdensity().")
    function(x, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))

      denss <- tf$stack(lapply(seq_len(k), function(i) {
        tf$where(x == k_constant(i), probs[, i], K$zero)
      }), axis = 1L)

      dens <- tf$reduce_sum(denss, axis = 1L)
      dens_safe <- tf$where(dens > 0.0, dens, 1.0)
      tf$where(dens > 0.0, log(dens_safe), K$neg_inf)
    }
  },
  tf_logprobability = function() {
    k <- as.integer(self$default_params$size)
    if (is.null(k)) stop("`size` must be fixed for tf_logprobability().")
    function(qmin, qmax, args) {
      probs <- args[["probs"]]
      probs <- tf$reshape(probs, list(-1L, k))

      all_probs <- tf$stack(lapply(seq_len(k), function(i) {
        tf$where(
          qmin <= i & qmax >= i,
          probs[, i],
          K$zero
        )
      }), axis = 1L)

      prob <- tf$reduce_sum(all_probs, axis = 1L)
      prob_safe <- tf$where(prob > 0.0, prob, 1.0)
      tf$where(prob > 0.0, log(prob_safe), K$neg_inf)
    }
  },
  tf_make_constants = function(with_params = list()) {
    check_installed("keras")
    params <- private$.make_params(with_params, 1)
    out <- list()
    if (!is.null(params$size)) {
      out$size <- keras::k_constant(params$size)
    }
    if (length(params$probs) && !is.null(params$probs[[1L]])) {
      probs <- as.numeric(params$probs)
      out$probs <- keras::k_constant(probs / sum(probs), shape = c(1L, length(probs)))
    }
    out
  },
  tf_compile_params = function(input, name_prefix = "") {
    ph <- self$get_placeholders()
    k <- length(ph$probs)
    if (length(ph$probs)) {
      out <- list(
        probs = keras::layer_dense(
          input, units = k, activation = "softmax",
          name = paste0(name_prefix, "probs")
        )
      )
    } else {
      out <- list()
    }

    list(
      outputs = out,
      output_inflater = eval(bquote(function(outputs) {
        if (!is.list(outputs)) outputs <- list(outputs)
        list(probs = outputs[[1L]])
      }))
    )
  }
)

#' @export
fit_dist_start.DiscreteDistribution <- function(dist, obs, ...) {
  obs <- as_trunc_obs(obs)
  res <- dist$get_placeholders()
  ph_probs <- length(res$probs) > 0L
  ph <- names(res)

  if ("size" %in% ph) {
    size <- max(obs$xmax)
    res$size <- size
  } else {
    size <- dist$get_params()$size
  }

  if (ph_probs) {
    densmat <- map_dbl_matrix(
      seq_len(size),
      function(i) {
        as.numeric(obs$xmin <= i & obs$xmax >= i)
      },
      nrow(obs)
    )
    res$probs <- local({
      p0 <- colSums(obs$w * (densmat / rowSums(densmat)))
      as.list(p0 / sum(p0))
    })
  }

  res
}

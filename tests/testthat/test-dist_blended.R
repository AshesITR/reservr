test_that("test dist_blended", {
  set.seed(1337L)

  dist <- dist_blended(
    list(
      dist_exponential(),
      dist_genpareto()
    )
  )

  params <- list(
    probs = list(0.9, 0.1),
    dists = list(
      list(rate = 2.0),
      list(u = 1.5, xi = 0.2, sigmau = 1.0)
    ),
    breaks = list(1.5),
    bandwidths = list(0.3)
  )

  x <- dist$sample(100L, with_params = params)

  dist$default_params$breaks <- params$breaks
  dist$default_params$bandwidths <- params$bandwidths
  expect_silent(fit(dist, x))
  expect_identical(dist$get_type(), "continuous")
  expect_length(dist$get_components(), 2L)

  p_lower <- pexp(params$breaks[[1L]], rate = params$dists[[1L]]$rate)

  x_lhs <- x[x < params$breaks[[1L]] - params$bandwidths[[1L]]]
  x_rhs <- x[x > params$breaks[[1L]] + params$bandwidths[[1L]]]

  expect_density(
    dist,
    function(x, log = FALSE, ...) {
      if (log) {
        log(list(...)$probs[[1L]]) +
          dexp(x, rate = list(...)$dists[[1L]]$rate, log = TRUE) -
          pexp(list(...)$breaks[[1L]], rate = list(...)$dists[[1L]]$rate, log = TRUE)
      } else {
        list(...)$probs[[1L]] *
          dexp(x, rate = list(...)$dists[[1L]]$rate) /
          pexp(list(...)$breaks[[1L]], rate = list(...)$dists[[1L]]$rate)
      }
    },
    params,
    x_lhs
  )

  expect_density(
    dist,
    function(x, log = FALSE, ...) {
      params_gpd <- list(...)$dists[[2L]]
      if (log) {
        log(list(...)$probs[[2L]]) +
          do.call(dgpd, c(list(x = x, log = TRUE), params_gpd)) -
          do.call(pgpd, c(list(q = list(...)$breaks[[1L]], lower.tail = FALSE, log = TRUE), params_gpd))
      } else {
        list(...)$probs[[2L]] *
          do.call(dgpd, c(list(x = x), params_gpd)) /
          do.call(pgpd, c(list(q = list(...)$breaks[[1L]], lower.tail = FALSE), params_gpd))
      }
    },
    params,
    x_rhs
  )

  expect_probability(
    dist,
    function(q, log.p = FALSE, lower.tail = TRUE, ...) {
      pr <- list(...)$probs[[1L]] *
        pexp(q, rate = list(...)$dists[[1L]]$rate) /
        pexp(list(...)$breaks[[1L]], rate = list(...)$dists[[1L]]$rate)
      if (!lower.tail) pr <- 1 - pr
      if (log.p) pr <- log(pr)
      pr
    },
    params,
    x_lhs
  )

  expect_probability(
    dist,
    function(q, log.p = FALSE, lower.tail = TRUE, ...) {
      params_gpd <- list(...)$dists[[2L]]
      pr <- (
        list(...)$probs[[1L]] +
          list(...)$probs[[2L]] *
          do.call(pgpd, c(list(q = q), params_gpd)) /
            do.call(pgpd, c(list(q = params$breaks[[1L]], lower.tail = FALSE), params_gpd))
      )
      if (!lower.tail) pr <- 1 - pr
      if (log.p) pr <- log(pr)
      pr
    },
    params,
    x_rhs
  )

  expect_identical(
    dist$is_in_support(x, with_params = params),
    rep_len(TRUE, length(x))
  )

  expect_tf_logdensity(dist, params, x)
  expect_tf_logprobability(dist, params, x, x + 1.0)
  expect_tf_logprobability(dist, params, 0, x)
  expect_tf_logprobability(dist, params, x, Inf)
})
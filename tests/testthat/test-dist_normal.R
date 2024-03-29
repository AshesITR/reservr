test_that("normal distribution works", {
  set.seed(1337L)
  dist <- dist_normal()
  params <- list(mean = 1, sd = 2)
  x <- dist$sample(100L, with_params = params)

  expect_silent(fit(dist, x))
  expect_identical(dist$get_type(), "continuous")
  expect_density(dist, dnorm, params, x)
  expect_probability(dist, pnorm, params, x)
  expect_quantile(dist, qnorm, params)
  expect_identical(dist$is_in_support(x), rep_len(TRUE, length(x)))
  expect_diff_density(dist, x, params)
  expect_diff_density(dist, x, list(mean = 3, sd = 5))
  expect_diff_probability(dist, x, params)
  expect_diff_probability(dist, x, list(mean = 3, sd = 5))
  expect_tf_logdensity(dist, params, x)
  expect_tf_logprobability(dist, params, x, x + 1.0)
  expect_tf_logprobability(dist, params, 0, x)
  expect_tf_logprobability(dist, params, x, Inf)

  expect_tf_fit(dist, params, I_REALS)
})

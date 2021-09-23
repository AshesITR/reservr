test_that("erlang mixture distribution works", {
  set.seed(1337L)
  dist <- dist_erlangmix(list(NULL, NULL, NULL))
  params <- list(
    shapes = list(1L, 4L, 12L),
    scale = 2.0,
    probs = list(0.5, 0.3, 0.2)
  )
  alt_params <- list(
    shapes = list(2L, 6L, 100L),
    scale = 0.1,
    probs = list(0.7, 0.2, 0.1)
  )
  x <- dist$sample(100L, with_params = params)

  expect_silent(fit(dist, x, init = "shapes",
                    shapes = as.numeric(params$shapes)))
  expect_silent(fit(dist, x, init = "fan",
                    spread = 3L))
  expect_silent(fit(dist, x, init = "kmeans"))
  expect_silent(fit(dist, x, init = "cmm"))

  expect_identical(dist$get_type(), "continuous")
  # TODO test density and probability for correctness
  # use gamma as reference
  expect_identical(dist$is_in_support(x), rep_len(TRUE, length(x)))
  expect_diff_density(dist, x, params)
  expect_diff_density(dist, x, alt_params)
  # TODO implement
  #> expect_diff_probability(emix, x, params)
  #> expect_diff_probability(emix, x, alt_params)
  expect_tf_logdensity(dist, params, x)
  # Extreme shapes cause greater numeric instability.
  expect_tf_logdensity(dist, alt_params, x, tolerance = 1.0e-5)
  expect_tf_logprobability(dist, params, x, x + 1.0)
  expect_tf_logprobability(dist, params, x, rep_len(Inf, 100L))
  expect_tf_logprobability(dist, params, rep_len(0, 100L), x)
  # Extreme outliers can't be handled, so we need a good sample
  x_alt <- dist$sample(100L, with_params = alt_params)
  expect_tf_logprobability(dist, alt_params, x_alt, x_alt + 1.0)
})

test_that("can use erlang mixtures with 1 component", {
  set.seed(1337L)
  dist <- dist_erlangmix(list(NULL))
  params <- list(
    shapes = list(3L),
    scale = 3.0,
    probs = list(1.0)
  )
  dist_equiv <- dist_gamma()
  params_equiv <- list(
    shape = params$shape[[1L]],
    rate = 1.0 / params$scale
  )

  x <- dist$sample(100L, with_params = params)

  expect_equal(
    dist$density(x, with_params = params),
    dist_equiv$density(x, with_params = params_equiv)
  )
  expect_equal(
    dist$probability(x, with_params = params),
    dist_equiv$probability(x, with_params = params_equiv)
  )
})

test_that("numerically unstable tf fitting works", {
  skip_if_no_tensorflow()
  set.seed(2350L)

  dist <- dist_erlangmix(list(1, 50))
  params <- list(probs = list(0.9, 0.1), scale = 20)
  N <- 1000L
  x <- dist$sample(N, with_params = params)
  tensorflow::tf$keras$backend$set_floatx("float32")
  on.exit({ tensorflow::tf$keras$backend$set_floatx("float64") })
  tmax <- runif(N, min = 20, max = 50)

  obs <- truncate_obs(x, tmin_new = 0, tmax_new = tmax)

  rand_input <- k_matrix(runif(nrow(obs)))

  tf_in <- keras::layer_input(1L)
  mod <- tf_compile_model(
    inputs = list(tf_in),
    intermediate_output = tf_in,
    dist = dist,
    optimizer = keras::optimizer_adam()
  )

  tf_initialise_model(mod, params, mode = "zero")

  expect_silent({
    tf_fit <- fit(
      object = mod,
      x = rand_input,
      y = obs,
      epochs = 10L,
      callbacks = list(
        callback_debug_dist_gradients(mod, rand_input, obs)
      )
    )
  })

})

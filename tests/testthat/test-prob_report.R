test_that("test prob_report", {
  expect_equal(
    prob_report(
      dist_dirac(0),
      data.frame(
        xmin = 0.0,
        xmax = 1.0,
        tmin = c(0.0, 1.0),
        tmax = c(1.0, 2.0)
      )
    ),
    c(1, 0)
  )

  expect_equal(
    prob_report(
      dist_dirac(0),
      data.frame(
        xmin = 0.0,
        xmax = 1.0,
        tmin = c(0.0, 1.0),
        tmax = c(1.0, 2.0)
      ),
      expo = identity
    ),
    c(1, 0)
  )

  expect_equal(
    prob_report(
      dist_uniform(0, 1),
      data.frame(
        xmin = 0,
        xmax = 1.0,
        tmin = c(0.0, 1.0, 0.0),
        tmax = c(1.0, 2.0, 2.0)
      )
    ),
    c(0.5, 0.5, 1.0)
  )
})

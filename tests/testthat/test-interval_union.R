test_that("interval_union computes correct unions", {
  expect_equal(
    as.list(interval_union(interval(0, 1), interval(1, 3))),
    list(interval(0, 1), interval(1, 3))
  )

  expect_equal(
    as.list(interval_union(interval(0, 1, closed = TRUE), interval(1, 3))),
    list(interval(0, 3, include_lowest = TRUE))
  )

  expect_equal(
    as.list(interval_union(interval(0, 1, closed = TRUE), interval(2, 3, closed = TRUE))),
    list(
      interval(0, 1, closed = TRUE), interval(2, 3, closed = TRUE)
    )
  )

  expect_equal(
    as.list(interval_union(
      interval(0, 1, closed = TRUE, integer = TRUE),
      interval(2, 3, closed = TRUE, integer = TRUE)
    )),
    list(interval(0, 3, closed = TRUE, integer = TRUE))
  )

  expect_equal(
    as.list(interval_union(interval(0, 1, closed = TRUE), interval(0.5, 2, include_lowest = TRUE))),
    list(interval(0, 2, include_lowest = TRUE))
  )

  expect_equal(
    as.list(interval_union(
      interval(0, 1, closed = TRUE, integer = TRUE),
      interval(0.5, 2, include_lowest = TRUE, integer = TRUE),
      2.5, 3.5, -1.0
    )),
    list(interval(-1, 1, closed = TRUE, integer = TRUE), interval(0.5, 3.5, closed = TRUE, integer = TRUE))
  )

  expect_equal(
    as.list(interval_union(0, 1, 2, 3:6)),
    list(interval(c(0, 6), closed = TRUE, integer = TRUE))
  )

  expect_equal(nrow(interval_union()), 0L)
})

test_that("interval_union prints correctly", {
  expect_output(print(interval_union()), "{}", fixed = TRUE)
  expect_output(print(interval_union(1)), "{{1}}", fixed = TRUE)
  expect_output(print(interval_union(1, 2)), "{[1, 2] (int)}", fixed = TRUE)
  expect_output(print(interval_union(interval(1, 2))), "{(1, 2)}", fixed = TRUE)
})

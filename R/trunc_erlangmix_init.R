# TODO embed into fit_dist_start.ErlangMixtureDistribution
.trunc_erlangmix_init <- function(x, init, num_components, spread = 1L, shapes) {
  init <- match.arg(init, c("shapes", "fan", "kmeans", "cmm"))

  switch(init,
    shapes = {
      assert_that(
        is_integerish(shapes),
        all(shapes >= 1L),
        msg = paste0(
          "`init` = \"shapes\" requires manual specification of ",
          "`shapes` >= 1L."
        )
      )

      shapes <- sort(unique(shapes))
      num_components <- length(shapes)
    },
    fan = {
      assert_that(
        is_scalar_integerish(num_components),
        num_components >= 1L,
        msg = paste0(
          "`init` = \"fan\" requires manual specification of ",
          "`num_components` >= 1L."
        )
      )

      shapes <- 1L + (seq_len(num_components) - 1L) * spread
    },
    kmeans = {
      assert_that(
        is_scalar_integerish(num_components),
        num_components >= 1L,
        msg = paste0(
          "`init` = \"kmeans\" requires manual specification of ",
          "`num_components` >= 1L."
        )
      )

      cts <- kmeans(x, num_components)$centers
      shapes <- as.integer(sort(cts) / min(cts, diff(sort(cts))))
    },
    cmm = {
      assert_that(
        is_scalar_integerish(num_components),
        num_components >= 1L,
        msg = paste0(
          "`init` = \"cmm\" requires manual specification of ",
          "`num_components` >= 1L."
        )
      )

      clust <- kmeans(x, num_components)
      probs <- clust$size / length(x)
      means <- as.vector(clust$centers)
      ord <- order(means)

      e_x <- mean(x)
      e_x_2 <- mean(x^2)

      scale <- min((e_x_2 - sum(probs * means^2)) / e_x, means)

      shapes <- ceiling(means / scale)

      return(list(
        shapes = shapes[ord],
        scale = scale,
        probs = probs[ord]
      ))
    }
  )

  n_obs <- length(x)

  scale <- max(x) / shapes[num_components]
  bin <- .bincode(x, c(0, scale * shapes))
  # fix cases where numerically max(x) / shapes[m] * shapes[m] < max(x)
  bin[is.na(bin)] <- num_components
  probs <- tabulate(bin) / n_obs
  # Better starting value for scale via method of moments
  scale <- mean(x) / weighted.mean(shapes, probs)

  list(
    shapes = shapes,
    scale = scale,
    probs = probs
  )
}

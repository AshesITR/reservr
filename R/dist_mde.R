#' Construct a MDE-Family
#'
#' Constructs a MDE-Family distribution with fixed number of components.
#'
#' @param n Number of dirac components, starting with a point mass at 0.
#' @param m Number of erlang components, translated by `n - 0.5`.
#'
#' @export
#'
#' @return A `MixtureDistribution` containing `n` `DiracDistribution`s at 0 .. n - 1 and a `TranslatedDistribution`
#' with offset `n - 0.5` of an `ErlangMixtureDistribution` with `m` shapes.
#'
#' @examples
#' params <- list(
#'   dists = list(
#'     list(),
#'     list(dist = list(scale = 1.0, probs = list(0.5, 0.3, 0.2), shapes = list(1L, 2L, 3L)))
#'   ),
#'   probs = list(0.1, 0.9)
#' )
#' dist <- dist_mde(1, 3)
#' x <- dist$sample(100, with_params = params)
#' d_emp <- dist_empirical(x)
#'
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = dist,
#'   with_params = list(
#'     theoretical = params
#'   ),
#'   .x = seq(1e-4, 5, length.out = 100)
#' )
#'
#' @family Distributions
dist_mde <- function(n, m) {
  assert_that(
    is_scalar_integerish(n, finite = TRUE),
    n >= 0,
    msg = "`n` must be a non-negative integer."
  )
  assert_that(
    is_scalar_integerish(m, finite = TRUE),
    m >= 0,
    msg = "`m` must be a non-negative integer."
  )

  diracs <- lapply(seq_len(n) - 1.0, dist_dirac)
  erlangs <- list(dist_translate(
    dist_erlangmix(shapes = vector("list", m)),
    offset = n - 0.5
  ))

  dist_mixture(dists = c(diracs, erlangs))
}

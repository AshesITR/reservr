#' Construct a BDEGP-Family
#'
#' Constructs a BDEGP-Family distribution with fixed number of components and
#' blending interval.
#'
#' @param n Number of dirac components, starting with a point mass at 0.
#' @param m Number of erlang components, translated by `n - 0.5`.
#' @param u Blending cut-off, must be a positive real.
#' @param epsilon Blending radius, must be a positive real less than `u`.
#' The blending interval will be `u - epsilon < x < u + epsilon`.
#'
#' @return
#'  - A `MixtureDistribution` of
#'    + `n` `DiracDistribution`s at 0 .. n - 1 and
#'    + a `BlendedDistribution` object with child Distributions
#'      - a `TranslatedDistribution` with offset `n - 0.5` of an `ErlangMixtureDistribution` with `m` shapes
#'      - and a `GeneralizedParetoDistribution` with shape parameter restricted to \[0, 1] and location parameter fixed
#'        at `u`
#'      With break `u` and bandwidth `epsilon`.
#'
#' @export
#'
#' @examples
#' dist <- dist_bdegp(n = 1, m = 2, u = 10, epsilon = 3)
#' params <- list(
#'   dists = list(
#'     list(
#'       dist = list(probs = list(1.0))
#'     ),
#'     list(
#'       dists = list(
#'         list(
#'           dist = list(
#'             shapes = list(1L, 2L),
#'             scale = 1.0,
#'             probs = list(0.7, 0.3)
#'           )
#'         ),
#'         list(
#'           sigmau = 1.0,
#'           xi = 0.1
#'         )
#'       ),
#'       probs = list(0.1, 0.9)
#'     )
#'   ),
#'   probs = list(0.95, 0.05)
#' )
#' x <- dist$sample(100, with_params = params)
#'
#' plot_distributions(
#'   theoretical = dist,
#'   empirical = dist_empirical(x),
#'   .x = seq(0, 20, length.out = 101),
#'   with_params = list(theoretical = params)
#' )
#'
#' @family Distributions
dist_bdegp <- function(n, m, u, epsilon) {
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
  assert_that(
    is_scalar_double(u),
    is.finite(u),
    u > 0,
    msg = "`u` must be a positive real."
  )
  assert_that(
    is_scalar_double(epsilon),
    epsilon > 0,
    epsilon < u,
    msg = "`epsilon` must be a positive real < `u`."
  )

  erlangs <- dist_translate(
    dist_erlangmix(shapes = vector("list", m)),
    offset = n - 0.5
  )

  dist <- dist_mixture(
    dists = list(
      dist_translate(
        dist_discrete(size = n),
        offset = -1.0
      ),
      dist_blended(
        dists = list(erlangs, dist_genpareto1(u = u)),
        breaks = list(u),
        bandwidths = list(epsilon)
      )
    )
  )
  class(dist) <- c("BDEGPDistribution", class(dist))
  dist
}

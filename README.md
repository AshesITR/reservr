
# reservr

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/AshesITR/reservr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/AshesITR/reservr?branch=master)
[![R-CMD-check](https://github.com/AshesITR/reservr/workflows/R-CMD-check/badge.svg)](https://github.com/AshesITR/reservr/actions)
<!-- badges: end -->

The goal of reservr is to provide a flexible interface for specifying distributions and fitting them to (randomly) 
truncated and possibly interval-censored data.
It provides custom fitting algorithms to fit distributions to i.i.d. samples as well as dynnamic TensorFlow
integration to allow training neural networks with arbitrary output distributions.
The latter can be used to include explanatory variables in the distributional fits.
Reservr also provides some tools relevant for working with its core functionality in an actuarial setting, namely the
functions `prob_report()` and `truncate_claims()`, both of which make assumptions on the type of random truncation
applied to the data.

Please refer to the vignettes `distributions.Rmd` and `tensorflow.Rmd` for detailed introductions.

## Installation

reservr is not yet on CRAN.
You can install the latest development version of reservr via
``` r
devtools::install_github("AshesITR/reservr")
```

You can install the released version of reservr from [CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("reservr")
```

If you want to use all of reservrs features, make sure to also install
[tensorflow](https://tensorflow.rstudio.com/installation/).

## Example

This is a basic example which shows how to fit a normal distribution to randomly truncated and censored data.

``` r
library(reservr)
set.seed(123)
mu <- 0
sigma <- 1
N <- 1000
p_cens <- 0.8

x <- rnorm(N, mean = mu, sd = sigma)
is_censored <- rbinom(N, size = 1L, prob = p_cens) == 1L
x_lower <- x
x_lower[is_censored] <- x[is_censored] - runif(sum(is_censored), min = 0, max = 0.5)
x_upper <- x
x_upper[is_censored] <- x[is_censored] + runif(sum(is_censored), min = 0, max = 0.5)

t_lower <- runif(N, min = -2, max = 0)
t_upper <- runif(N, min = 0, max = 2)

is_observed <- t_lower <= x & x <= t_upper

obs <- trunc_obs(
  xmin = pmax(x_lower, t_lower)[is_observed],
  xmax = pmin(x_upper, t_upper)[is_observed],
  tmin = t_lower[is_observed],
  tmax = t_upper[is_observed]
)

# Summary of the simulation
cat(sprintf(
  "simulated samples: %d\nobserved samples: %d\ncensored samples: %d\n", 
  N, nrow(obs), sum(is.na(obs$x))
))

# Define outcome distribution and perform fit to truncated and (partially) censored sample
dist <- dist_normal()
the_fit <- fit(dist, obs)

# Visualize resulting parameters and show a kernel density estimate of the samples.
# We replace interval-censored samples with their midpoint for the kernel density estimate.
plot_distributions(
  true = dist,
  fitted = dist, 
  empirical = dist_empirical(0.5 * (obs$xmin + obs$xmax)), 
  .x = seq(-5, 5, length.out = 201), 
  plots = "density", 
  with_params = list(
    true = list(mean = mu, sd = sigma), 
    fitted = the_fit$params
  )
)
```

## Code of Conduct

Please note that the reservr project is released with a
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

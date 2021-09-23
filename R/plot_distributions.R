#' Plot several distributions
#'
#' @param ... distribution objects (must be named)
#' @param distributions Named list of distribution objects.
#' This is concatenated with `...`.
#' @param .x Numeric vector of points to evaluate at.
#' @param plots Plots to be created. May be abbreviated.
#' The plots will be stacked in the order given from top to bottom.
#' @param with_params list of distribution parameters to be given to each
#' distribution using `with_params`. If named, the names are matched to the
#' distribution names. Otherwise, they are allocated positionally, index 1
#' corresponding to the first element of `distributions`, then all other
#' elements from `distributions` followed by the arguments in `...` in order.
#' @param as_list return a list of ggplots instead of a patchwork?
#'
#' @return A stacked patchwork of the requested ggplots
#' @export
#'
#' @examples
#' rate <- 1
#' x <- rexp(20, rate)
#' d_emp <- dist_empirical(x, positive = TRUE)
#' d_exp <- dist_exponential()
#' plot_distributions(
#'   empirical = d_emp,
#'   theoretical = d_exp,
#'   estimated = d_exp,
#'   with_params = list(
#'     theoretical = list(rate = rate),
#'     estimated = list(rate = 1 / mean(x))
#'   ),
#'   .x = seq(1e-4, 5, length.out = 100)
#' )
plot_distributions <- function(..., distributions = list(), .x,
                               plots = c("density", "probability", "hazard"),
                               with_params = list(), as_list = FALSE) {
  check_installed(c("ggplot2", "patchwork", "colorspace"))

  if (rlang::dots_n(...) > 0) {
    dot_nms <- names(rlang::quos_auto_name(rlang::quos(...)))
    dot_dists <- list(...)
    names(dot_dists) <- dot_nms
    distributions <- c(distributions, dot_dists)
  }
  if (!rlang::is_named(distributions)) {
    stop("`distributions` must be named.")
  }

  if (rlang::is_named(with_params)) {
    wp_new <- vector(mode = "list", length = length(distributions))
    names(wp_new) <- names(distributions)
    wp_new[names(with_params)] <- with_params
    with_params <- wp_new
  } else if (length(with_params)) {
    check_lengths(
      with_params,
      .len = length(distributions),
      .msg = "number of distributions"
    )
    with_params <- rep_len(with_params, length.out = length(distributions))
    names(with_params) <- names(distributions)
  } else {
    with_params <- lapply(distributions, function(.) list())
  }

  plot_opts <- c("density", "probability", "hazard")
  plots <- plot_opts[pmatch(plots, plot_opts, duplicates.ok = FALSE)]

  plot_data <- expand.grid(
    x = .x,
    distribution = names(distributions),
    stringsAsFactors = FALSE
  )

  res <- list()

  plot_data$discrete <- vapply(
    names(distributions),
    function(dist_name) {
      distributions[[dist_name]]$is_discrete_at(
        x = .x,
        with_params = with_params[[dist_name]]
      )
    },
    logical(length(.x))
  )[seq_len(nrow(plot_data))]

  if ("density" %in% plots || ("probability" %in% plots && any(plot_data$discrete))) {
    plot_data$density <- vapply(
      names(distributions),
      function(dist_name) {
        distributions[[dist_name]]$density(
          x = .x,
          with_params = with_params[[dist_name]]
        )
      },
      numeric(length(.x))
    )[seq_len(nrow(plot_data))]
  }

  if ("probability" %in% plots) {
    plot_data$probability <- vapply(
      names(distributions),
      function(dist_name) {
        distributions[[dist_name]]$probability(
          q = .x,
          with_params = with_params[[dist_name]]
        )
      },
      numeric(length(.x))
    )[seq_len(nrow(plot_data))]

    plot_data$probability_lower <- plot_data$probability

    if (any(plot_data$discrete)) {
      plot_data$probability_lower[plot_data$discrete] <- with(plot_data, probability[discrete] - density[discrete])
    }
  }

  if ("hazard" %in% plots) {
    plot_data$hazard <- vapply(
      names(distributions),
      function(dist_name) {
        distributions[[dist_name]]$hazard(
          x = .x,
          with_params = with_params[[dist_name]]
        )
      },
      numeric(length(.x))
    )[seq_len(nrow(plot_data))]
  }

  for (plot in plots) {
    plot_sym <- rlang::sym(plot)

    # appease R CMD check:
    x <- NULL
    distribution <- NULL
    discrete <- NULL
    probability <- NULL
    probability_lower <- NULL
    group <- NULL

    p <- switch(
      plot,
      density =,
      hazard = {
        plot_data_cont <- plot_data
        plot_data_cont[[plot_sym]][plot_data$discrete] <- 0.0
        plot_data_disc <- plot_data[plot_data$discrete, ]

        ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = {{ plot_sym }}, color = distribution)) +
          ggplot2::geom_line(data = plot_data_cont) +
          ggplot2::geom_point(data = plot_data_disc, show.legend = FALSE) +
          ggplot2::geom_linerange(
            ggplot2::aes(ymax = {{ plot_sym }}, ymin = rep_len(0, nrow(plot_data_disc))), # fixes continuous distributions
            data = plot_data_disc, linetype = 2L, show.legend = FALSE
          ) +
          ggplot2::theme_bw() +
          ggplot2::guides(color = ggplot2::guide_legend(direction = "horizontal")) +
          colorspace::scale_color_discrete_qualitative()
      },
      probability = {
        `%>%` <- dplyr::`%>%`

        plot_data_cont <- plot_data %>%
          dplyr::group_by(distribution) %>%
          dplyr::arrange(x) %>%
          dplyr::mutate(group = cumsum(discrete)) %>%
          dplyr::bind_rows(
            plot_data %>%
              dplyr::filter(discrete) %>%
              dplyr::group_by(distribution) %>%
              dplyr::mutate(group = dplyr::row_number() - 1, probability = probability_lower)
          )

        plot_data_jumps <- plot_data %>%
          dplyr::filter(discrete) %>%
          dplyr::mutate(group = dplyr::row_number())
        plot_data_jumps <- dplyr::bind_rows(
          plot_data_jumps,
          dplyr::mutate(plot_data_jumps, probability = probability_lower)
        )

        ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = probability, color = distribution,
                                               group = factor(paste(distribution, group)))) +
          ggplot2::geom_line(data = plot_data_cont, size = 1) +
          ggplot2::geom_line(data = plot_data_jumps, linetype = 2L) +
          ggplot2::theme_bw() +
          ggplot2::guides(color = ggplot2::guide_legend(direction = "horizontal")) +
          colorspace::scale_color_discrete_qualitative()
      }
    )

    res <- c(res, list(p))
  }

  if (!as_list) {
    res <- Reduce(`+`, res) +
      patchwork::guide_area() +
      patchwork::plot_layout(
        ncol = 1,
        guides = "collect",
        heights = c(rep(5, length(plots)), 1)
      )
  }

  res
}

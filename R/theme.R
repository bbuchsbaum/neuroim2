#' A minimal, publication-friendly theme for image slices
#'
#' Quiet axes, thin panel border, no grid, generous margins, slim legend.
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @importFrom ggplot2 %+replace% rel
#' @export
theme_neuro <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "grey70", fill = NA, linewidth = .4),
      plot.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 6)),
      plot.caption = ggplot2::element_text(hjust = 0, size = rel(.85)),
      legend.key.height = grid::unit(3, "cm"),
      legend.key.width  = grid::unit(.4, "cm"),
      legend.title = ggplot2::element_text(size = rel(.9)),
      legend.position = "right",
      plot.margin = grid::unit(c(6, 6, 6, 6), "pt")
    )
}

#' A minimal, publication-friendly theme for image slices
#'
#' Quiet axes, thin panel border, no grid, generous margins, slim legend.
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param style Visual style, either \code{"light"} or \code{"dark"}.
#' @importFrom ggplot2 %+replace% rel
#' @export
theme_neuro <- function(base_size = 10, base_family = "", style = c("light", "dark")) {
  style <- match.arg(style)
  if (style == "dark") {
    fg <- "grey92"
    muted <- "grey72"
    border <- "grey35"
    bg <- "grey8"
    panel_bg <- "black"
  } else {
    fg <- "grey12"
    muted <- "grey35"
    border <- "grey70"
    bg <- "white"
    panel_bg <- "white"
  }

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = panel_bg, colour = NA),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = border, fill = NA, linewidth = .35),
      plot.title = ggplot2::element_text(
        colour = fg, face = "bold", margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(colour = muted, margin = ggplot2::margin(b = 6)),
      plot.caption = ggplot2::element_text(colour = muted, hjust = 0, size = rel(.85)),
      strip.text = ggplot2::element_text(colour = fg, face = "bold", margin = ggplot2::margin(b = 4)),
      legend.key.height = grid::unit(3, "cm"),
      legend.key.width  = grid::unit(.4, "cm"),
      legend.title = ggplot2::element_text(colour = fg, size = rel(.9)),
      legend.text = ggplot2::element_text(colour = muted),
      legend.background = ggplot2::element_rect(fill = bg, colour = NA),
      legend.position = "right",
      plot.margin = grid::unit(c(6, 6, 6, 6), "pt")
    )
}

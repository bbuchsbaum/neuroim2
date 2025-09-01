#' Neuroimaging color palettes and helpers
#'
#' Lightweight, perceptually-uniform color tools with safe fallbacks.
#' @param name Palette name (e.g., "grays", "viridis", "inferno", "magma",
#'   "plasma", "turbo", "cividis"). Case-insensitive. If you pass a vector of
#'   colors, it's returned unchanged.
#' @param n Number of colors to generate.
#' @return A character vector of hex colors.
#' @export
resolve_cmap <- function(name = "grays", n = 256) {
  if (length(name) > 1L) return(name)              # already a color vector
  if (is.null(name)) name <- "grays"
  nm <- tolower(name)

  # Try grDevices::hcl.colors first (base R, no extra deps)
  hc_try <- function(pal) {
    pal2 <- try(grDevices::hcl.colors(n, palette = pal), silent = TRUE)
    if (inherits(pal2, "try-error")) NULL else pal2
  }

  # Common aliases -> hcl.colors palette names (case-insensitive in practice)
  alias <- c(
    grays   = "Grays", grey = "Grays", gray = "Grays",
    viridis = "Viridis",
    inferno = "Inferno",
    magma   = "Magma",
    plasma  = "Plasma",
    turbo   = "Turbo",
    cividis = "Cividis"
  )

  if (nm %in% names(alias)) {
    pal <- hc_try(alias[[nm]])
    if (!is.null(pal)) return(pal)
  }

  # Safe fallbacks
  if (nm %in% c("grays","grey","gray")) {
    return(grDevices::gray(seq(0, 1, length.out = n)))
  }

  # Generic fallback: interpolate a conservative viridis-like ramp
  base_col <- c("#440154", "#482777", "#3E4989", "#31688E",
                "#26828E", "#1F9E89", "#35B779", "#6CCE59",
                "#B4DE2C", "#FDE725")
  ramp <- grDevices::colorRampPalette(base_col, space = "Lab")
  ramp(n)
}

#' A ggplot2 fill scale with neuroimaging-friendly defaults
#'
#' @param cmap Palette name or vector of colors. See [resolve_cmap()].
#' @param range Either "robust" (quantiles) or "data" (min/max) to determine
#'   the default scale limits when `limits` is not provided.
#' @param probs Two-length numeric vector of quantiles for `range="robust"`.
#' @param limits Optional numeric limits (min, max). Overrides `range`.
#' @param na.value Color for NA.
#' @param guide Legend guide (default "colorbar").
#' @return A ggplot2 scale object.
#' @export
scale_fill_neuro <- function(
  cmap   = "grays",
  range  = c("robust","data"),
  probs  = c(.02, .98),
  limits = NULL,
  na.value = "transparent",
  guide  = "colorbar"
) {
  range <- match.arg(range)
  pal <- resolve_cmap(cmap, 256)
  ggplot2::scale_fill_gradientn(colours = pal, limits = limits,
                                na.value = na.value, guide = guide)
}


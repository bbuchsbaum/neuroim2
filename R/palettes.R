#' Neuroimaging color palettes and helpers
#'
#' Lightweight, perceptually-uniform color tools with safe fallbacks.
#' @param name Palette name (e.g., "grays", "viridis", "inferno", "magma",
#'   "plasma", "turbo", "cividis", "coldhot", "blue-red", "coolwarm"). Any
#'   palette in \code{grDevices::hcl.pals()} (e.g. "RdBu", "Spectral", "Reds")
#'   is also accepted. Case- and punctuation-insensitive. Unknown names emit a
#'   warning and fall back to a viridis-like ramp. If you pass a vector of
#'   colors, it's returned unchanged.
#' @param n Number of colors to generate.
#' @return A character vector of hex colors.
#' @export
resolve_cmap <- function(name = "grays", n = 256) {
  if (length(name) > 1L) return(name)              # already a color vector
  if (is.null(name)) name <- "grays"
  nm <- tolower(name)

  diverging <- .diverging_palettes()
  if (nm %in% names(diverging)) {
    ramp <- grDevices::colorRampPalette(diverging[[nm]], space = "Lab")
    return(ramp(n))
  }

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

  # Explicit gray fallback (covers the very common case directly)
  if (nm %in% c("grays", "grey", "gray")) {
    return(grDevices::gray(seq(0, 1, length.out = n)))
  }

  # Try hcl.colors for ANY other name. hcl.colors() matches the ~110 palettes in
  # grDevices::hcl.pals() case- and punctuation-insensitively (so "RdBu",
  # "Spectral", "Reds", "Blue-Red", ... all resolve) and errors on a true
  # mismatch, so try() cleanly distinguishes "recognized" from "unknown".
  pal <- hc_try(name)
  if (!is.null(pal) && !all(is.na(pal))) {
    return(pal)
  }

  # Unrecognized: warn (do not silently mis-color) and fall back to a
  # viridis-like ramp.
  warning(sprintf("Unknown palette '%s'; falling back to a viridis-like ramp. See grDevices::hcl.pals() for valid names.", name),
          call. = FALSE)
  base_col <- c("#440154", "#482777", "#3E4989", "#31688E",
                "#26828E", "#1F9E89", "#35B779", "#6CCE59",
                "#B4DE2C", "#FDE725")
  ramp <- grDevices::colorRampPalette(base_col, space = "Lab")
  ramp(n)
}

#' Built-in diverging palettes (anchor colors), keyed by lowercase name.
#' @keywords internal
#' @noRd
.diverging_palettes <- function() {
  list(
    coldhot    = c("#2447a8", "#59b7ff", "#111111", "#ffbf4d", "#d7191c"),
    "blue-red" = c("#3b4cc0", "#91bfdb", "#f7f7f7", "#fc8d59", "#b40426"),
    bluered    = c("#3b4cc0", "#91bfdb", "#f7f7f7", "#fc8d59", "#b40426"),
    coolwarm   = c("#3b4cc0", "#91bfdb", "#f7f7f7", "#fc8d59", "#b40426")
  )
}

#' Heuristic: is a palette name a diverging palette?
#'
#' Recognizes the package's built-in diverging palettes plus the common
#' diverging palettes shipped with \code{grDevices::hcl.pals("diverging")} and
#' \code{hcl.pals("divergingx")} (e.g. \code{"RdBu"}, \code{"Blue-Red"},
#' \code{"Spectral"}). Returns \code{FALSE} for a color vector or sequential map.
#'
#' @param name Palette name (scalar character) or color vector.
#' @keywords internal
#' @noRd
is_diverging_cmap <- function(name) {
  if (is.null(name) || length(name) != 1L || !is.character(name)) return(FALSE)
  nm <- tolower(name)
  if (nm %in% names(.diverging_palettes())) return(TRUE)
  div <- tryCatch(
    c(grDevices::hcl.pals("diverging"), grDevices::hcl.pals("divergingx")),
    error = function(e) character(0)
  )
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  norm(nm) %in% norm(div)
}

#' @keywords internal
#' @noRd
squish_oob <- function(x, range = c(0, 1), only.finite = TRUE) {
  keep <- if (only.finite) is.finite(x) else !is.na(x)
  x[keep & x < range[[1L]]] <- range[[1L]]
  x[keep & x > range[[2L]]] <- range[[2L]]
  x
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
#' @details Values outside `limits` are clamped to the nearest color endpoint
#'   rather than converted to `NA`. This keeps robust display limits from
#'   creating transparent holes in bright or dark image regions.
#' @importFrom grid unit
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
  if (identical(guide, "colorbar") || identical(guide, "colourbar")) {
    guide <- ggplot2::guide_colourbar(barheight = grid::unit(3, "cm"))
  }
  ggplot2::scale_fill_gradientn(colours = pal, limits = limits, oob = squish_oob,
                                na.value = na.value, guide = guide)
}

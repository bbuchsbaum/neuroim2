# Shared figure engine for the plot_* family.
#
# Every slice figure is built from the same pieces: a borderless dark "tile"
# per slice (cropped to the head, labelled in world coordinates, with
# orientation letters), a compact fixed-size colorbar, an optional one-line
# key, and a card that carries the title block. Keeping these pieces in one
# place is what makes plot_overlay(), plot_ortho(), plot_montage(),
# plot_edge_overlay() and plot_checkerboard() look like one family.

# ---------------------------------------------------------------------------
# Style and argument checking
# ---------------------------------------------------------------------------

#' Style tokens for the plotting family
#'
#' @param style "light", "dark" or "report".
#' @return A named list of colours.
#' @keywords internal
#' @noRd
neuro_style_tokens <- function(style = c("light", "dark", "report")) {
  style <- match.arg(style)
  base <- list(tile = "#000000", tile_fg = "#f2f2f2", tile_muted = "#9a9a9a",
               cross = "#5fd0ff", gutter = 2)
  # Separator between adjacent tiles: the page colour (black on the dark
  # page, so the tiles merge into it).
  utils::modifyList(base, switch(style,
    light  = list(card = "#ffffff", fg = "#1a1a1a", muted = "#5f5f5f",
                  faint = "#e7e7e7", seam = "#ffffff"),
    report = list(card = "#f7f7f5", fg = "#161616", muted = "#5a5a58",
                  faint = "#e6e6e2", seam = "#f7f7f5"),
    dark   = list(card = "#000000", fg = "#ededed", muted = "#a0a0a0",
                  faint = "#2c2c2c", seam = "#000000")
  ))
}

#' Match a string argument against its choices with a cli error
#' @keywords internal
#' @noRd
match_choice <- function(x, choices, arg = "style") {
  if (is.null(x) || identical(x, choices)) return(choices[[1L]])
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    cli::cli_abort("{.arg {arg}} must be one of {.val {choices}}, not {.val {x}}.",
                   call = NULL)
  }
  x
}

#' Validate overlay threshold and opacity; returns the threshold
#' @keywords internal
#' @noRd
check_overlay_args <- function(ov_thresh, ov_alpha) {
  if (is.null(ov_thresh)) ov_thresh <- 0
  if (!is.numeric(ov_thresh) || length(ov_thresh) != 1L || !is.finite(ov_thresh) || ov_thresh < 0) {
    cli::cli_abort("{.arg ov_thresh} must be a single non-negative number (values with |v| below it are hidden).",
                   call = NULL)
  }
  if (!is.numeric(ov_alpha) || length(ov_alpha) != 1L || !is.finite(ov_alpha) ||
      ov_alpha < 0 || ov_alpha > 1) {
    cli::cli_abort("{.arg ov_alpha} must be a single number between 0 and 1.", call = NULL)
  }
  ov_thresh
}

#' Validate a positive whole-number count
#' @keywords internal
#' @noRd
check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x)) {
    cli::cli_abort("{.arg {arg}} must be a positive whole number.", call = NULL)
  }
  invisible(as.integer(x))
}

#' One-time hint that slice positions are voxel indices by default
#'
#' Panels are labelled in world coordinates (mm) but, for backward
#' compatibility, \code{zlevels}/\code{coord} are voxel indices unless
#' \code{unit = "mm"}. When a caller supplies positions without choosing a unit
#' on an image whose indices differ from its mm coordinates, say so (at most
#' once every eight hours per session).
#' @keywords internal
#' @noRd
unit_hint <- function(pos, unit_missing, vol, what = "zlevels") {
  if (is.null(pos) || !unit_missing) return(invisible())
  tr <- tryCatch(trans(space(vol)), error = function(e) NULL)
  if (is.null(tr) || isTRUE(all.equal(unname(tr[1:3, 4]), c(0, 0, 0)))) return(invisible())
  cli::cli_inform(
    c("i" = "{.arg {what}} are read as voxel indices; panels are labelled in mm.",
      " " = "Pass {.code unit = \"mm\"} to give positions in world coordinates."),
    .frequency = "regularly", .frequency_id = paste0("neuroim2_unit_", what)
  )
}

#' World-coordinate bounds of a volume (rounded mm, 2 x 3)
#' @keywords internal
#' @noRd
world_bounds <- function(vol) {
  d <- dim(vol)[1:3]
  corners <- as.matrix(expand.grid(c(1, d[1]), c(1, d[2]), c(1, d[3])))
  w <- grid_to_coord(space(vol), corners)
  round(apply(w, 2L, range))
}

#' Validate a length-2 label vector
#' @keywords internal
#' @noRd
check_labels2 <- function(labels) {
  if (!is.character(labels) || length(labels) != 2L || anyNA(labels)) {
    cli::cli_abort("{.arg labels} must be a character vector of length 2.", call = NULL)
  }
  labels
}

# ---------------------------------------------------------------------------
# Slices, coordinates, layout
# ---------------------------------------------------------------------------

#' Resolve the display slice indices for a slice-grid figure
#'
#' @param zlevels User slice positions (or NULL for automatic selection).
#' @param vol Volume that defines the grid.
#' @param along Native slicing axis.
#' @param unit "index" or "mm": how to interpret \code{zlevels}.
#' @param n Number of slices for automatic selection.
#' @param support Optional logical array whose extent drives automatic
#'   selection (e.g. the supra-threshold overlay).
#' @keywords internal
#' @noRd
resolve_slice_levels <- function(zlevels, vol, along, unit = c("index", "mm"),
                                 n = 9L, support = NULL) {
  unit <- match.arg(unit)
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along) || along < 1L || along > 3L) {
    cli::cli_abort("{.arg along} must be one of 1, 2, or 3.", call = NULL)
  }
  if (is.null(zlevels)) {
    return(default_slice_levels(vol, along = along, n = n, support = support))
  }
  if (unit == "mm") {
    idx <- mm_to_slice_index(vol, zlevels, along)
    bad <- is.na(idx) | idx < 1L | idx > dim(vol)[[along]]
    if (any(bad)) {
      rng <- slice_world_range(vol, along)
      cli::cli_abort(c(
        "{.arg zlevels} must lie inside the image along this axis.",
        "i" = "Valid range is {format(round(rng[1]))} to {format(round(rng[2]))} mm; got {.val {zlevels[bad]}}."
      ), call = NULL)
    }
    zlevels <- idx
  }
  zlevels
}

#' World-coordinate range (mm) of the slices along a native axis
#' @keywords internal
#' @noRd
slice_world_range <- function(vol, along) {
  d <- dim(vol)[1:3]
  ends <- vapply(c(1, d[[along]]), function(z) {
    g <- (d + 1) / 2
    g[[along]] <- z
    w <- as.numeric(grid_to_coord(space(vol), matrix(g, nrow = 1L)))
    w[[which.max(abs(perm_mat(axes(space(vol)))[, along]))]]
  }, numeric(1))
  sort(ends)
}

#' Convert world coordinates (mm) along a slicing axis to slice indices
#' @keywords internal
#' @noRd
mm_to_slice_index <- function(vol, mm, along) {
  sp <- space(vol)
  d <- dim(vol)[1:3]
  directions <- perm_mat(axes(sp))
  world_axis <- which.max(abs(directions[, along]))
  centre <- as.numeric(grid_to_coord(sp, matrix((d + 1) / 2, nrow = 1L)))
  vapply(as.numeric(mm), function(m) {
    w <- centre
    w[[world_axis]] <- m
    g <- as.numeric(coord_to_grid(sp, matrix(w, nrow = 1L)))
    as.integer(round(g[[along]]))
  }, integer(1))
}

#' Choose an ncol that makes a grid of tiles fill the canvas
#'
#' @param n Number of tiles.
#' @param tile_aspect Tile width / height.
#' @param canvas_aspect Available width / height for the tile block.
#' @keywords internal
#' @noRd
auto_ncol <- function(n, tile_aspect = 1, canvas_aspect = 4 / 3) {
  if (n <= 1L) return(1L)
  best <- 1L; best_score <- -Inf
  for (nc in seq_len(n)) {
    nr <- ceiling(n / nc)
    grid_aspect <- (nc * tile_aspect) / nr
    # fraction of the canvas covered by tile pixels
    fill <- if (grid_aspect > canvas_aspect) canvas_aspect / grid_aspect else grid_aspect / canvas_aspect
    used <- fill * n / (nc * nr)
    # prefer layouts without a ragged last row
    score <- used - 0.08 * ((nc * nr) - n)
    if (score > best_score + 1e-9) { best <- nc; best_score <- score }
  }
  best
}

#' Size (inches) of the target canvas
#'
#' An explicit \code{canvas = c(width, height)} wins; otherwise the size of the
#' open graphics device (so knitr chunks and interactive devices are fitted
#' exactly); otherwise 10 x 7.5 in.
#' @keywords internal
#' @noRd
canvas_size <- function(canvas = NULL, default = c(10, 7.5)) {
  if (!is.null(canvas)) {
    if (!is.numeric(canvas) || length(canvas) != 2L || any(!is.finite(canvas)) ||
        any(canvas <= 0)) {
      cli::cli_abort("{.arg canvas} must be two positive numbers, c(width, height) in inches.",
                     call = NULL)
    }
    return(as.numeric(canvas))
  }
  if (grDevices::dev.cur() > 1L) {
    ds <- tryCatch(grDevices::dev.size("in"), error = function(e) NULL)
    if (!is.null(ds) && all(is.finite(ds)) && all(ds > 0)) return(ds)
  }
  default
}

#' Round a positive number up to two significant digits
#' @keywords internal
#' @noRd
nice_ceiling <- function(x) {
  if (!is.finite(x) || x <= 0) return(x)
  p <- 10^(floor(log10(x)) - 1)
  ceiling(x / p - 1e-9) * p
}

#' Format numbers for colorbar ticks (at most 3 significant digits)
#' @keywords internal
#' @noRd
format_tick <- function(x) {
  out <- formatC(signif(x, 3), format = "fg", digits = 3, flag = "#")
  out <- sub("\\.$", "", sub("(\\.[0-9]*?)0+$", "\\1", out))
  out <- sub("\\.$", "", out)
  out[abs(x) < 1e-12] <- "0"
  trimws(out)
}

# ---------------------------------------------------------------------------
# Overlay colour scale
# ---------------------------------------------------------------------------

#' Resolve the colour scale of a statistical overlay
#'
#' Computed once from all finite non-zero voxels of the \emph{whole} overlay,
#' so the same map gets the same scale whichever slices are shown.
#'
#' @param values Overlay values (whole volume).
#' @param ov_range "robust", "data" or numeric \code{c(lo, hi)}.
#' @param thresh Threshold (>= 0).
#' @param ov_cmap Palette name/vector or NULL (automatic).
#' @param ov_symmetric NULL (auto), TRUE or FALSE.
#' @param ov_cap Optional magnitude cap.
#' @return list(lim, pal, cmap, diverging, signed, over_hi, over_lo)
#' @keywords internal
#' @noRd
overlay_scale <- function(values, ov_range = "robust", probs = c(.02, .98),
                          thresh = 0, ov_cmap = NULL, ov_symmetric = NULL,
                          ov_cap = NULL) {
  fin <- values[is.finite(values)]
  signed <- length(fin) > 0L && min(fin) < 0 && max(fin) > 0
  symmetric <- if (is.null(ov_symmetric)) signed else isTRUE(ov_symmetric)
  lim <- overlay_display_limits(ov_range, values, probs = probs, thresh = thresh)
  auto_range <- !is.numeric(ov_range)
  cap <- if (!is.null(ov_cap)) abs(ov_cap) else max(abs(lim), na.rm = TRUE)
  if (auto_range && is.null(ov_cap)) cap <- nice_ceiling(cap)
  if (isTRUE(thresh > 0) && is.finite(cap) && cap <= thresh) {
    cli::cli_warn(c("The overlay colour scale does not exceed {.arg ov_thresh}.",
                    "i" = "Extending it to {format_tick(thresh * 1.25)} so the threshold is visible."))
    cap <- thresh * 1.25
  }
  if (symmetric && is.finite(cap) && cap > 0) {
    lim <- c(-cap, cap)
  } else if (!is.null(ov_cap) || auto_range) {
    lim[2] <- cap
    if (auto_range && lim[1] > 0) lim[1] <- 0
  }
  if (is.null(ov_cmap)) {
    ov_cmap <- if (symmetric) "cold_hot" else "hot"
  } else if (signed && length(ov_cmap) == 1L && !is_diverging_cmap(ov_cmap)) {
    cli::cli_warn(c("The overlay has positive and negative values but {.arg ov_cmap} is not diverging.",
                    "i" = "Negative values may render poorly; consider {.val cold_hot} or {.val RdBu}."))
  }
  list(lim = lim, pal = resolve_cmap(ov_cmap, 256), cmap = ov_cmap,
       diverging = symmetric, signed = signed,
       over_hi = length(fin) > 0L && max(fin) > lim[2] * (1 + 1e-9),
       over_lo = length(fin) > 0L && min(fin) < lim[1] * (1 + 1e-9) && lim[1] < 0)
}

#' Map overlay values to palette positions (threshold-aware)
#'
#' Visible values are spread over the whole visible part of the palette: with a
#' threshold, \code{|v| = thresh} maps to the first saturated colour (the
#' neutral centre of a diverging palette, or the near-black start of a
#' sequential one, is skipped) and the cap maps to the palette end. This gives
#' clusters internal gradation instead of flat saturated blobs.
#'
#' @param v Numeric values.
#' @param lim Display limits (symmetric when \code{diverging}).
#' @param thresh Threshold (0 for none).
#' @param diverging Logical; treat the palette as two-sided.
#' @param skip Fraction of each palette arm skipped at the threshold.
#' @return Palette positions in [0, 1] (NA for sub-threshold / missing).
#' @keywords internal
#' @noRd
overlay_positions <- function(v, lim, thresh = 0, diverging = TRUE, skip = 0.3) {
  thresh <- if (isTRUE(thresh > 0)) thresh else 0
  sk <- if (thresh > 0) skip else 0
  pos <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  if (diverging) {
    cap <- max(abs(lim))
    m <- abs(v)
    ok <- ok & m >= thresh
    den <- max(cap - thresh, .Machine$double.eps)
    s <- pmin(pmax((m[ok] - thresh) / den, 0), 1)
    pos[ok] <- 0.5 + sign(v[ok]) * 0.5 * (sk + (1 - sk) * s)
  } else {
    lo <- if (thresh > 0) max(lim[1], thresh) else lim[1]
    hi <- lim[2]
    ok <- ok & (if (thresh > 0) abs(v) >= thresh else TRUE)
    den <- max(hi - lo, .Machine$double.eps)
    s <- pmin(pmax((v[ok] - lo) / den, 0), 1)
    pos[ok] <- sk + (1 - sk) * s
  }
  pos
}

#' Default soft-alpha floor: 60\% when a threshold is set, else none
#' @keywords internal
#' @noRd
default_alpha_floor <- function(thresh) if (isTRUE(thresh > 0)) 0.6 else 0

#' Build the per-voxel opacity function for an overlay
#'
#' Returns a function of absolute value giving alpha in [0, 1]. The same
#' function shades the colorbar, so the key always matches the picture. With a
#' threshold, every supra-threshold voxel keeps at least \code{floor} opacity
#' in the graded modes, so nothing that passes the threshold disappears.
#' @keywords internal
#' @noRd
overlay_alpha_fun <- function(mode, thresh, cap, soft = NULL, floor = 0.6) {
  thresh <- if (isTRUE(thresh > 0)) thresh else 0
  cap <- if (is.finite(cap) && cap > 0) cap else 1
  fl <- if (thresh > 0) floor else 0
  graded <- function(t) {
    t <- pmin(pmax(t, 0), 1)
    fl + (1 - fl) * t
  }
  switch(mode,
    binary = function(m) ifelse(is.finite(m) & m >= thresh & m > 0, 1, 0),
    proportional = function(m) {
      a <- if (thresh > 0) graded((m - thresh) / max(cap - thresh, 1e-12)) else pmin(m / cap, 1)
      ifelse(is.finite(m) & m >= thresh & m > 0, a, 0)
    },
    ramp = function(m) {
      a <- graded((m - thresh) / max(cap - thresh, 1e-12))
      if (thresh == 0) a <- pmin(pmax(m / cap, 0), 1)
      ifelse(is.finite(m) & m >= thresh & m > 0, a, 0)
    },
    soft = function(m) {
      t <- (m - soft$lo) / max(soft$hi - soft$lo, 1e-12)
      sf <- if (is.null(soft$alpha_floor)) fl else soft$alpha_floor
      a <- sf + (1 - sf) * pmin(pmax(t, 0), 1)^soft$gamma
      ifelse(is.finite(m) & m >= max(thresh, soft$lo) & m > 0, a, 0)
    }
  )
}

#' Build an RGBA raster array from palette positions
#' @keywords internal
#' @noRd
positions_to_rgba <- function(pos, nr, nc, cols, alpha = 1, alpha_map = NULL) {
  idx <- 1L + floor(pos * (length(cols) - 1L))
  idx[!is.finite(idx)] <- 1L
  idx <- pmax(1L, pmin(length(cols), idx))
  rgb <- grDevices::col2rgb(cols[idx]) / 255
  eff <- if (is.null(alpha_map)) rep(alpha, length(pos)) else as.numeric(alpha_map) * alpha
  eff[!is.finite(eff) | !is.finite(pos)] <- 0
  eff <- pmax(0, pmin(1, eff))
  out <- array(0, dim = c(nr, nc, 4L))
  out[, , 1] <- rgb[1, ]
  out[, , 2] <- rgb[2, ]
  out[, , 3] <- rgb[3, ]
  out[, , 4] <- eff
  out
}

#' Raster grob for one overlay slice
#' @keywords internal
#' @noRd
overlay_slice_grob <- function(overlay, z, along, scale, thresh, alpha_fun,
                               alpha = 1, downsample = 1L) {
  mov <- volume_slice_matrix(overlay, z, along = along)
  amap <- matrix(alpha_fun(abs(mov)), nrow(mov), ncol(mov))
  o <- orient_volume_slice_for_raster(overlay, z, along = along, mat = mov,
                                      alpha_map = amap, downsample = downsample)
  pos <- overlay_positions(as.numeric(o$mat), scale$lim, thresh = thresh,
                           diverging = scale$diverging)
  grid::rasterGrob(positions_to_rgba(pos, nrow(o$mat), ncol(o$mat), scale$pal,
                                     alpha = alpha, alpha_map = o$alpha_map),
                   interpolate = FALSE)
}

# ---------------------------------------------------------------------------
# Tiles
# ---------------------------------------------------------------------------

#' A single borderless slice tile
#'
#' @param oriented Output of \code{orient_volume_slice_for_raster()} for the
#'   background.
#' @param bg_lim,bg_cmap Background display limits and palette.
#' @param layers List of grobs (overlay rasters) drawn over the background.
#' @param window Optional \code{list(xlim, ylim)} crop window.
#' @param label Slice label drawn inside the tile (top-left), or NULL.
#' @param orient Character vector of orientation letters
#'   (\code{left/right/top/bottom}) to draw, or NULL.
#' @param tokens Style tokens.
#' @param interpolate Smooth the background raster.
#' @keywords internal
#' @noRd
neuro_tile <- function(oriented, bg_lim, bg_cmap = "grays", layers = list(),
                       window = NULL, label = NULL, orient = NULL,
                       tokens = neuro_style_tokens("light"),
                       interpolate = FALSE, label_size = 3.1, extra = list(),
                       floor_y = -Inf) {
  df <- oriented_raster_df(oriented)
  xr <- raster_extent_from_centers(oriented$x)
  yr <- raster_extent_from_centers(oriented$y)
  if (is.null(window)) window <- list(xlim = xr, ylim = yr)
  # The tile field (and any padding beyond the image) takes the palette's
  # lowest colour, so it is continuous with the image's own background
  # whatever the colour map.
  tokens$tile <- resolve_cmap(bg_cmap, 2L)[[1L]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_raster(interpolate = interpolate) +
    scale_fill_neuro(cmap = bg_cmap, limits = bg_lim, guide = "none",
                     na.value = tokens$tile)
  for (g in layers) {
    p <- p + ggplot2::annotation_custom(g, xmin = xr[1], xmax = xr[2],
                                        ymin = yr[1], ymax = yr[2])
  }
  for (e in extra) p <- p + e
  p <- p +
    ggplot2::coord_fixed(xlim = window$xlim, ylim = window$ylim, expand = FALSE) +
    neuro_tile_theme(tokens)
  add_tile_text(p, window, label = label, orient = orient, tokens = tokens,
                label_size = label_size, floor_y = floor_y)
}

#' Theme for a borderless tile
#' @keywords internal
#' @noRd
neuro_tile_theme <- function(tokens) {
  g <- tokens$gutter / 2
  ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = tokens$tile, colour = NA),
      plot.background  = ggplot2::element_rect(fill = tokens$seam, colour = NA),
      plot.margin = grid::unit(c(g, g, g, g), "pt"),
      legend.position = "none"
    )
}

#' Text with a dark halo so it stays legible over bright tissue
#' @keywords internal
#' @noRd
halo_text <- function(x, y, label, hjust, vjust, colour, halo, size,
                      fontface = "plain", offset) {
  dx <- offset * c(-1, 1, 0, 0, -0.7, 0.7, -0.7, 0.7)
  dy <- offset * c(0, 0, -1, 1, -0.7, -0.7, 0.7, 0.7)
  list(
    ggplot2::annotate("text", x = rep(x, each = 8L) + rep(dx, length(x)),
                      y = rep(y, each = 8L) + rep(dy, length(y)),
                      label = rep(label, each = 8L), hjust = rep(hjust, each = 8L),
                      vjust = rep(vjust, each = 8L), colour = halo, size = size,
                      fontface = fontface, alpha = 0.85),
    ggplot2::annotate("text", x = x, y = y, label = label, hjust = hjust,
                      vjust = vjust, colour = colour, size = size, fontface = fontface)
  )
}

#' Add the in-tile slice label and orientation letters
#' @keywords internal
#' @noRd
add_tile_text <- function(p, window, label = NULL, orient = NULL, tokens,
                          label_size = 3.1, floor_y = -Inf) {
  xw <- diff(window$xlim); yw <- diff(window$ylim)
  s <- max(xw, yw)
  pad <- 0.035 * s
  off <- 0.0045 * s
  if (!is.null(label) && nzchar(label)) {
    p <- p + halo_text(window$xlim[1] + pad, window$ylim[2] - pad, label, 0, 1,
                       tokens$tile_fg, tokens$tile, label_size, offset = off)
  }
  if (!is.null(orient) && length(orient)) {
    pos <- list(
      left   = c(window$xlim[1] + pad * 0.8, mean(window$ylim), 0, 0.5),
      right  = c(window$xlim[2] - pad * 0.8, mean(window$ylim), 1, 0.5),
      top    = c(mean(window$xlim), window$ylim[2] - pad * 0.6, 0.5, 1),
      # Keep the inferior letter with the image when the view is padded
      # below the field of view.
      bottom = c(mean(window$xlim), max(window$ylim[1], floor_y) + pad * 0.6, 0.5, 0)
    )
    sides <- intersect(names(orient), names(pos))
    q <- do.call(rbind, pos[sides])
    p <- p + halo_text(q[, 1], q[, 2], unname(orient[sides]), q[, 3], q[, 4],
                       tokens$tile_fg, tokens$tile, label_size * 1.1,
                       fontface = "bold", offset = off)
  }
  p
}

#' Orientation letters to draw on a tile
#'
#' Axial and coronal tiles get L/R (the convention that matters most for
#' lateralised results); sagittal tiles get A/P.
#' @keywords internal
#' @noRd
orientation_letters <- function(oriented, all = FALSE) {
  lab <- oriented$labels
  if (isTRUE(all)) return(lab)
  lab[c("left", "right")]
}

#' Width / height of a tile after cropping
#' @keywords internal
#' @noRd
tile_aspect_of <- function(window, oriented) {
  if (is.null(window)) {
    window <- list(xlim = raster_extent_from_centers(oriented$x),
                   ylim = raster_extent_from_centers(oriented$y))
  }
  diff(window$xlim) / diff(window$ylim)
}

#' Crop window around the head (and any supra-threshold overlay)
#'
#' @param bgvol Background volume.
#' @param zlevels,along Slices on display.
#' @param extra Optional list of volumes whose non-zero (above
#'   \code{extra_thresh}) voxels must stay inside the window.
#' @param margin Fractional padding.
#' @keywords internal
#' @noRd
foreground_crop_window <- function(bgvol, zlevels, along, extra = list(),
                                   extra_thresh = 0, margin = 0.07,
                                   bg_thresh = NULL) {
  if (is.null(bg_thresh)) {
    vals <- unlist(lapply(zlevels, function(z) as.numeric(volume_slice_matrix(bgvol, z, along))))
    bg_thresh <- foreground_threshold(vals)
  }
  if (!is.finite(bg_thresh)) return(NULL)
  xs <- NULL; ys <- NULL; ext_x <- NULL; ext_y <- NULL
  add <- function(o, keep) {
    ext_x <<- range(c(ext_x, raster_extent_from_centers(o$x)))
    ext_y <<- range(c(ext_y, raster_extent_from_centers(o$y)))
    if (!length(keep)) return(invisible())
    nx <- length(o$x)
    xi <- ((keep - 1L) %% nx) + 1L
    yi <- ((keep - 1L) %/% nx) + 1L
    xs <<- range(c(xs, o$x[xi]))
    ys <<- range(c(ys, o$y[yi]))
  }
  for (z in zlevels) {
    o <- orient_volume_slice_for_raster(bgvol, z, along = along)
    v <- c(t(o$mat))
    add(o, which(is.finite(v) & v > bg_thresh))
    for (ev in extra) {
      oe <- orient_volume_slice_for_raster(ev, z, along = along)
      ve <- c(t(oe$mat))
      add(oe, which(is.finite(ve) & abs(ve) > extra_thresh))
    }
  }
  if (length(xs) < 2L || length(ys) < 2L || diff(xs) == 0 || diff(ys) == 0) return(NULL)
  pad <- max(diff(xs), diff(ys)) * margin
  list(xlim = c(max(xs[1] - pad, ext_x[1]), min(xs[2] + pad, ext_x[2])),
       ylim = c(max(ys[1] - pad, ext_y[1]), min(ys[2] + pad, ext_y[2])))
}

#' Gapped crosshair that stops short of the tile edges
#' @keywords internal
#' @noRd
add_crosshair <- function(p, at, window, tokens, gap = 0.045, edge = 0.12) {
  xw <- diff(window$xlim); yw <- diff(window$ylim)
  s <- max(xw, yw)
  g <- gap * s; e <- edge * s
  segs <- data.frame(
    x    = c(window$xlim[1] + e, at[1] + g, at[1], at[1]),
    xend = c(at[1] - g, window$xlim[2] - e, at[1], at[1]),
    y    = c(at[2], at[2], window$ylim[1] + e, at[2] + g),
    yend = c(at[2], at[2], at[2] - g, window$ylim[2] - e)
  )
  segs <- segs[(segs$xend - segs$x) >= 0 & (segs$yend - segs$y) >= 0, , drop = FALSE]
  p + ggplot2::annotate("segment", x = segs$x, xend = segs$xend, y = segs$y,
                        yend = segs$yend, colour = tokens$cross, alpha = 0.8,
                        linewidth = 0.4)
}

# ---------------------------------------------------------------------------
# Colorbar and key
# ---------------------------------------------------------------------------

#' Compact fixed-size colorbar
#'
#' @param lim Display limits.
#' @param cols Palette colours.
#' @param thresh Threshold (neutral band and ticks when > 0).
#' @param diverging Whether the overlay mapping is two-sided.
#' @param positions Logical; use \code{overlay_positions()} for colours (overlay
#'   bars) rather than a plain linear ramp.
#' @param alpha_fun Optional opacity function of |value| (overlay bars), so the
#'   bar fades exactly as the overlay does.
#' @param over_hi,over_lo Logical; data extend beyond the upper / lower limit
#'   (the end tick is then prefixed with a comparison sign).
#' @param tokens Style tokens.
#' @param title Title drawn above the bar.
#' @keywords internal
#' @noRd
neuro_colorbar <- function(lim, cols, thresh = 0, diverging = FALSE,
                           positions = TRUE, tokens, title = NULL,
                           alpha_fun = NULL, alpha = 1, over_hi = FALSE,
                           over_lo = FALSE) {
  n <- 256L
  edges <- seq(lim[1], lim[2], length.out = n + 1L)
  mids <- (edges[-1] + edges[-(n + 1L)]) / 2
  pos <- if (isTRUE(positions)) {
    overlay_positions(mids, lim, thresh = thresh, diverging = diverging)
  } else {
    (mids - lim[1]) / diff(lim)
  }
  idx <- 1L + floor(pos * (length(cols) - 1L))
  fill <- ifelse(is.finite(idx), cols[pmax(1L, pmin(length(cols), idx))], tokens$faint)
  if (!is.null(alpha_fun)) {
    # Fade toward a tissue grey, not the page: faded voxels in the figure sit
    # on brain, so this is the colour they actually take on.
    a <- alpha_fun(abs(mids)) * alpha
    ok <- is.finite(idx)
    fill[ok] <- blend_colours(fill[ok], "#b8b8b8", a[ok])
  }
  # The sub-threshold band is drawn as a narrow neutral strip inside the bar
  # outline, so it reads as "not shown" rather than as a colour.
  band <- !is.finite(idx)
  df <- data.frame(ymin = edges[-(n + 1L)], ymax = edges[-1], fill = fill,
                   xmin = ifelse(band, 0.42, 0), xmax = ifelse(band, 0.58, 1))

  breaks <- if (isTRUE(thresh > 0)) {
    b <- if (diverging) c(lim[1], -thresh, thresh, lim[2]) else c(max(lim[1], thresh), lim[2])
    if (!diverging && lim[1] < thresh) b <- c(lim[1], b)
    sort(unique(b[b >= lim[1] & b <= lim[2]]))
  } else {
    # Always label both ends of the bar (the window limits), plus round
    # interior values that are not crowded against the ends.
    inner <- pretty(lim, n = 4)
    inner <- inner[inner > lim[1] + 0.12 * diff(lim) & inner < lim[2] - 0.12 * diff(lim)]
    b <- sort(unique(c(lim[1], inner, lim[2])))
    if (diverging && !any(abs(b) < 1e-12)) b <- sort(c(b, 0))
    b
  }
  labels <- format_tick(breaks)
  if (!isTRUE(thresh > 0) && !isTRUE(positions)) {
    # Window limits are arbitrary numbers; show them at two significant
    # digits so they sit comfortably beside the round interior ticks.
    # Round inward, so the label never claims more range than the window has.
    lo_i <- which.min(breaks); hi_i <- which.max(breaks)
    p10 <- function(x) 10^(floor(log10(abs(x))) - 1)
    lo_v <- breaks[lo_i]; hi_v <- breaks[hi_i]
    if (lo_v != 0) labels[lo_i] <- format_tick(ceiling(lo_v / p10(lo_v)) * p10(lo_v))
    if (hi_v != 0) labels[hi_i] <- format_tick(floor(hi_v / p10(hi_v)) * p10(hi_v))
  }
  # End ticks read ">= cap" / "<= -cap" when the data extend beyond the scale.
  # plotmath keeps the source ASCII and renders the signs on every device
  # (a literal Unicode sign fails on the base pdf() device).
  mark <- rep("", length(labels))
  if (isTRUE(over_hi) && length(breaks)) {
    top <- which.max(breaks)
    if (isTRUE(all.equal(breaks[top], lim[2]))) mark[top] <- ">="
  }
  if (isTRUE(over_lo) && length(breaks)) {
    bot <- which.min(breaks)
    if (isTRUE(all.equal(breaks[bot], lim[1]))) mark[bot] <- "<="
  }
  if (any(nzchar(mark))) {
    labels <- parse(text = ifelse(nzchar(mark),
                                  sprintf("phantom() %s '%s'", mark, labels),
                                  sprintf("'%s'", labels)))
  }

  # Single-symbol quantities (t, z, F, r) are set in italic, as in journals.
  short_title <- !is.null(title) && nchar(title) <= 2L
  ggplot2::ggplot(df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = df$fill, colour = NA) +
    ggplot2::scale_y_continuous(position = "right", breaks = breaks,
                                labels = labels,
                                expand = ggplot2::expansion(0)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      # Transparent: the card behind it is painted by the figure, and a bar
      # can then never paint over neighbouring tiles.
      plot.background = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_text(colour = tokens$fg, size = 8.5,
                                                margin = ggplot2::margin(l = 3)),
      axis.ticks.y.right = ggplot2::element_line(colour = tokens$fg, linewidth = 0.3),
      axis.ticks.length.y.right = grid::unit(3, "pt"),
      plot.title = ggplot2::element_text(
        colour = tokens$fg, size = if (short_title) 10 else 9, hjust = 0,
        face = if (short_title) "italic" else "plain",
        margin = ggplot2::margin(b = 5)),
      plot.title.position = "plot",
      plot.margin = grid::unit(c(2, 2, 2, 6), "pt")
,
      # A fixed physical bar width, whatever the label lengths.
      panel.widths = grid::unit(0.3, "in")
    )
}

#' Blend colours over a background at the given opacities
#' @keywords internal
#' @noRd
blend_colours <- function(cols, bg, alpha) {
  alpha <- pmin(pmax(alpha, 0), 1)
  fg <- grDevices::col2rgb(cols) / 255
  b <- as.numeric(grDevices::col2rgb(bg) / 255)
  out <- fg * rep(alpha, each = 3L) + b * (1 - rep(alpha, each = 3L))
  grDevices::rgb(out[1, ], out[2, ], out[3, ])
}

#' One-line key strip (swatches + text), packed from the left
#'
#' @param items List of \code{list(fill = colour, label = text)}; \code{fill}
#'   may be NA for a text-only item.
#' @param tokens Style tokens.
#' @keywords internal
#' @noRd
neuro_key_strip <- function(items, tokens, fontsize = 8.5) {
  # All text shares one gpar set on the parent tree, so grobWidth() is
  # evaluated in the font the text is drawn in and positions stay exact.
  # Items are separated by a short vertical rule (no font metrics involved).
  grobs <- list()
  x <- grid::unit(2, "pt")
  for (k in seq_along(items)) {
    it <- items[[k]]
    if (k > 1L) {
      grobs <- c(grobs, list(grid::segmentsGrob(
        x0 = x - grid::unit(13, "pt"), x1 = x - grid::unit(13, "pt"),
        y0 = grid::unit(0.5, "npc") - grid::unit(4.5, "pt"),
        y1 = grid::unit(0.5, "npc") + grid::unit(4.5, "pt"),
        gp = grid::gpar(col = tokens$faint, lwd = 1))))
    }
    if (!is.null(it$fill) && !is.na(it$fill)) {
      grobs <- c(grobs, list(grid::rectGrob(
        x = x, y = 0.5, width = grid::unit(9, "pt"), height = grid::unit(9, "pt"),
        just = c("left", "centre"), gp = grid::gpar(fill = it$fill, col = NA))))
      x <- x + grid::unit(13, "pt")
    }
    tg <- grid::textGrob(it$label, x = x, y = 0.5, just = c("left", "centre"))
    grobs <- c(grobs, list(tg))
    x <- x + grid::grobWidth(tg) + grid::unit(26, "pt")
  }
  tree <- do.call(grid::grobTree, c(grobs, list(gp = grid::gpar(col = tokens$muted, fontsize = fontsize))))
  patchwork::wrap_elements(full = tree) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = tokens$card, colour = NA),
                   plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
}

#' Compact key for an overlay figure
#' @keywords internal
#' @noRd
overlay_key <- function(thresh, diverging, quantity, tokens, plane = "Axial") {
  items <- list(list(fill = NA, label = sprintf("%s slices, neurological view (L = left)", plane)))
  if (isTRUE(thresh > 0)) {
    q <- if (is.null(quantity) || identical(quantity, "value")) "value" else quantity
    thr <- format_tick(thresh)
    # Single-symbol quantities are italic; |.| is drawn with plotmath group()
    # delimiters, spaced so the bars read as bars rather than as the letter l.
    qe <- if (nchar(q) <= 2L) bquote(italic(.(q))) else q
    lab <- if (diverging) bquote(group("|", .(qe), "|") ~ phantom() >= .(thr) ~ "shown") else bquote(.(qe) >= .(thr) ~ "shown")
    items <- c(items, list(list(fill = NA, label = lab)))
  }
  neuro_key_strip(items, tokens)
}

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

#' Assemble tiles + colorbar + key into a titled figure
#'
#' The layout is fitted to the target canvas (see \code{canvas_size()}): the
#' tile block is sized to fill it at the tiles' aspect ratio, the colorbar sits
#' immediately right of the tiles, and title, key and caption are
#' left-aligned with the tiles. The whole block is centred.
#'
#' @param plots List of tile ggplots.
#' @param ncol Columns (NULL = automatic from tile and canvas aspect).
#' @param tile_aspect Tile width / height (for automatic ncol).
#' @param colorbar A colorbar ggplot or NULL.
#' @param key A key strip or NULL.
#' @param widths Optional relative widths of the tiles (single-row layouts),
#'   in units of the common tile height.
#' @keywords internal
#' @noRd
neuro_assemble <- function(plots, ncol = NULL, tile_aspect = 1, colorbar = NULL,
                           key = NULL, tokens, title = NULL, subtitle = NULL,
                           caption = NULL, widths = NULL, canvas = NULL) {
  n <- length(plots)
  if (!is.null(ncol)) {
    ok <- is.numeric(ncol) && length(ncol) == 1L && is.finite(ncol) &&
      ncol >= 1 && ncol == round(ncol)
    if (!ok) {
      cli::cli_abort("{.arg ncol} must be a positive integer (or NULL for automatic).",
                     call = NULL)
    }
  }
  dev <- canvas_size(canvas)
  cb_in <- if (is.null(colorbar)) 0 else 0.8
  key_in <- if (is.null(key)) 0 else 0.28
  head_in <- title_block_height(title, subtitle)
  cap_in <- if (is.null(caption)) 0 else 0.28
  margin_in <- 0.25
  avail_w <- max(dev[1] - 2 * margin_in - cb_in, 1)
  avail_h <- max(dev[2] - 2 * margin_in - head_in - key_in - cap_in, 1)

  if (!is.null(widths)) {
    ncol <- n
    nrow <- 1L
    grid_aspect <- sum(widths)
  } else {
    if (is.null(ncol)) ncol <- auto_ncol(n, tile_aspect, avail_w / avail_h)
    ncol <- max(1L, min(as.integer(ncol), n))
    nrow <- ceiling(n / ncol)
    grid_aspect <- ncol * tile_aspect / nrow
  }
  block_w <- min(avail_w, avail_h * grid_aspect)
  block_h <- block_w / grid_aspect
  spare_w <- max(avail_w - block_w, 0)
  spare_h <- max(avail_h - block_h, 0)

  # Scale in-tile text (slice labels, orientation letters) with the tile
  # height, so small figures (a 3.5 in journal column) stay legible.
  tile_h <- block_h / nrow
  k <- min(1, max(0.8, tile_h / 2.6))
  if (k < 0.999) plots <- lapply(plots, scale_tile_text, k = k)
  grid_block <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow, widths = widths)

  # Main row: tiles plus the colorbar immediately to their right.
  main <- grid_block
  main_w <- block_w
  if (!is.null(colorbar)) {
    bar_h <- min(block_h * if (nrow == 1L) 0.7 else 0.5, 3)
    pad_h <- (block_h - bar_h) / 2
    # wrap_elements() keeps the colorbar's title and tick labels out of
    # patchwork's panel alignment, which would otherwise shrink the tiles.
    cb_col <- patchwork::wrap_plots(card_spacer(tokens),
                                    patchwork::wrap_elements(full = colorbar) +
                                      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"),
                                                     plot.background = ggplot2::element_blank()),
                                    card_spacer(tokens),
                                    ncol = 1L, heights = c(pad_h, bar_h, pad_h))
    main <- patchwork::wrap_plots(grid_block, cb_col, nrow = 1L,
                                  widths = grid::unit.c(grid::unit(1, "null"), grid::unit(cb_in, "in")))
    main_w <- block_w + cb_in
  }

  # Column: title block, main row, key, caption -- all left-aligned with the
  # tiles. Non-tile rows use absolute heights: with fixed-aspect tiles,
  # patchwork collapses rows sized in relative units.
  parts <- list(); heights <- list()
  if (head_in > 0) {
    parts <- c(parts, list(title_element(title, subtitle, tokens)))
    heights <- c(heights, list(grid::unit(head_in, "in")))
  }
  parts <- c(parts, list(main)); heights <- c(heights, list(grid::unit(1, "null")))
  if (!is.null(key)) {
    parts <- c(parts, list(key)); heights <- c(heights, list(grid::unit(key_in, "in")))
  }
  if (!is.null(caption)) {
    parts <- c(parts, list(caption_element(caption, tokens)))
    heights <- c(heights, list(grid::unit(cap_in, "in")))
  }
  column <- if (length(parts) > 1L) {
    patchwork::wrap_plots(parts, ncol = 1L, heights = do.call(grid::unit.c, heights))
  } else {
    main
  }

  # Centre the column on the canvas. Patchwork reserves absolute sizes
  # (colorbar width, title/key/caption heights) before sharing out relative
  # ones, so the spacers are weighed against the relative part only.
  if (spare_w > 0.02) {
    column <- patchwork::wrap_plots(card_spacer(tokens), column, card_spacer(tokens),
                                    nrow = 1L, widths = c(spare_w / 2, block_w, spare_w / 2))
  }
  if (spare_h > 0.02) {
    # Keep the title block tight to the tiles: vertical slack goes mostly
    # below the figure (a little above, so it does not hug the top edge).
    column <- patchwork::wrap_plots(card_spacer(tokens), column, card_spacer(tokens),
                                    ncol = 1L,
                                    heights = c(spare_h * 0.3, block_h, spare_h * 0.7))
  }
  column + patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = tokens$card, colour = NA),
      plot.margin = grid::unit(rep(margin_in * 72, 4), "pt")
    )
  )
}

#' Height (inches) of the title block: line heights plus one fixed gap
#' @keywords internal
#' @noRd
title_block_height <- function(title, subtitle) {
  h <- 0.25 * (!is.null(title)) + 0.2 * (!is.null(subtitle))
  if (h > 0) h + 0.14 else 0
}

#' Left-aligned title / subtitle block as a patchwork element
#'
#' Anchored to the bottom of its cell so the gap to the tiles is fixed.
#' @keywords internal
#' @noRd
title_element <- function(title, subtitle, tokens) {
  grobs <- list()
  base <- grid::unit(0.14, "in")
  if (!is.null(subtitle)) {
    grobs <- c(grobs, list(grid::textGrob(
      subtitle, x = 0, y = base, hjust = 0, vjust = 0,
      gp = grid::gpar(col = tokens$muted, fontsize = 10.5))))
  }
  if (!is.null(title)) {
    y <- if (is.null(subtitle)) base else base + grid::unit(0.2, "in")
    grobs <- c(grobs, list(grid::textGrob(
      title, x = 0, y = y, hjust = 0, vjust = 0,
      gp = grid::gpar(col = tokens$fg, fontsize = 14, fontface = "bold"))))
  }
  patchwork::wrap_elements(full = do.call(grid::grobTree, grobs)) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = tokens$card, colour = NA),
                   plot.margin = grid::unit(c(0, 0, 0, tokens$gutter / 2), "pt"))
}

#' Scale the size of every text layer of a tile (without touching the
#' original layers, which are reused when the figure is re-fitted)
#' @keywords internal
#' @noRd
scale_tile_text <- function(p, k) {
  p$layers <- lapply(p$layers, function(l) {
    if (!inherits(l$geom, "GeomText") || is.null(l$aes_params$size)) return(l)
    l2 <- ggplot2::ggproto(NULL, l)
    l2$aes_params <- utils::modifyList(l$aes_params, list(size = l$aes_params$size * k))
    l2
  })
  p
}

#' Left-aligned caption as a patchwork element
#' @keywords internal
#' @noRd
caption_element <- function(caption, tokens) {
  patchwork::wrap_elements(full = grid::textGrob(
    caption, x = 0, y = 0.5, hjust = 0, vjust = 0.5,
    gp = grid::gpar(col = tokens$muted, fontsize = 8.5))) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = tokens$card, colour = NA),
                   plot.margin = grid::unit(c(0, 0, 0, 2), "pt"))
}

#' A blank, transparent patchwork cell
#' @keywords internal
#' @noRd
card_spacer <- function(tokens) {
  patchwork::plot_spacer() +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
}

#' Finish a figure according to draw/assemble
#'
#' Follows the ggplot idiom: the figure is returned visibly (so it prints at
#' the console or in knitr, and can be passed to \code{ggsave()}). With
#' \code{draw = TRUE} it is also printed immediately and returned invisibly.
#' @keywords internal
#' @noRd
neuro_finish <- function(fig, plots, draw, assemble, title, subtitle, caption,
                         style, panel_names = NULL) {
  if (!is.null(panel_names) && is.null(names(plots))) names(plots) <- panel_names
  attr(plots, "labels") <- list(title = title, subtitle = subtitle, caption = caption)
  if (isTRUE(assemble)) {
    if (isTRUE(draw)) {
      print(fig)
      return(invisible(fig))
    }
    return(fig)
  }
  if (!isTRUE(draw)) return(invisible(plots))
  draw_plot_panel_grid(plots, ncol = min(3L, length(plots)), title = title,
                       subtitle = subtitle, caption = caption,
                       style = if (identical(style, "dark")) "dark" else "light")
}

# ---------------------------------------------------------------------------
# Device-fitted figures
# ---------------------------------------------------------------------------

#' Wrap a figure so its layout is re-fitted to the device when drawn
#'
#' The tiles are built once; only the arrangement (columns, centring,
#' colorbar placement) depends on the canvas. \code{build(canvas)} recreates
#' that arrangement for a canvas size in inches. Printing (console, knitr) and
#' \code{ggsave()} (via \code{grid.draw()}) call it with the size of the
#' device actually drawn on, so \code{p <- plot_overlay(...);
#' ggsave("f.png", p, width = 6, height = 9)} is laid out for 6 x 9 in. If the
#' figure has been modified after creation (e.g. \code{p + theme(...)}), it is
#' drawn as-is so the modification is kept.
#' @keywords internal
#' @noRd
neuro_figure <- function(fig, build, canvas, tokens = NULL) {
  env <- new.env(parent = emptyenv())
  env$tokens <- tokens
  env$fig <- fig
  env$build <- build
  env$canvas <- canvas
  attr(fig, "neuro_figure") <- env
  class(fig) <- unique(c("neuro_fig", class(fig)))
  fig
}

#' Resolve a device-fitted figure to a plain patchwork for drawing
#' @keywords internal
#' @noRd
neuro_fig_resolve <- function(x) {
  env <- attr(x, "neuro_figure")
  y <- neuro_fig_strip(x)
  if (!is.environment(env) || grDevices::dev.cur() <= 1L) return(y)
  if (!identical(y, env$fig)) return(y)
  ds <- tryCatch(grDevices::dev.size("in"), error = function(e) NULL)
  if (is.null(ds) || any(!is.finite(ds)) || any(ds <= 0)) return(y)
  if (all(abs(ds - env$canvas) < 0.01 * env$canvas)) return(y)
  env$build(ds)
}

#' Strip the neuro_fig wrapper from a figure
#' @keywords internal
#' @noRd
neuro_fig_strip <- function(x) {
  attr(x, "neuro_figure") <- NULL
  class(x) <- setdiff(class(x), "neuro_fig")
  x
}

#' Add a ggplot/patchwork component to a device-fitted figure
#'
#' The component is applied now and recorded, so that when the figure is
#' re-fitted to another canvas the same additions are replayed on the new
#' arrangement (\code{p + labs(title = "x")} then \code{ggsave(width = 6)}
#' keeps both the title and a fitted layout).
#' @keywords internal
#' @noRd
neuro_fig_add <- function(e1, e2, op) {
  env <- attr(e1, "neuro_figure")
  # labs()/ggtitle() on a composed figure would land on an inner cell; turn
  # figure-level titles into a patchwork annotation styled like the family.
  if (inherits(e2, "ggplot2::labels") || inherits(e2, "labels")) {
    lab <- tryCatch(as.list(unclass(e2)), error = function(e) list())
    fl <- lab[intersect(names(lab), c("title", "subtitle", "caption"))]
    if (length(fl)) {
      tk <- if (is.environment(env) && !is.null(env$tokens)) env$tokens else neuro_style_tokens("light")
      e2 <- do.call(patchwork::plot_annotation, c(fl, list(theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, colour = tk$fg, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 10.5, colour = tk$muted, hjust = 0),
        plot.caption = ggplot2::element_text(size = 8.5, colour = tk$muted, hjust = 0),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.background = ggplot2::element_rect(fill = tk$card, colour = NA)))))
    }
  }
  if (is_plot_like(e2)) {
    # Composing with another plot: treat the whole figure as one element so
    # its internal grid does not collide with patchwork's layout.
    return(op(patchwork::wrap_elements(full = neuro_fig_resolve(e1)), e2))
  }
  res <- op(neuro_fig_strip(e1), e2)
  if (!is.environment(env)) return(res)
  build0 <- env$build
  neuro_figure(res, function(cv) op(build0(cv), e2), env$canvas, env$tokens)
}

#' Arithmetic on device-fitted figures
#'
#' \code{+} and \code{&} work as for any patchwork; the additions are replayed
#' when the figure is re-fitted to the device it is drawn on.
#' @param e1 A figure returned by a \code{plot_*} function.
#' @param e2 A ggplot2 or patchwork component.
#' @return A figure of the same class.
#' @keywords internal
#' @method + neuro_fig
#' @export
`+.neuro_fig` <- function(e1, e2) neuro_fig_add(e1, e2, `+`)

#' @rdname plus-.neuro_fig
#' @method & neuro_fig
#' @export
`&.neuro_fig` <- function(e1, e2) neuro_fig_add(e1, e2, `&`)

#' @rdname plus-.neuro_fig
#' @method | neuro_fig
#' @export
`|.neuro_fig` <- function(e1, e2) neuro_fig_add(e1, e2, `|`)

#' @rdname plus-.neuro_fig
#' @method / neuro_fig
#' @export
`/.neuro_fig` <- function(e1, e2) neuro_fig_add(e1, e2, `/`)

#' Is an object a plot (as opposed to a theme, scale or label component)?
#' @keywords internal
#' @noRd
is_plot_like <- function(x) {
  inherits(x, c("ggplot", "patchwork", "patch", "grob", "gtable"))
}

#' Print a device-fitted figure
#'
#' Re-fits the arrangement of the already-built tiles to the current device,
#' then prints it like any patchwork.
#' @param x A figure returned by a \code{plot_*} function.
#' @param ... Passed to the patchwork print method.
#' @return \code{x}, invisibly.
#' @keywords internal
#' @exportS3Method base::print
print.neuro_fig <- function(x, ...) {
  print(neuro_fig_resolve(x), ...)
  invisible(x)
}


#' @exportS3Method grid::grid.draw
grid.draw.neuro_fig <- function(x, recording = TRUE) {
  grid::grid.draw(neuro_fig_resolve(x), recording = recording)
}


#' Brain part of a bright-tissue slice mask
#'
#' On T1 images bright scalp fat passes the same intensity threshold as brain.
#' Scalp and skull form a rim of roughly constant thickness around the head,
#' so the brain is taken as the bright tissue inside the head mask eroded by
#' that rim (about 7\% of the head's width). Unlike region growing, this keeps
#' both hemispheres on high slices where the fissure separates them.
#' @param m Logical matrix (bright tissue).
#' @param head Logical matrix (head foreground).
#' @keywords internal
#' @noRd
brain_component <- function(m, head = m) {
  m[is.na(m)] <- FALSE
  head[is.na(head)] <- FALSE
  if (!any(m) || !any(head)) return(m)
  # Fill the head mask (air inside the sinuses or ventricles must not erode).
  idx <- which(head, arr.ind = TRUE)
  width <- max(diff(range(idx[, 1])), diff(range(idx[, 2])))
  r <- max(2L, as.integer(round(0.08 * width)))
  core <- !dilate_mask(!fill_rows_cols(head), r)
  out <- m & core
  # Morphological opening drops thin skull/marrow fragments left in the core.
  out <- !dilate_mask(!out, 2L)
  out <- dilate_mask(out, 2L) & m
  if (sum(out) < 0.2 * sum(m)) m else out
}

#' Fill a mask along rows and columns (a cheap hole fill for convex-ish heads)
#' @keywords internal
#' @noRd
fill_rows_cols <- function(m) {
  fill_line <- function(v) {
    w <- which(v)
    if (length(w)) v[min(w):max(w)] <- TRUE
    v
  }
  a <- t(apply(m, 1L, fill_line))
  b <- apply(m, 2L, fill_line)
  a & b
}

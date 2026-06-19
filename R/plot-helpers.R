# Internal utilities used by plotting helpers

#' @keywords internal
#' @noRd
utils::globalVariables(c("x", "y", "value", "z", "fill"))

#' Coerce a NeuroSlice (or matrix-like) to a numeric matrix
#' @keywords internal
#' @noRd
slice_to_matrix <- function(slc) {
  plain_matrix <- function(x) {
    dx <- dim(x)
    if (length(dx) != 2L) {
      stop("Expected a 2D object when coercing to a plain matrix.", call. = FALSE)
    }
    matrix(
      as.numeric(x),
      nrow = dx[1],
      ncol = dx[2],
      dimnames = dimnames(x)
    )
  }

  if (is.matrix(slc) && !isS4(slc)) {
    return(plain_matrix(slc))
  }
  m <- try(as.matrix(slc), silent = TRUE)
  if (!inherits(m, "try-error") && is.matrix(m)) return(plain_matrix(m))
  a <- try(as.array(slc), silent = TRUE)
  if (!inherits(a, "try-error")) {
    if (length(dim(a)) == 2) return(plain_matrix(a))
    if (length(dim(a)) == 3) return(plain_matrix(a[, , 1, drop = TRUE]))
  }
  # try common slots / accessors without committing to class internals
  data <- try(slc@data, silent = TRUE)
  if (!inherits(data, "try-error") && is.matrix(data)) return(plain_matrix(data))
  stop("Cannot coerce 'slc' to a numeric matrix. Provide a 2D matrix or a NeuroSlice with as.matrix().")
}

#' Convert a slice to a tidy data.frame ready for ggplot2::geom_raster
#' @param slc A NeuroSlice or 2D matrix
#' @param downsample Integer decimation factor (>=1)
#' @keywords internal
#' @noRd
slice_df <- function(slc, downsample = 1L) {
  m <- slice_to_matrix(slc)
  if (!is.numeric(m)) m <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m)))
  if (downsample > 1L) {
    idxr <- seq(1L, nrow(m), by = downsample)
    idxc <- seq(1L, ncol(m), by = downsample)
    m <- m[idxr, idxc, drop = FALSE]
  }
  # Build grid with x=cols, y=rows. y reversed later to match radiological display
  df <- expand.grid(x = seq_len(ncol(m)), y = seq_len(nrow(m)))
  df$value <- c(t(m))
  df
}

#' Orient slice-aligned matrices for raster rendering in world coordinates
#' @param slc A NeuroSlice describing the 2D slice geometry.
#' @param mat Numeric matrix aligned with \code{slc}.
#' @param alpha_map Optional numeric matrix aligned with \code{slc}.
#' @keywords internal
#' @noRd
orient_slice_for_raster <- function(slc, mat, alpha_map = NULL) {
  mat <- slice_to_matrix(mat)
  if (!is.null(alpha_map)) {
    alpha_map <- slice_to_matrix(alpha_map)
    if (!identical(dim(alpha_map), dim(mat))) {
      stop("'alpha_map' must have the same dimensions as 'mat'.", call. = FALSE)
    }
  }

  coords <- index_to_coord(space(slc), seq_len(length(slc)))
  xvals <- sort(unique(coords[, 1]))
  yvals <- sort(unique(coords[, 2]), decreasing = TRUE)

  out_mat <- matrix(NA_real_, nrow = length(yvals), ncol = length(xvals))
  out_alpha <- if (is.null(alpha_map)) NULL else {
    matrix(NA_real_, nrow = length(yvals), ncol = length(xvals))
  }

  col_idx <- match(coords[, 1], xvals)
  row_idx <- match(coords[, 2], yvals)
  fill_idx <- cbind(row_idx, col_idx)

  out_mat[fill_idx] <- as.numeric(mat)
  if (!is.null(out_alpha)) {
    out_alpha[fill_idx] <- as.numeric(alpha_map)
  }

  xr <- raster_extent_from_centers(xvals)
  yr <- raster_extent_from_centers(yvals)

  list(
    mat = out_mat,
    alpha_map = out_alpha,
    x = xvals,
    y = yvals,
    xmin = xr[1],
    xmax = xr[2],
    ymin = yr[1],
    ymax = yr[2]
  )
}

#' Compute raster annotation extents from pixel centers
#' @param centers Numeric vector of pixel centers.
#' @keywords internal
#' @noRd
raster_extent_from_centers <- function(centers) {
  vals <- sort(unique(centers))
  if (!length(vals)) return(c(0, 1))
  if (length(vals) == 1L) return(c(vals[1] - 0.5, vals[1] + 0.5))
  step <- stats::median(diff(vals))
  c(vals[1] - step / 2, vals[length(vals)] + step / 2)
}

#' Validate volumes that must already occupy the same image grid
#' @keywords internal
#' @noRd
assert_same_neuro_grid <- function(reference, ..., reference_name = "bgvol") {
  ref_dim <- dim(reference)
  if (length(ref_dim) != 3L) {
    stop(sprintf("`%s` must be a 3D volume.", reference_name), call. = FALSE)
  }

  ref_space <- space(reference)
  volumes <- list(...)
  volume_names <- names(volumes)
  if (is.null(volume_names)) {
    volume_names <- rep("", length(volumes))
  }

  for (i in seq_along(volumes)) {
    volume_name <- volume_names[[i]]
    if (!nzchar(volume_name)) {
      volume_name <- paste0("volume", i)
    }

    if (!identical(dim(volumes[[i]]), ref_dim)) {
      stop(
        sprintf("`%s` must have the same dimensions as `%s`.", volume_name, reference_name),
        call. = FALSE
      )
    }
    if (!identical(space(volumes[[i]]), ref_space)) {
      stop(
        sprintf("`%s` must be on the same NeuroSpace grid as `%s`.", volume_name, reference_name),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Validate shared slice-panel layout arguments
#' @keywords internal
#' @noRd
validate_slice_panel_args <- function(zlevels, along, dims, ncol) {
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along) || along < 1L || along > length(dims)) {
    stop("`along` must be one of 1, 2, or 3.", call. = FALSE)
  }

  zlevels <- as.integer(zlevels)
  if (!length(zlevels)) {
    stop("`zlevels` must contain at least one slice index.", call. = FALSE)
  }
  if (anyNA(zlevels) || any(zlevels < 1L | zlevels > dims[[along]])) {
    stop(
      sprintf("`zlevels` must be valid slice indices along axis %d.", along),
      call. = FALSE
    )
  }

  ncol <- as.integer(ncol)
  if (length(ncol) != 1L || is.na(ncol) || ncol < 1L) {
    stop("`ncol` must be a positive integer.", call. = FALSE)
  }

  list(zlevels = zlevels, along = along, ncol = ncol)
}

#' Draw a list of ggplot panels with optional layout-level labels
#' @keywords internal
#' @noRd
draw_plot_panel_grid <- function(plots, ncol, title = NULL, subtitle = NULL, caption = NULL,
                                 style = c("light", "dark")) {
  style <- match.arg(style)
  n <- length(plots)
  ncol <- min(ncol, n)
  nrow <- ceiling(n / ncol)

  has_title <- !is.null(title)
  has_subtitle <- !is.null(subtitle)
  has_caption <- !is.null(caption)
  top_rows <- as.integer(has_title) + as.integer(has_subtitle)
  bottom_rows <- as.integer(has_caption)

  heights <- grid::unit(rep(1, nrow + top_rows + bottom_rows), "null")
  if (top_rows > 0L || bottom_rows > 0L) {
    label_rows <- c(seq_len(top_rows), top_rows + nrow + seq_len(bottom_rows))
    heights[label_rows] <- grid::unit(1.2, "lines")
  }

  if (style == "dark") {
    bg <- "grey8"
    fg <- "grey94"
    muted <- "grey72"
  } else {
    bg <- "white"
    fg <- "grey10"
    muted <- "grey35"
  }

  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = bg, col = NA))
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow + top_rows + bottom_rows, ncol, heights = heights)
  ))
  on.exit(grid::upViewport(0), add = TRUE)

  label_row <- 0L
  if (has_title) {
    label_row <- label_row + 1L
    grid::grid.text(
      title,
      vp = grid::viewport(layout.pos.row = label_row, layout.pos.col = seq_len(ncol)),
      gp = grid::gpar(col = fg, fontface = "bold", fontsize = 12)
    )
  }
  if (has_subtitle) {
    label_row <- label_row + 1L
    grid::grid.text(
      subtitle,
      vp = grid::viewport(layout.pos.row = label_row, layout.pos.col = seq_len(ncol)),
      gp = grid::gpar(col = muted, fontsize = 10)
    )
  }

  for (i in seq_len(n)) {
    r <- ((i - 1L) %/% ncol) + 1L + top_rows
    c <- ((i - 1L) %% ncol) + 1L
    print(plots[[i]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
  }

  if (has_caption) {
    grid::grid.text(
      caption,
      vp = grid::viewport(layout.pos.row = top_rows + nrow + 1L, layout.pos.col = seq_len(ncol)),
      gp = grid::gpar(col = muted, fontsize = 9)
    )
  }

  invisible(plots)
}

#' Compute robust or data-based limits
#' @keywords internal
#' @noRd
compute_limits <- function(x, mode = c("robust","data"), probs = c(.02,.98)) {
  mode <- match.arg(mode)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  if (mode == "data") {
    rng <- range(x, finite = TRUE)
  } else {
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    rng <- c(min(q), max(q))
    if (!is.finite(rng[1]) || !is.finite(rng[2])) rng <- range(x, finite = TRUE)
    if (rng[1] == rng[2]) rng <- range(x, finite = TRUE)
  }
  rng
}

#' Shared color tokens for the plotting styles
#'
#' \code{"report"} renders dark brain tiles on a light "card" with bold/italic
#' typography (see the activation-map look); \code{"dark"}/\code{"light"} are the
#' classic single-tone styles.
#'
#' @keywords internal
#' @noRd
.plot_style_colors <- function(style = c("light", "dark", "report")) {
  style <- match.arg(style)
  switch(style,
    dark   = list(card = "grey8",   fg = "grey94", muted = "grey72", panel = "dark"),
    report = list(card = "#f6f6f4", fg = "grey10", muted = "grey40", panel = "dark"),
    light  = list(card = "white",   fg = "grey12", muted = "grey35", panel = "light")
  )
}

#' patchwork plot_annotation theme for the assembled figure (title/card)
#' @keywords internal
#' @noRd
annotation_theme <- function(style) {
  cols <- .plot_style_colors(style)
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = cols$card, colour = NA),
    plot.title    = ggplot2::element_text(face = "bold", colour = cols$fg, size = 16,
                                           margin = ggplot2::margin(b = 2)),
    plot.subtitle = ggplot2::element_text(face = "italic", colour = cols$muted, size = 11,
                                           margin = ggplot2::margin(b = 6)),
    plot.caption  = ggplot2::element_text(colour = cols$muted, hjust = 0),
    plot.margin   = grid::unit(c(12, 12, 10, 12), "pt")
  )
}

#' Assemble a panel block into a single figure (colorbar, optional legend, card)
#'
#' Shared by \code{plot_overlay()}, \code{plot_ortho()}, and \code{plot_montage()}
#' so the assembled look (right-hand colorbar, light/dark/report card, bold
#' title) is identical across the family.
#'
#' @param panel_block A patchwork/ggplot holding the slice panels.
#' @param lim Numeric length-2 limits for the colorbar.
#' @param cmap Palette name for the colorbar.
#' @param thresh Threshold to mark on the colorbar (0 to omit).
#' @param style "light"/"dark"/"report".
#' @param colorbar Logical; append the right-hand colorbar.
#' @param cbar_title Colorbar title.
#' @param legend Optional ggplot legend strip placed under the panels.
#' @param title,subtitle,caption Layout labels.
#' @return A patchwork object.
#' @keywords internal
#' @noRd
assemble_figure <- function(panel_block, lim, cmap, thresh = 0,
                            style = "light", colorbar = TRUE,
                            cbar_title = "value", legend = NULL,
                            title = NULL, subtitle = NULL, caption = NULL) {
  top <- panel_block
  if (isTRUE(colorbar)) {
    cbar <- make_overlay_colorbar(lim, cmap, thresh = thresh, style = style,
                                  title = cbar_title)
    top <- patchwork::wrap_plots(top, cbar, nrow = 1L, widths = c(1, 0.08))
  }
  combined <- top
  if (!is.null(legend)) {
    combined <- patchwork::wrap_plots(top, legend, ncol = 1L, heights = c(1, 0.11))
  }
  combined +
    patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption,
                               theme = annotation_theme(style))
}

#' Per-panel theme tweaks for borderless "report" brain tiles
#' @keywords internal
#' @noRd
report_tile_theme <- function() {
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    plot.margin  = grid::unit(c(5, 5, 5, 5), "pt")
  )
}

#' Faceted-montage theme for the "report" style (dark tiles on a light card)
#' @keywords internal
#' @noRd
report_facet_theme <- function() {
  sc <- .plot_style_colors("report")
  ggplot2::theme_minimal(base_size = 10) %+replace% ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = sc$card, colour = NA),
    panel.background = ggplot2::element_rect(fill = "black", colour = NA),
    panel.grid       = ggplot2::element_blank(),
    panel.spacing    = grid::unit(5, "pt"),
    axis.title = ggplot2::element_blank(),
    axis.text  = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "black", colour = NA),
    strip.text = ggplot2::element_text(colour = "grey92", face = "bold",
                                       margin = ggplot2::margin(t = 2, b = 3)),
    legend.position = "none",
    plot.margin = grid::unit(c(4, 4, 4, 4), "pt")
  )
}

#' Bounding-box crop window (world coords) around the brain + suprathreshold overlay
#'
#' Returns \code{list(xlim, ylim)} spanning, across all displayed slices, the
#' background voxels above \code{bg_thresh} unioned with overlay voxels above
#' \code{ov_thresh} (so clusters are never clipped), padded by \code{margin}.
#'
#' @keywords internal
#' @noRd
compute_crop_window <- function(bgvol, overlay, zlevels, along, bg_thresh,
                                ov_thresh = 0, margin = 0.06) {
  xs <- NULL; ys <- NULL
  accumulate <- function(o, keep) {
    if (!length(keep)) return(invisible())
    nx <- length(o$x)
    xi <- ((keep - 1L) %% nx) + 1L
    yi <- ((keep - 1L) %/% nx) + 1L
    xs <<- range(c(xs, o$x[xi]))
    ys <<- range(c(ys, o$y[yi]))
  }
  for (z in zlevels) {
    sl <- slice(bgvol, z, along = along)
    o  <- orient_slice_for_raster(sl, slice_to_matrix(sl))
    v  <- c(t(o$mat))
    accumulate(o, which(is.finite(v) & v > bg_thresh))

    slo <- slice(overlay, z, along = along)
    oo  <- orient_slice_for_raster(slo, slice_to_matrix(slo))
    vo  <- c(t(oo$mat))
    accumulate(oo, which(is.finite(vo) & abs(vo) > ov_thresh))
  }
  if (length(xs) < 2L || length(ys) < 2L || diff(xs) == 0 || diff(ys) == 0) {
    return(NULL)
  }
  padx <- diff(xs) * margin; pady <- diff(ys) * margin
  list(xlim = c(xs[1] - padx, xs[2] + padx),
       ylim = c(ys[1] - pady, ys[2] + pady))
}

#' Horizontal legend strip for an assembled overlay figure
#'
#' @param thresh Threshold value (0 to omit the threshold entry).
#' @param pos_col,neg_col Swatch colors for positive/negative activation.
#' @param symmetric Logical; show both pos/neg (signed) or a single entry.
#' @param style "light"/"dark"/"report".
#' @param plane Plane label ("Axial"/"Coronal"/"Sagittal").
#' @return A ggplot object (blank panel with swatches + labels).
#' @keywords internal
#' @noRd
make_overlay_legend <- function(thresh, pos_col, neg_col, symmetric = TRUE,
                                style = c("light", "dark", "report"),
                                plane = "Axial") {
  style <- match.arg(style)
  cols <- .plot_style_colors(style)
  fg <- cols$fg; muted <- cols$muted

  swatch <- function(p, x, fill, label, sub) {
    p +
      ggplot2::annotate("rect", xmin = x, xmax = x + 0.018, ymin = 0.42, ymax = 0.78,
                        fill = fill, colour = NA) +
      ggplot2::annotate("text", x = x + 0.028, y = 0.74, label = label, hjust = 0,
                        vjust = 1, colour = fg, fontface = "bold", size = 3.5) +
      ggplot2::annotate("text", x = x + 0.028, y = 0.40, label = sub, hjust = 0,
                        vjust = 1, colour = muted, fontface = "italic", size = 3.0)
  }

  p <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  if (isTRUE(symmetric)) {
    p <- swatch(p, 0.02, pos_col, "Positive activation", "higher than threshold")
    p <- swatch(p, 0.30, neg_col, "Negative activation", "lower than threshold")
  } else {
    p <- swatch(p, 0.02, pos_col, "Activation", "above threshold")
  }

  if (is.finite(thresh) && thresh > 0) {
    p <- p +
      ggplot2::annotate("segment", x = 0.58, xend = 0.625, y = 0.60, yend = 0.60,
                        colour = fg, linetype = "dashed", linewidth = 0.6) +
      ggplot2::annotate("text", x = 0.635, y = 0.60,
                        label = sprintf("Threshold (+/-%g)", thresh),
                        hjust = 0, vjust = 0.5, colour = fg, size = 3.3)
  }

  p <- p +
    ggplot2::annotate("text", x = 0.86, y = 0.74, label = plane, hjust = 0, vjust = 1,
                      colour = fg, fontface = "bold", size = 3.5) +
    ggplot2::annotate("text", x = 0.86, y = 0.40, label = "orientation", hjust = 0,
                      vjust = 1, colour = muted, fontface = "italic", size = 3.0)

  p +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = cols$card, colour = NA),
      plot.margin = grid::unit(c(4, 10, 4, 10), "pt")
    )
}

#' Self-tuning nonlinear alpha mapping for statistical overlays
#'
#' Parameters for the \code{"soft"} alpha mode, where per-voxel opacity is
#' \deqn{alpha(m) = clamp((m - lo) / (hi - lo), 0, 1)^{gamma}.}
#' The knee \code{lo} (below which alpha is 0) defaults to the threshold, or --
#' when no threshold is set -- to the median in-mask magnitude, a robust
#' noise-floor proxy, so near-zero values stay transparent. \code{gamma} is
#' tuned so the median \emph{displayed} (supra-knee) magnitude maps to a faint
#' \code{alpha_mid}; this pushes the noisy bulk toward transparency while the
#' upper tail saturates to opaque, and adapts automatically to the value
#' distribution. \code{gamma} is clamped to \code{[gamma_min, gamma_max]} so the
#' curve stays convex (low values are never boosted).
#'
#' @param mags Numeric vector of overlay magnitudes (typically \code{abs(values)}).
#' @param thresh Hard/soft threshold; used as the knee when \code{> 0}.
#' @param cap Upper anchor (mapped to alpha 1); defaults to \code{max(mags)}.
#' @param gamma Optional fixed exponent; \code{NULL} auto-tunes it.
#' @param alpha_mid Target alpha for the median displayed magnitude.
#' @param gamma_min,gamma_max Clamp range for the tuned exponent (>= 1 keeps the
#'   curve convex; the lower bound guarantees a visible nonlinearity).
#' @return A list with \code{lo}, \code{hi}, and \code{gamma}.
#' @keywords internal
#' @noRd
soft_alpha_params <- function(mags, thresh = 0, cap = NULL, gamma = NULL,
                              alpha_mid = 0.2, gamma_min = 1.5, gamma_max = 5) {
  mags <- mags[is.finite(mags) & mags > 0]
  knee <- if (isTRUE(thresh > 0)) {
    thresh
  } else if (length(mags)) {
    stats::median(mags)
  } else {
    0
  }
  hi <- if (!is.null(cap) && is.finite(cap)) cap else if (length(mags)) max(mags) else knee + 1
  if (!is.finite(hi) || hi <= knee) hi <- knee + 1

  if (is.null(gamma)) {
    supra <- mags[mags > knee]
    if (length(supra) >= 10L) {
      t_med <- stats::median((supra - knee) / (hi - knee))
      t_med <- min(max(t_med, 1e-3), 0.999)
      gamma <- log(alpha_mid) / log(t_med)
    } else {
      gamma <- 2
    }
    gamma <- min(max(gamma, gamma_min), gamma_max)
  }
  list(lo = knee, hi = hi, gamma = gamma)
}

#' Resolve a display-range argument to numeric limits
#'
#' Accepts either a mode string (\code{"robust"}/\code{"data"}, resolved against
#' \code{values} via [compute_limits()]) or an explicit numeric \code{c(lo, hi)}.
#' This lets plotting functions pin a fixed scale (e.g. \code{ov_range = c(-6, 6)})
#' for cross-panel / cross-subject comparability.
#'
#' @keywords internal
#' @noRd
resolve_display_limits <- function(range_arg, values, probs = c(.02, .98)) {
  if (is.numeric(range_arg)) {
    if (length(range_arg) != 2L || any(!is.finite(range_arg)) ||
        range_arg[1] == range_arg[2]) {
      stop("A numeric range must be two distinct finite values, c(lo, hi).",
           call. = FALSE)
    }
    return(range(range_arg))
  }
  mode <- match.arg(range_arg, c("robust", "data"))
  compute_limits(values, mode = mode, probs = probs)
}

#' Build a standalone colorbar plot for an overlay statistic
#'
#' Renders a thin value->color strip with a right-hand axis, and marks the
#' threshold (and its mirror, for signed maps). Used to give the composited
#' overlay an explicit scale, since it is drawn as a raster annotation rather
#' than a mapped ggplot layer.
#'
#' @param limits Numeric length-2 display limits.
#' @param cmap Overlay palette name or color vector.
#' @param thresh Threshold to mark (0 for none).
#' @param style "light" or "dark".
#' @param title Colorbar title.
#' @return A ggplot object.
#' @keywords internal
#' @noRd
make_overlay_colorbar <- function(limits, cmap, thresh = 0,
                                  style = c("light", "dark", "report"),
                                  title = "value") {
  style <- match.arg(style)
  sc <- .plot_style_colors(style)
  cols <- resolve_cmap(cmap, 256)
  yy <- seq(limits[1], limits[2], length.out = 256)
  df <- data.frame(x = 1, y = yy, fill = yy)
  fg <- sc$fg
  bg <- sc$card

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = fill)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = cols, limits = limits, guide = "none") +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::scale_x_continuous(breaks = NULL,
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),
      axis.text.y = ggplot2::element_text(colour = fg),
      axis.ticks.y = ggplot2::element_line(colour = fg),
      axis.text.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = fg, size = 9, hjust = 0.5)
    )

  if (is.finite(thresh) && thresh > 0) {
    marks <- c(thresh, -thresh)
    marks <- marks[marks >= limits[1] & marks <= limits[2]]
    if (length(marks)) {
      p <- p + ggplot2::geom_hline(yintercept = marks, colour = fg,
                                   linetype = "dashed", linewidth = 0.3)
    }
  }
  p
}

#' Apply the `enhance` argument of a plotting function to a volume
#'
#' Coerces the user-facing `enhance` argument into a call to
#' \code{\link{enhance_stat_map}}. Accepts \code{FALSE}/\code{NULL} (no-op),
#' \code{TRUE} (default enhancement), or a named \code{list} of arguments
#' forwarded to \code{enhance_stat_map()}.
#'
#' @keywords internal
#' @noRd
apply_enhance_arg <- function(vol, enhance, mask = NULL) {
  if (is.null(enhance) || isFALSE(enhance)) return(vol)
  args <- list(vol = vol)
  if (!is.null(mask)) args$mask <- mask
  if (is.list(enhance)) {
    args <- utils::modifyList(args, enhance)
  } else if (!isTRUE(enhance)) {
    stop("`enhance` must be TRUE, FALSE, or a named list of enhance_stat_map() arguments.",
         call. = FALSE)
  }
  do.call(enhance_stat_map, args)
}

#' Simple rescale
#' @keywords internal
#' @noRd
rescale01 <- function(x, from) {
  if (is.null(from) || any(!is.finite(from))) return(x)
  a <- from[1]; b <- from[2]
  if (b == a) return(rep(.5, length(x)))
  (x - a) / (b - a)
}

#' Map numeric matrix to RGBA colors
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 vector of display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional numeric matrix (same dims as \code{mat}, values in
#'   \code{[0,1]}) providing per-pixel alpha. Effective alpha per pixel is
#'   \code{alpha_map[i] * alpha}. When \code{NULL} the scalar \code{alpha} is
#'   applied to all pixels.
#' @keywords internal
#' @noRd
matrix_to_colors <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                              alpha_map = NULL) {
  cols <- resolve_cmap(cmap, 256)
  if (is.null(limits)) limits <- range(mat[is.finite(mat)])
  s <- rescale01(as.numeric(mat), limits)
  s[!is.finite(s)] <- NA_real_
  # index into palette
  idx <- 1 + floor(s * (length(cols) - 1))
  col <- cols[pmax(1, pmin(length(cols), idx))]
  # apply alpha
  if (is.null(alpha_map)) {
    grDevices::adjustcolor(col, alpha.f = alpha)
  } else {
    eff <- as.numeric(alpha_map) * alpha
    mapply(function(c, a) grDevices::adjustcolor(c, alpha.f = a),
           col, eff, USE.NAMES = FALSE)
  }
}

#' Build a numeric RGBA array for a matrix
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional per-pixel alpha matrix (same dims as \code{mat},
#'   values in \code{[0,1]}). Effective alpha is \code{alpha_map * alpha}.
#' @keywords internal
#' @noRd
matrix_to_rgba <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                           alpha_map = NULL) {
  cols <- resolve_cmap(cmap, 256)
  if (is.null(limits)) limits <- range(mat[is.finite(mat)])

  s <- rescale01(as.numeric(mat), limits)
  idx <- 1 + floor(s * (length(cols) - 1))
  idx[!is.finite(idx)] <- 1L
  idx <- pmax(1L, pmin(length(cols), idx))

  rgb <- grDevices::col2rgb(cols[idx], alpha = TRUE) / 255

  if (is.null(alpha_map)) {
    eff_alpha <- rep(alpha, length(s))
  } else {
    eff_alpha <- as.numeric(alpha_map) * alpha
  }
  eff_alpha[!is.finite(eff_alpha)] <- 0
  eff_alpha <- pmax(0, pmin(1, eff_alpha))
  eff_alpha[!is.finite(s)] <- 0

  rgba <- array(0, dim = c(nrow(mat), ncol(mat), 4L))
  rgba[, , 1] <- matrix(rgb[1, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 2] <- matrix(rgb[2, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 3] <- matrix(rgb[3, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 4] <- matrix(eff_alpha, nrow = nrow(mat), ncol = ncol(mat))

  rgba
}

#' Create a rasterGrob from a numeric matrix using a palette
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional per-pixel alpha matrix (same dims as \code{mat},
#'   values in \code{[0,1]}). Effective alpha is \code{alpha_map * alpha}.
#' @keywords internal
#' @noRd
matrix_to_raster_grob <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                                   alpha_map = NULL) {
  rgba <- matrix_to_rgba(
    mat = mat,
    cmap = cmap,
    limits = limits,
    alpha = alpha,
    alpha_map = alpha_map
  )
  grid::rasterGrob(rgba, interpolate = FALSE)
}

#' Coordinate helper: fixed aspect and reversed y for radiological convention
#' @keywords internal
#' @noRd
coord_neuro_fixed <- function() {
  list(ggplot2::coord_fixed(), ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = 0)))
}

#' Add L/R and A/P/S/I annotations (optional)
#'
#' @param plane "axial", "coronal", or "sagittal"
#' @param dims c(nrow, ncol) of the slice matrix
#' @param gp grid::gpar style
#' @return A ggplot2 layer with annotation_custom grobs
#' @export
annotate_orientation <- function(plane = c("axial","coronal","sagittal"),
                                 dims, gp = grid::gpar(col = "white", cex = .9, fontface = "bold")) {
  plane <- match.arg(plane)
  nr <- dims[1]; nc <- dims[2]
  lab_lr <- list(left  = "L", right = "R")
  lab_tb <- switch(plane,
                   axial   = list(top = "A", bottom = "P"),
                   coronal = list(top = "S", bottom = "I"),
                   sagittal= list(top = "S", bottom = "I"))
  layers <- list(
    ggplot2::annotation_custom(grid::textGrob(lab_lr$left,  gp = gp),
                               xmin = 0.5, xmax = 0.5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(lab_lr$right, gp = gp),
                               xmin = nc + .5, xmax = nc + .5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(lab_tb$top,   gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = .5, ymax = .5),
    ggplot2::annotation_custom(grid::textGrob(lab_tb$bottom,gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = nr + .5, ymax = nr + .5)
  )
  layers
}

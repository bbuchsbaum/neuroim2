#' Montage of slices through a volume
#'
#' Draws a grid of slices through one volume (typically a structural image) in
#' the same style as \code{\link{plot_overlay}}: black tiles cropped to the
#' head, world-coordinate labels, L/R markers on the first tile, and a grid
#' that fills the canvas. Also accepts a list of \code{NeuroSlice} objects or
#' plain matrices.
#'
#' @param x A 3D volume, or a list of \code{NeuroSlice} objects / matrices.
#' @param zlevels Slices to plot when \code{x} is a volume: indices along
#'   \code{along} (\code{unit = "index"}, the default) or world coordinates
#'   (\code{unit = "mm"}). \code{NULL} (default) picks \code{n_slices} slices
#'   spread over the brain.
#' @param along Native voxel-grid axis along which to slice. For canonically
#'   ordered images, 1 = sagittal, 2 = coronal, and 3 = axial. Display
#'   orientation is inferred from the image affine.
#' @param cmap Palette name or vector (see [resolve_cmap()]).
#' @param range "robust" (quantile-based), "data" (min/max), or an explicit
#'   numeric \code{c(lo, hi)}.
#' @param probs Quantiles for `range="robust"`.
#' @param ncol Number of columns in the facet layout.
#' @param downsample Integer decimation for speed.
#' @param title,subtitle,caption Optional ggplot labels.
#' @param style Visual style: \code{"light"}, \code{"dark"}, or \code{"report"}
#'   (light card, dark cropped tiles, typography, and a colorbar -- matching
#'   \code{\link{plot_overlay}}'s report look).
#' @param cbar_title Character; the quantity label drawn above the colorbar.
#'   Supplying it explicitly also turns the colorbar on.
#' @param crop,interpolate Logical; crop to the head bounding box (volume
#'   input) / smooth the raster. Both default to \code{TRUE}.
#' @param colorbar Logical or \code{NULL}. \code{NULL} (default) shows a slim
#'   colorbar only for non-grayscale palettes or when \code{cbar_title} is
#'   supplied; arbitrary structural intensity units carry no information.
#' @param unit \code{"index"} (default) or \code{"mm"}: how \code{zlevels} is
#'   interpreted for volume input. Panels are always labelled in world
#'   coordinates (mm).
#' @param annotate Logical; draw L/R (or A/P) orientation letters on the first
#'   panel.
#' @param n_slices Number of slices chosen automatically when \code{zlevels}
#'   is \code{NULL}. Automatic slices are spread over the extent of the brain
#'   (bright tissue), skipping neck, scalp-only and empty planes.
#' @param draw Logical; if \code{TRUE}, also print the figure immediately (and
#'   return it invisibly). By default it is returned visibly, like any ggplot.
#' @param canvas Optional \code{c(width, height)} in inches to fit the layout
#'   to (default: the open device, else 10 x 7.5 in).
#' @return A figure (class \code{neuro_fig}, a \pkg{patchwork} wrapping one
#'   faceted ggplot with one facet per slice), returned visibly; invisibly when
#'   \code{draw = TRUE}. Its layout is re-fitted to the device it is drawn on,
#'   so \code{ggsave()} at any size gives a filled, centred grid.
#' @examples
#' \donttest{
#' bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
#' p <- plot_montage(bg, title = "MNI152 (downsampled)")
#' ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
#' }
#' @family plot_neuro
#' @export
plot_montage <- function(
  x, zlevels = NULL, along = 3L,
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  ncol = NULL, downsample = 1L,
  title = NULL, subtitle = NULL, caption = NULL,
  style = c("light", "dark", "report"),
  crop = TRUE, interpolate = TRUE, cbar_title = "value", colorbar = NULL,
  unit = c("index", "mm"), annotate = TRUE, n_slices = 12L,
  draw = FALSE, canvas = NULL
) {
  unit_missing <- missing(unit)
  cbar_title_given <- !missing(cbar_title)
  cbar_title <- validate_cbar_title(cbar_title)
  style <- match_choice(style, c("light", "dark", "report"))
  unit <- match_choice(unit, c("index", "mm"), "unit")
  check_count(n_slices, "n_slices")
  if (inherits(x, "NeuroVol")) unit_hint(zlevels, unit_missing, x)
  tokens <- neuro_style_tokens(style)
  is_report <- identical(style, "report")
  do_crop   <- if (is.null(crop)) TRUE else isTRUE(crop)
  interp_bg <- isTRUE(interpolate)

  neuro_volume <- inherits(x, "NeuroVol")
  neuro_slice_list <- is.list(x) && length(x) &&
    all(vapply(x, inherits, logical(1), what = "NeuroSlice"))
  labels <- NULL
  first_oriented <- NULL
  if (neuro_volume) {
    zlevels <- resolve_slice_levels(zlevels, x, along, unit = unit, n = n_slices)
    zlevels <- validate_slice_panel_args(zlevels, along, dim(x), 1L)$zlevels
    # Anatomically orient the regular native grid. Raw oblique world
    # coordinates are not a valid geom_raster grid and would be shifted.
    dfl <- lapply(zlevels, function(z) {
      oriented <- orient_volume_slice_for_raster(
        x, z, along = along, downsample = downsample
      )
      if (is.null(first_oriented)) first_oriented <<- oriented
      df <- oriented_raster_df(oriented)
      df$z <- z
      df
    })
    labels <- vapply(zlevels, function(z) slice_world_label(x, z, along), character(1))
  } else if (neuro_slice_list) {
    dfl <- Map(function(sl, idx) {
      oriented <- orient_slice_for_raster(
        sl, slice_to_matrix(sl), downsample = downsample
      )
      if (is.null(first_oriented)) first_oriented <<- oriented
      df <- oriented_raster_df(oriented)
      df$z <- idx
      df
    }, x, seq_along(x))
  } else {
    # Fallback for plain matrices and other pixel-grid slice providers.
    make_slice <- function(z) {
      sl <- slice(x, z, along = along)
      df <- slice_df(sl, downsample = downsample)
      df$z <- z
      df
    }
    if (is.list(x) && is.null(zlevels)) {
      dfl <- Map(function(sl, idx) { df <- slice_df(sl, downsample); df$z <- idx; df },
                 x, seq_along(x))
    } else {
      if (is.null(zlevels)) {
        zlevels <- unique(round(seq(1, dim(x)[along], length.out = n_slices)))
      }
      dfl <- lapply(zlevels, make_slice)
    }
  }
  df <- do.call(rbind, dfl)
  facets <- unique(df$z)
  if (is.null(labels)) labels <- paste("slice", facets)
  lim <- if (neuro_volume || neuro_slice_list) {
    background_display_limits(range, df$value, probs = probs)
  } else {
    resolve_display_limits(range, df$value, probs = probs)
  }

  # Crop to the head bounding box (volume path only).
  window <- NULL
  if (isTRUE(do_crop) && neuro_volume) {
    window <- foreground_crop_window(x, zlevels, along)
  }
  if (is.null(window)) {
    window <- list(xlim = raster_extent_from_centers(df$x),
                   ylim = raster_extent_from_centers(df$y))
  }
  oriented_input <- neuro_volume || neuro_slice_list

  show_cbar <- if (is.null(colorbar)) {
    length(cmap) > 1L || !(tolower(cmap[1]) %in% c("grays", "gray", "grey", "greys")) ||
      cbar_title_given
  } else isTRUE(colorbar)

  user_ncol <- ncol
  xw <- diff(window$xlim); yw <- diff(window$ylim)
  pad <- 0.035 * max(xw, yw)
  off <- 0.0045 * max(xw, yw)
  lab_df <- data.frame(z = facets, label = labels,
                       x = window$xlim[1] + pad, y = window$ylim[2] - pad)
  cbar <- if (show_cbar) {
    neuro_colorbar(lim, resolve_cmap(cmap, 256), positions = FALSE,
                   tokens = tokens, title = cbar_title)
  }

  # Everything that depends on the canvas: the number of columns, the margins
  # that centre the tile block (keeping the title aligned with its left
  # edge), and the inset colorbar position.
  build <- function(ds) {
    head_in <- title_block_height(title, subtitle)
    cap_in <- if (is.null(caption)) 0 else 0.28
    cb_in <- if (show_cbar) 0.8 else 0
    avail_w <- max(ds[1] - 0.5 - cb_in, 1)
    avail_h <- max(ds[2] - 0.5 - head_in - cap_in, 1)
    tile_aspect <- xw / yw
    nc <- if (is.null(user_ncol)) auto_ncol(length(facets), tile_aspect, avail_w / avail_h) else user_ncol
    nc <- max(1L, min(as.integer(nc), length(facets)))
    nr <- ceiling(length(facets) / nc)
    grid_aspect <- nc * tile_aspect / nr
    block_w <- min(avail_w, avail_h * grid_aspect)
    block_h <- block_w / grid_aspect
    spare_w <- max(avail_w - block_w, 0)
    spare_h <- max(avail_h - block_h, 0)
    margins <- 72 * c(0.25 + spare_h / 2, 0.25 + spare_w / 2 + cb_in,
                      0.25 + spare_h / 2, 0.25 + spare_w / 2)

    p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = value), interpolate = interp_bg) +
      ggplot2::facet_wrap(~ z, ncol = nc, scales = "fixed") +
      scale_fill_neuro(cmap = cmap, limits = lim, na.value = tokens$tile, guide = "none")

    if (oriented_input) {
      p <- p +
        ggplot2::geom_text(data = halo_offsets(lab_df, off),
                           ggplot2::aes(x = x, y = y, label = label),
                           inherit.aes = FALSE, hjust = 0, vjust = 1,
                           colour = tokens$tile, alpha = 0.85, size = 3.1) +
        ggplot2::geom_text(data = lab_df, ggplot2::aes(x = x, y = y, label = label),
                           inherit.aes = FALSE, hjust = 0, vjust = 1,
                           colour = tokens$tile_fg, size = 3.1)
      if (isTRUE(annotate) && !is.null(first_oriented)) {
        ol <- first_oriented$labels
        or_df <- data.frame(
          z = facets[[1L]],
          x = c(window$xlim[1] + pad * 0.8, window$xlim[2] - pad * 0.8),
          y = mean(window$ylim),
          label = c(ol[["left"]], ol[["right"]]),
          hjust = c(0, 1)
        )
        p <- p +
          ggplot2::geom_text(data = halo_offsets(or_df, off),
                             ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
                             inherit.aes = FALSE, vjust = 0.5, fontface = "bold",
                             colour = tokens$tile, alpha = 0.85, size = 3.4) +
          ggplot2::geom_text(data = or_df,
                             ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
                             inherit.aes = FALSE, vjust = 0.5, fontface = "bold",
                             colour = tokens$tile_fg, size = 3.4)
      }
      p <- p + ggplot2::coord_fixed(xlim = window$xlim, ylim = window$ylim, expand = FALSE)
    } else {
      p <- p + coord_neuro_fixed()
    }

    p <- p +
      montage_theme(tokens, labelled = oriented_input, has_subtitle = !is.null(subtitle)) +
      ggplot2::theme(plot.margin = grid::unit(margins, "pt")) +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
    # The family's shared colorbar, inset just right of the tile block
    # (positions in inches on the canvas, converted to npc).
    if (show_cbar) {
      x_r <- 0.25 + spare_w / 2 + block_w
      y_b <- 0.25 + spare_h / 2 + cap_in
      bar_h <- min(block_h * 0.5, 3)
      p <- p + patchwork::inset_element(
        cbar,
        left = (x_r + 0.05) / ds[1], right = (x_r + cb_in) / ds[1],
        bottom = (y_b + (block_h - bar_h) / 2) / ds[2],
        top = (y_b + (block_h + bar_h) / 2 + 0.2) / ds[2],
        align_to = "full", clip = FALSE
      )
    }
    # Always a patchwork, so the card colour fills the whole canvas (a
    # fixed-aspect ggplot only paints behind its own plot table) and the
    # return type does not depend on the options.
    p + patchwork::plot_annotation(theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = tokens$card, colour = NA)))
  }

  cv0 <- canvas_size(canvas)
  p <- build(cv0)
  if (is.null(canvas)) p <- neuro_figure(p, build, cv0, tokens)
  if (isTRUE(draw)) {
    print(p)
    return(invisible(p))
  }
  p
}

#' Theme for faceted montages
#' @keywords internal
#' @noRd
montage_theme <- function(tokens, labelled = TRUE, has_subtitle = FALSE) {
  ggplot2::theme_void(base_size = 10) %+replace% ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = tokens$card, colour = NA),
    panel.background = ggplot2::element_rect(fill = tokens$tile, colour = NA),
    panel.spacing    = grid::unit(tokens$gutter, "pt"),
    strip.text       = if (labelled) ggplot2::element_blank() else
      ggplot2::element_text(colour = tokens$fg, size = 9, margin = ggplot2::margin(b = 3)),
    legend.title = ggplot2::element_text(colour = tokens$fg, size = 9,
                                         margin = ggplot2::margin(b = 5)),
    legend.text  = ggplot2::element_text(colour = tokens$fg, size = 8.5),
    legend.background = ggplot2::element_rect(fill = tokens$card, colour = NA),
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold", colour = tokens$fg, size = 14,
                                       hjust = 0, margin = ggplot2::margin(b = if (has_subtitle) 4 else 11)),
    plot.subtitle = ggplot2::element_text(colour = tokens$muted, size = 10.5, hjust = 0,
                                          margin = ggplot2::margin(b = 11)),
    plot.caption = ggplot2::element_text(colour = tokens$muted, size = 8.5, hjust = 0,
                                         margin = ggplot2::margin(t = 6)),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = grid::unit(c(10, 10, 8, 10), "pt")
  )
}

#' Eight offset copies of a label data frame (for a text halo)
#' @keywords internal
#' @noRd
halo_offsets <- function(df, off) {
  dx <- off * c(-1, 1, 0, 0, -0.7, 0.7, -0.7, 0.7)
  dy <- off * c(0, 0, -1, 1, -0.7, -0.7, 0.7, 0.7)
  out <- df[rep(seq_len(nrow(df)), each = 8L), , drop = FALSE]
  out$x <- out$x + rep(dx, nrow(df))
  out$y <- out$y + rep(dy, nrow(df))
  out
}

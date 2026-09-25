#' Orthogonal three-plane view with optional crosshairs and overlay
#'
#' Draws sagittal, coronal and axial sections through one point, optionally
#' with a thresholded statistical map on top (the classic "stat map" view).
#' All three views share one physical scale and one head-bounding-box crop, so
#' the crosshair lines up across views.
#'
#' @param vol A 3D background volume.
#' @param coord Length-3 coordinate of the target point: voxel indices
#'   (\code{unit = "index"}, default) or world coordinates in mm
#'   (\code{unit = "mm"}). \code{NULL} (default) uses the peak absolute value
#'   of \code{overlay} when one is given, otherwise the volume centre.
#' @param unit \code{"index"} or \code{"mm"}: how \code{coord} is interpreted.
#' @param cmap Palette for the background.
#' @param range Background intensity limits shared by all panels:
#'   \code{"robust"} (computed over head voxels), \code{"data"}, or numeric
#'   \code{c(lo, hi)}.
#' @param probs Quantiles for robust scaling.
#' @param crosshair Logical; draw the (gapped) crosshair.
#' @param annotate Logical; draw orientation letters on every view.
#' @param downsample Integer decimation for speed.
#' @param title,subtitle,caption Optional figure labels.
#' @param draw Logical; if \code{TRUE}, also print the figure immediately (and
#'   return it invisibly). By default the figure is returned visibly.
#' @param style Visual style: \code{"light"}, \code{"dark"}, or
#'   \code{"report"}.
#' @param enhance Display-only enhancement of an unsmoothed statistical
#'   \code{vol}; see \code{\link{plot_overlay}}.
#' @param cbar_title Character; the quantity label drawn above the colorbar.
#'   Supplying it explicitly also turns the colorbar on.
#' @param crop Logical; crop views to the head bounding box.
#' @param interpolate Logical; smooth the background raster (default
#'   \code{TRUE}).
#' @param colorbar Logical or \code{NULL}. \code{NULL} (default) shows a
#'   colorbar when it carries information: an \code{overlay} is given, the
#'   background uses a non-grayscale palette, or \code{cbar_title} is supplied.
#' @param assemble Logical; if \code{TRUE} (default) return one assembled
#'   \pkg{patchwork} figure; if \code{FALSE} return the named list of the
#'   \code{axial}, \code{coronal} and \code{sagittal} ggplots.
#' @param overlay Optional 3D statistical volume on the same grid as \code{vol},
#'   drawn over all three views.
#' @param ov_thresh,ov_cmap,ov_range,ov_alpha,ov_alpha_mode,ov_symmetric,ov_cap
#'   Overlay threshold, palette, scaling, opacity and opacity mode; identical
#'   in meaning to the same arguments of \code{\link{plot_overlay}}.
#' @param canvas Optional \code{c(width, height)} in inches to fit the layout
#'   to (default: the open device).
#' @return A figure (class \code{neuro_fig}, a \pkg{patchwork} whose
#'   layout is re-fitted to the device it is drawn on) when
#'   \code{assemble = TRUE} or a named list
#'   of ggplots (\code{assemble = FALSE}); invisibly when \code{draw = TRUE}.
#' @details The affine determines which native voxel axis is nearest each
#'   anatomical plane and how that plane must be permuted or flipped for
#'   display. Oblique images are shown on their regular native voxel planes;
#'   values are not silently resampled. Use \code{deoblique()} or
#'   \code{resample_to()} first when true cardinal-plane sections are required.
#'
#'   Each view is labelled with the world coordinate of its plane (mm).
#' @examples
#' \donttest{
#' bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
#' p <- plot_ortho(bg, coord = c(24, 26, 26))
#' ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 9, height = 3.5)
#' }
#' @family plot_neuro
#' @export
plot_ortho <- function(
  vol, coord = NULL, unit = c("index","mm"),
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  crosshair = TRUE, annotate = TRUE, downsample = 1L,
  title = NULL, subtitle = NULL, caption = NULL,
  draw = FALSE, style = c("light", "dark", "report"), enhance = FALSE,
  crop = TRUE, interpolate = TRUE, cbar_title = "value",
  colorbar = NULL, assemble = TRUE,
  overlay = NULL, ov_thresh = 0, ov_cmap = NULL,
  ov_range = c("robust", "data"), ov_alpha = 1,
  ov_alpha_mode = c("binary", "proportional", "ramp", "soft"),
  ov_symmetric = NULL, ov_cap = NULL, canvas = NULL
) {
  unit_missing <- missing(unit)
  cbar_title_given <- !missing(cbar_title)
  cbar_title <- validate_cbar_title(cbar_title)
  unit <- match_choice(unit, c("index", "mm"), "unit")
  style <- match_choice(style, c("light", "dark", "report"))
  ov_alpha_mode <- match_choice(ov_alpha_mode, c("binary", "proportional", "ramp", "soft"),
                                "ov_alpha_mode")
  tokens <- neuro_style_tokens(style)
  is_report <- identical(style, "report")
  do_crop   <- if (is.null(crop)) TRUE else isTRUE(crop)
  interp_bg <- isTRUE(interpolate)
  if (!is.null(overlay)) assert_same_neuro_grid(vol, overlay = overlay, reference_name = "vol")
  ov_thresh <- check_overlay_args(ov_thresh, ov_alpha)
  unit_hint(coord, unit_missing, vol, what = "coord")

  # Optional display-only enhancement of an unsmoothed statistical volume.
  vol <- apply_enhance_arg(vol, enhance)
  d <- dim(vol)[1:3]

  if (is.null(coord)) {
    coord <- if (!is.null(overlay)) {
      a <- abs(as.array(overlay))
      a[!is.finite(a)] <- 0
      as.integer(arrayInd(which.max(a), dim(a)))
    } else {
      round(d / 2)
    }
  } else if (unit == "mm") {
    if (length(coord) != 3L || any(!is.finite(coord))) {
      cli::cli_abort("{.arg coord} must be a length-3 world coordinate (mm).", call = NULL)
    }
    mm <- coord
    coord <- as.integer(round(coord_to_grid(space(vol), matrix(coord, ncol = 3))))
    if (any(coord < 1L | coord > d)) {
      rng <- world_bounds(vol)
      cli::cli_abort(c(
        "{.arg coord} = ({paste(format(mm), collapse = ', ')}) mm lies outside the image.",
        "i" = "The image spans x {rng[1,1]} to {rng[2,1]}, y {rng[1,2]} to {rng[2,2]}, z {rng[1,3]} to {rng[2,3]} mm."
      ), call = NULL)
    }
  }
  if (length(coord) != 3L || anyNA(coord) || any(!is.finite(coord))) {
    cli::cli_abort("{.arg coord} must be a length-3 coordinate.", call = NULL)
  }
  coord <- as.integer(round(coord))
  if (any(coord < 1L | coord > d)) {
    cli::cli_abort(c(
      "{.arg coord} must be a valid voxel coordinate.",
      "i" = "Valid indices are 1 to {d[1]}, 1 to {d[2]} and 1 to {d[3]}.",
      "i" = "To give the point in world coordinates, use {.code unit = \"mm\"}."
    ), call = NULL)
  }

  # Find the native grid axis nearest each anatomical normal. This preserves
  # axial/coronal/sagittal semantics when a valid NIfTI affine permutes the
  # voxel axes, and it keeps oblique native planes on regular raster grids.
  directions <- perm_mat(axes(space(vol)))
  native_axis_for_world <- vapply(
    seq_len(3L),
    function(world_axis) which.max(abs(directions[world_axis, ])),
    integer(1)
  )
  if (anyDuplicated(native_axis_for_world)) {
    cli::cli_abort("Volume axes do not define three distinct anatomical directions.",
                   call = NULL)
  }

  planes <- list(axial = 3L, coronal = 2L, sagittal = 1L)
  info <- lapply(planes, function(world_normal) {
    along <- native_axis_for_world[[world_normal]]
    oriented <- orient_volume_slice_for_raster(vol, z = coord[[along]], along = along,
                                               downsample = downsample)
    list(along = along, z = coord[[along]], oriented = oriented,
         cross = slice_grid_to_display(oriented, coord[-along]))
  })

  # Shared intensity limits across the three views.
  all_vals <- unlist(lapply(info, function(i) as.numeric(i$oriented$mat)))
  lim <- background_display_limits(range, all_vals, probs = probs)

  # Overlay scale: the same helper plot_overlay() uses.
  scale <- NULL
  if (!is.null(overlay)) {
    ov_all <- as.numeric(as.array(overlay))
    scale <- overlay_scale(ov_all, ov_range = ov_range, probs = probs,
                           thresh = ov_thresh, ov_cmap = ov_cmap,
                           ov_symmetric = ov_symmetric, ov_cap = ov_cap)
    cap <- max(abs(scale$lim))
    soft <- if (ov_alpha_mode == "soft") soft_alpha_params(abs(ov_all), thresh = ov_thresh, cap = cap)
    alpha_fun <- overlay_alpha_fun(ov_alpha_mode, ov_thresh, cap, soft = soft)
  }

  # One 3D head bounding box, projected into every view, so the views share a
  # scale and the crosshair lines up; views are then padded to a common height.
  bbox <- if (isTRUE(do_crop)) volume_foreground_bbox(vol, overlay, ov_thresh) else NULL
  windows <- lapply(info, function(i) {
    full <- list(xlim = raster_extent_from_centers(i$oriented$x),
                 ylim = raster_extent_from_centers(i$oriented$y))
    if (is.null(bbox)) return(full)
    kept <- setdiff(seq_len(3L), i$along)
    corners <- as.matrix(expand.grid(bbox[1:2, kept[1]], bbox[1:2, kept[2]]))
    disp <- t(apply(corners, 1L, function(g) slice_grid_to_display(i$oriented, g)))
    half <- i$oriented$display_spacing / 2
    # Stay one voxel inside the image horizontally, so the raster always
    # covers the full width of the view (no 1-px step where padding meets it).
    inset <- 2 * half[1]
    list(xlim = c(max(min(disp[, 1]) - half[1], full$xlim[1] + inset),
                  min(max(disp[, 1]) + half[1], full$xlim[2] - inset)),
         ylim = c(max(min(disp[, 2]) - half[2], full$ylim[1]), min(max(disp[, 2]) + half[2], full$ylim[2])))
  })
  # A common height gives every view the same scale. Windows are grown
  # upward from the bottom of the field of view where needed: padding above
  # the head is black on black air (invisible), whereas padding below the
  # field of view would show as an empty band.
  extents <- lapply(info, function(i) raster_extent_from_centers(i$oriented$y))
  common_h <- max(vapply(windows, function(w) diff(w$ylim), numeric(1)))
  windows <- Map(function(w, ext) {
    lo <- mean(w$ylim) - common_h / 2
    lo <- max(lo, ext[1])
    w$ylim <- c(lo, lo + common_h)
    w
  }, windows, extents)

  world <- as.numeric(grid_to_coord(space(vol), matrix(coord, nrow = 1L)))
  plane_label <- c(axial = "z", coronal = "y", sagittal = "x")
  plane_value <- c(axial = world[3], coronal = world[2], sagittal = world[1])

  is_gray <- length(cmap) == 1L && tolower(cmap) %in% c("grays", "gray", "grey", "greys")
  make_panel <- function(name) {
    i <- info[[name]]
    w <- windows[[name]]
    layers <- if (!is.null(scale)) {
      list(overlay_slice_grob(overlay, i$z, i$along, scale, ov_thresh, alpha_fun,
                              alpha = ov_alpha, downsample = downsample))
    } else list()
    label <- sprintf("%s = %s mm", plane_label[[name]], format(round(plane_value[[name]]), trim = TRUE))
    p <- neuro_tile(i$oriented, bg_lim = lim, bg_cmap = cmap, layers = layers,
                    window = w, label = label,
                    orient = if (isTRUE(annotate)) orientation_letters(i$oriented, all = TRUE) else NULL,
                    tokens = tokens, interpolate = interp_bg,
                    floor_y = raster_extent_from_centers(i$oriented$y)[1],
                    # The inferior crop usually cuts through the neck; fade it
                    # out rather than end the image on a hard edge.
                    extra = if (!is.null(bbox) && name != "axial" && is_gray) {
                      list(fade_bottom_layer(w, resolve_cmap(cmap, 2L)[[1L]],
                                             floor_y = raster_extent_from_centers(i$oriented$y)[1]),
                           fade_top_layer(w, resolve_cmap(cmap, 2L)[[1L]],
                                          ceil_y = raster_extent_from_centers(i$oriented$y)[2]))
                    } else list())
    if (isTRUE(crosshair) && length(i$cross) == 2L) {
      p <- add_crosshair(p, i$cross, w, tokens)
    }
    p
  }
  plots <- lapply(names(planes), make_panel)
  names(plots) <- names(planes)

  show_cbar <- if (is.null(colorbar)) {
    !is.null(overlay) || !is_gray || cbar_title_given
  } else isTRUE(colorbar)
  cbar <- NULL
  if (show_cbar) {
    cbar <- if (!is.null(scale)) {
      neuro_colorbar(scale$lim, scale$pal, thresh = ov_thresh, diverging = scale$diverging,
                     tokens = tokens, title = cbar_title,
                     alpha_fun = if (ov_alpha_mode == "binary") NULL else alpha_fun,
                     alpha = ov_alpha, over_hi = scale$over_hi, over_lo = scale$over_lo)
    } else {
      neuro_colorbar(lim, resolve_cmap(cmap, 256), positions = FALSE,
                     tokens = tokens, title = cbar_title)
    }
  }
  # Conventional reading order: sagittal, coronal, axial.
  shown <- c("sagittal", "coronal", "axial")
  widths <- vapply(windows[shown], function(w) diff(w$xlim) / common_h, numeric(1))
  build <- function(cv) neuro_assemble(unname(plots[shown]), colorbar = cbar, tokens = tokens,
                        title = title, subtitle = subtitle, caption = caption,
                        widths = unname(widths), canvas = cv)
  cv0 <- canvas_size(canvas)
  fig <- build(cv0)
  if (is.null(canvas)) fig <- neuro_figure(fig, build, cv0, tokens)
  neuro_finish(fig, plots, draw = draw, assemble = assemble, title = title,
               subtitle = subtitle, caption = caption, style = style)
}

#' Voxel bounding box of the head (and supra-threshold overlay)
#'
#' @return A 2 x 3 matrix of (min, max) voxel indices per axis, padded by 5\%,
#'   or NULL when no foreground is found.
#' @keywords internal
#' @noRd
volume_foreground_bbox <- function(vol, overlay = NULL, ov_thresh = 0, margin = 0.13) {
  arr <- as.array(vol)
  v <- as.numeric(arr)
  thr <- foreground_threshold(v)
  if (!is.finite(thr)) return(NULL)
  # Frame the brain (bright tissue) plus a margin that takes in the scalp,
  # rather than the whole head: this keeps the neck out of the S-I extent.
  thr2 <- foreground_threshold(v[is.finite(v) & v > thr])
  fg <- arr > (if (is.finite(thr2)) thr2 else thr)
  fg[is.na(fg)] <- FALSE
  if (!is.null(overlay)) {
    ov <- as.array(overlay)
    extra <- is.finite(ov) & ov != 0 & abs(ov) >= ov_thresh
    fg <- fg | extra
  }
  if (!any(fg)) return(NULL)
  d <- dim(arr)
  vapply(seq_len(3L), function(ax) {
    prof <- apply(fg, ax, sum)
    r <- range(which(prof > 0.08 * max(prof)))
    pad <- ceiling(diff(r) * margin)
    c(max(1, r[1] - pad), min(d[[ax]], r[2] + pad))
  }, numeric(2))
}

#' Gradient layer fading the bottom of a tile into the tile colour
#' @keywords internal
#' @noRd
fade_bottom_layer <- function(window, colour, frac = 0.15, floor_y = -Inf) {
  # Fade from wherever the image actually ends (the field of view may stop
  # above the window's lower edge) over the bottom part of the view.
  h <- diff(window$ylim) * frac
  # Start a little below the image edge so the last (interpolated) row is
  # inside the ramp and no hairline survives.
  y0 <- max(window$ylim[1], floor_y - 0.02 * h)
  rgb <- grDevices::col2rgb(colour) / 255
  n <- 64L
  a <- seq(0, 1, length.out = n)   # row 1 (top) transparent, bottom opaque
  ras <- array(0, dim = c(n, 1L, 4L))
  ras[, , 1] <- rgb[1]; ras[, , 2] <- rgb[2]; ras[, , 3] <- rgb[3]
  ras[, , 4] <- a
  ggplot2::annotation_custom(grid::rasterGrob(ras, width = grid::unit(1, "npc"),
                                              height = grid::unit(1, "npc"),
                                              interpolate = TRUE),
                             xmin = window$xlim[1], xmax = window$xlim[2],
                             ymin = y0, ymax = y0 + h)
}

#' Short gradient hiding the image's top edge where the view is padded above it
#' @keywords internal
#' @noRd
fade_top_layer <- function(window, colour, ceil_y, frac = 0.04) {
  if (!is.finite(ceil_y) || ceil_y >= window$ylim[2]) return(NULL)
  h <- diff(window$ylim) * frac
  rgb <- grDevices::col2rgb(colour) / 255
  n <- 32L
  a <- rev(seq(0, 1, length.out = n))  # row 1 (top) opaque, bottom transparent
  ras <- array(0, dim = c(n, 1L, 4L))
  ras[, , 1] <- rgb[1]; ras[, , 2] <- rgb[2]; ras[, , 3] <- rgb[3]
  ras[, , 4] <- a
  ggplot2::annotation_custom(grid::rasterGrob(ras, width = grid::unit(1, "npc"),
                                              height = grid::unit(1, "npc"),
                                              interpolate = TRUE),
                             xmin = window$xlim[1], xmax = window$xlim[2],
                             ymin = ceil_y - h, ymax = ceil_y + 0.02 * h)
}

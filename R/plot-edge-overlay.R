#' Overlay fixed and moving edge maps on a background volume
#'
#' Displays a structural/reference background with two edge channels rendered in
#' distinct colors. This is intended for registration QC where fixed/template
#' edges and moving/result edges should coincide.
#'
#' @param bgvol Background 3D volume.
#' @param fixed_edges Edge map for the fixed/reference image on the same
#'   NeuroSpace grid as `bgvol`.
#' @param moving_edges Edge map for the moving/result image on the same
#'   NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot: indices along `along` (\code{unit = "index"})
#'   or world coordinates (\code{unit = "mm"}). \code{NULL} (default) picks
#'   \code{n_slices} slices spread over the brain.
#' @param along Native voxel-grid axis for slicing. Display orientation is
#'   inferred from the image affine.
#' @param bg_cmap Background palette.
#' @param fixed_color,moving_color Overlay colors for the two edge maps.
#' @param agree_color Colour for contour pixels present in both edge maps
#'   (where the images agree), keyed "both"; \code{NA} turns the agreement
#'   layer and its key entry off.
#' @param bg_range,edge_range "robust" or "data" intensity scaling.
#' @param probs Quantiles for robust scaling.
#' @param edge_thresh Values below this edge magnitude are transparent.
#' @param edge_alpha Global alpha for edge overlays.
#' @param ncol Number of columns (\code{NULL} = chosen to fill the canvas).
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if \code{TRUE}, also print the figure immediately (and
#'   return it invisibly). By default the figure is returned visibly.
#' @param style Visual style: \code{"light"}, \code{"dark"} or \code{"report"}.
#' @param labels Length-2 character vector naming the fixed and moving edge
#'   maps in the key (e.g. \code{c("MNI template", "registered EPI")}).
#' @param legend Logical; draw the colour key under the panels.
#' @param bg_dim Brightness multiplier (0--1) for the background so the edge
#'   colours read clearly.
#' @param crop Logical; crop panels to the head bounding box.
#' @param unit \code{"index"} or \code{"mm"}: how \code{zlevels} is
#'   interpreted. Panels are labelled in world coordinates.
#' @param annotate Logical; draw orientation letters on the first panel.
#' @param assemble Logical; return one \pkg{patchwork} figure (default) or the
#'   list of panel ggplots.
#' @param n_slices Number of automatically chosen slices.
#' @param thin Logical; thin edge bands to their ridge lines (1--2 pixels) so
#'   the two contours can be compared precisely. Default \code{TRUE}.
#' @param compute_edges Logical; if \code{TRUE}, \code{fixed_edges} and
#'   \code{moving_edges} are ordinary intensity images (e.g. the template and
#'   the registered image) and their edges are computed per slice as the
#'   in-plane gradient magnitude, keeping the strongest 12\% of head voxels.
#' @param canvas Optional \code{c(width, height)} in inches to fit the layout
#'   to (default: the open device).
#' @param interpolate Logical; smooth the background and anti-alias the edge
#'   contours (default \code{TRUE}).
#' @param focus_brain Logical; fade edges that lie outside the (dilated)
#'   brain mask of \code{bgvol} (bright tissue connected to the centre of the
#'   head, which excludes scalp fat), so scalp and skull contours do not
#'   dominate the comparison. Default \code{TRUE}.
#' @return A figure (class \code{neuro_fig}, a \pkg{patchwork} whose
#'   layout is re-fitted to the device it is drawn on) when
#'   \code{assemble = TRUE} or a named list
#'   of panel ggplots (\code{assemble = FALSE}); invisibly when
#'   \code{draw = TRUE}.
#' @examples
#' \donttest{
#' fixed <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
#' moving <- fixed
#' moving[] <- fixed[c(2:dim(fixed)[1], 1), , ]   # a one-voxel shift
#' p <- plot_edge_overlay(fixed, fixed, moving, compute_edges = TRUE,
#'                        labels = c("template", "registered"))
#' ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
#' }
#' @family plot_neuro
#' @export
plot_edge_overlay <- function(
    bgvol, fixed_edges, moving_edges, zlevels = NULL, along = 3L,
    bg_cmap = "grays", fixed_color = "#2bb8f0", moving_color = "#ff3b30",
    agree_color = "#ffd23f",
    bg_range = c("robust", "data"), edge_range = c("robust", "data"),
    probs = c(.02, .98), edge_thresh = 0, edge_alpha = .9,
    ncol = NULL, title = NULL, subtitle = NULL, caption = NULL, draw = FALSE,
    style = c("light", "dark", "report"), labels = c("fixed", "moving"),
    legend = TRUE, bg_dim = 0.75, crop = TRUE, unit = c("index", "mm"),
    annotate = TRUE, assemble = TRUE, n_slices = 6L, thin = TRUE,
    compute_edges = FALSE, canvas = NULL, interpolate = TRUE, focus_brain = TRUE) {
  unit_missing <- missing(unit)
  assert_same_neuro_grid(bgvol, fixed_edges = fixed_edges, moving_edges = moving_edges)
  style <- match_choice(style, c("light", "dark", "report"))
  unit <- match_choice(unit, c("index", "mm"), "unit")
  tokens <- neuro_style_tokens(style)
  labels <- check_labels2(labels)
  show_agree <- !(length(agree_color) == 1L && is.na(agree_color))
  if (show_agree) {
    ok <- tryCatch({ grDevices::col2rgb(agree_color); length(agree_color) == 1L },
                   error = function(e) FALSE)
    if (!ok) cli::cli_abort("{.arg agree_color} must be a single colour, or NA to not mark agreement.",
                            call = NULL)
  }
  check_count(n_slices, "n_slices")
  unit_hint(zlevels, unit_missing, bgvol)
  zlevels <- resolve_slice_levels(zlevels, bgvol, along, unit = unit, n = n_slices)
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), 1L)
  zlevels <- panel_args$zlevels
  along <- panel_args$along

  sel <- function(vol) unlist(lapply(zlevels, function(z) as.numeric(volume_slice_matrix(vol, z, along))))
  bg_vals <- sel(bgvol)
  bg_lim <- background_display_limits(bg_range, bg_vals, probs = probs)
  # Bright-tissue threshold (two-level Otsu) for the brain focus mask.
  brain_thr <- NA_real_
  if (isTRUE(focus_brain)) {
    t1 <- foreground_threshold(bg_vals)
    if (is.finite(t1)) brain_thr <- foreground_threshold(bg_vals[is.finite(bg_vals) & bg_vals > t1])
  }

  # Edge slices (computed from intensities if asked), thinned to ridges.
  edge_slice <- function(vol, z) {
    m <- volume_slice_matrix(vol, z, along = along)
    if (isTRUE(compute_edges)) m <- slice_gradient_magnitude(m)
    abs(m)
  }
  edge_cut <- 0
  if (isTRUE(compute_edges)) {
    g <- unlist(lapply(zlevels, function(z) c(as.numeric(edge_slice(fixed_edges, z)),
                                              as.numeric(edge_slice(moving_edges, z)))))
    g <- g[is.finite(g) & g > 0]
    if (length(g)) edge_cut <- stats::quantile(g, 0.88, names = FALSE)
  }
  prepare <- function(vol, z) {
    m <- edge_slice(vol, z)
    m[!is.finite(m) | m <= edge_cut] <- 0
    if (isTRUE(thin)) m <- thin_edges(m)
    m
  }
  edge_vals <- unlist(lapply(zlevels, function(z) c(as.numeric(prepare(fixed_edges, z)),
                                                    as.numeric(prepare(moving_edges, z)))))
  edge_vals <- edge_vals[is.finite(edge_vals) & edge_vals > edge_thresh]
  edge_lim <- if (length(edge_vals)) resolve_display_limits(edge_range, edge_vals, probs = probs) else c(0, 1)
  # Dim the background ramp so the two edge colours dominate.
  bg_pal <- resolve_cmap(bg_cmap, 256)
  bg_pal <- grDevices::adjustcolor(bg_pal, red.f = bg_dim, green.f = bg_dim, blue.f = bg_dim)

  edge_alpha_map <- function(mat) {
    vals <- abs(mat)
    alpha <- if (any(!is.finite(edge_lim)) || edge_lim[[1L]] == edge_lim[[2L]]) {
      as.numeric(vals > edge_thresh & vals > 0)
    } else {
      pmax(0, pmin(1, rescale01(as.numeric(vals), edge_lim)))^0.7
    }
    alpha <- matrix(alpha, nrow = nrow(mat), ncol = ncol(mat))
    alpha[!is.finite(alpha)] <- 0
    alpha[vals <= edge_thresh] <- 0
    alpha
  }
  edge_grob <- function(vol, z, colour, keep = NULL) {
    m <- prepare(vol, z)
    am <- edge_alpha_map(m)
    if (!is.null(keep)) am[!keep] <- 0
    if (is.finite(brain_thr)) {
      bg <- volume_slice_matrix(bgvol, z, along = along)
      head <- is.finite(bg) & bg > foreground_threshold(bg_vals)
      brain <- dilate_mask(brain_component(is.finite(bg) & bg > brain_thr, head), 3L)
      am[!brain] <- 0
    }
    o <- orient_volume_slice_for_raster(vol, z, along = along, mat = m,
                                        alpha_map = am)
    pos <- ifelse(is.finite(as.numeric(o$mat)), 1, NA_real_)
    grid::rasterGrob(positions_to_rgba(pos, nrow(o$mat), ncol(o$mat), c(colour, colour),
                                       alpha = edge_alpha, alpha_map = o$alpha_map),
                     interpolate = isTRUE(interpolate))
  }

  window <- if (isTRUE(crop)) foreground_crop_window(bgvol, zlevels, along) else NULL
  plots <- lapply(seq_along(zlevels), function(k) {
    z <- zlevels[[k]]
    o <- orient_volume_slice_for_raster(bgvol, z, along = along)
    neuro_tile(o, bg_lim = bg_lim, bg_cmap = bg_pal,
               layers = {
                 # Where the two contours coincide, draw them once in the
                 # agreement colour, so aligned boundaries stay visible
                 # instead of one channel hiding the other.
                 # Agreement is judged on the edges actually drawn (after the
                 # threshold), not on the raw maps.
                 both <- prepare(fixed_edges, z) > edge_thresh &
                   prepare(moving_edges, z) > edge_thresh
                 if (!show_agree) both[] <- FALSE
                 c(list(edge_grob(fixed_edges, z, fixed_color, keep = !both),
                        edge_grob(moving_edges, z, moving_color, keep = !both)),
                   if (show_agree) list(edge_grob(fixed_edges, z, agree_color, keep = both)))
               },
               window = window, label = slice_world_label(bgvol, z, along),
               orient = if (isTRUE(annotate) && k == 1L) orientation_letters(o) else NULL,
               tokens = tokens, interpolate = isTRUE(interpolate))
  })
  key <- if (isTRUE(legend)) {
    neuro_key_strip(Filter(Negate(is.null), list(list(fill = fixed_color, label = labels[[1L]]),
                         list(fill = moving_color, label = labels[[2L]]),
                         if (show_agree) list(fill = agree_color, label = "both"),
                         list(fill = NA, label = "contours should coincide where the images are aligned"))),
                    tokens)
  }
  o1 <- orient_volume_slice_for_raster(bgvol, zlevels[[1L]], along = along)
  build <- function(cv) neuro_assemble(plots, ncol = ncol, tile_aspect = tile_aspect_of(window, o1),
                        key = key, tokens = tokens, title = title,
                        subtitle = subtitle, caption = caption, canvas = cv)
  cv0 <- canvas_size(canvas)
  fig <- build(cv0)
  if (is.null(canvas)) fig <- neuro_figure(fig, build, cv0, tokens)
  neuro_finish(fig, plots, draw = draw, assemble = assemble, title = title,
               subtitle = subtitle, caption = caption, style = style,
               panel_names = vapply(zlevels, function(z) slice_world_label(bgvol, z, along), ""))
}

#' In-plane gradient magnitude of a slice (central differences)
#' @keywords internal
#' @noRd
slice_gradient_magnitude <- function(m) {
  m[!is.finite(m)] <- 0
  nr <- nrow(m); nc <- ncol(m)
  gx <- matrix(0, nr, nc); gy <- matrix(0, nr, nc)
  if (nr > 2L) gx[2:(nr - 1L), ] <- (m[3:nr, , drop = FALSE] - m[1:(nr - 2L), , drop = FALSE]) / 2
  if (nc > 2L) gy[, 2:(nc - 1L)] <- (m[, 3:nc, drop = FALSE] - m[, 1:(nc - 2L), drop = FALSE]) / 2
  sqrt(gx^2 + gy^2)
}

#' Thin an edge-magnitude slice to its ridge lines
#'
#' Keeps a pixel when it is a local maximum of the magnitude along at least
#' one image axis, which reduces thick gradient bands to 1--2 pixel contours.
#' @keywords internal
#' @noRd
thin_edges <- function(m) {
  nr <- nrow(m); nc <- ncol(m)
  if (nr < 3L || nc < 3L) return(m)
  pad <- matrix(0, nr + 2L, nc + 2L)
  pad[2:(nr + 1L), 2:(nc + 1L)] <- m
  up <- pad[1:nr, 2:(nc + 1L)]; down <- pad[3:(nr + 2L), 2:(nc + 1L)]
  left <- pad[2:(nr + 1L), 1:nc]; right <- pad[2:(nr + 1L), 3:(nc + 2L)]
  keep <- m > 0 & ((m >= up & m >= down) | (m >= left & m >= right))
  m[!keep] <- 0
  m
}

#' Dilate a logical matrix by r pixels (square structuring element)
#' @keywords internal
#' @noRd
dilate_mask <- function(m, r = 1L) {
  m[is.na(m)] <- FALSE
  nr <- nrow(m); nc <- ncol(m)
  out <- m
  for (k in seq_len(r)) {
    pad <- matrix(FALSE, nr + 2L, nc + 2L)
    pad[2:(nr + 1L), 2:(nc + 1L)] <- out
    out <- out | pad[1:nr, 2:(nc + 1L)] | pad[3:(nr + 2L), 2:(nc + 1L)] |
      pad[2:(nr + 1L), 1:nc] | pad[2:(nr + 1L), 3:(nc + 2L)]
  }
  out
}

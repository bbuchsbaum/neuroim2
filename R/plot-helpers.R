# Internal utilities used by plotting helpers

#' @keywords internal
#' @noRd
utils::globalVariables(c("x", "y", "value", "z", "fill", "xmin", "xmax", "ymin",
                         "ymax", "label", "hjust"))

#' Coerce a NeuroSlice (or matrix-like) to a numeric matrix
#' @keywords internal
#' @noRd
slice_to_matrix <- function(slc) {
  plain_matrix <- function(x) {
    dx <- dim(x)
    if (length(dx) != 2L) {
      cli::cli_abort("Expected a 2D object when coercing to a plain matrix.", call = NULL)
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

#' Extract a native voxel plane without constructing a NeuroSlice
#' @param vol A 3D NeuroVol.
#' @param z Slice index along \code{along}.
#' @param along Grid axis to hold fixed.
#' @keywords internal
#' @noRd
volume_slice_matrix <- function(vol, z, along = 3L) {
  along <- as.integer(along)
  z <- as.integer(z)
  retained <- setdiff(seq_len(3L), along)
  out <- switch(
    as.character(along),
    "1" = vol[z, , ],
    "2" = vol[, z, ],
    "3" = vol[, , z],
    cli::cli_abort("`along` must be one of 1, 2, or 3.", call = NULL)
  )
  matrix(
    as.numeric(out),
    nrow = dim(vol)[retained[[1L]]],
    ncol = dim(vol)[retained[[2L]]]
  )
}

#' Anatomically orient a regular native slice grid
#'
#' The axis metadata provides a signed permutation from the two voxel axes to
#' their nearest anatomical axes. Values are permuted and flipped so screen x/y
#' increase toward R/A/S, while coordinates remain a regular native-pixel grid.
#' This avoids feeding oblique or sheared world coordinates to geom_raster(),
#' which silently shifts pixels onto an axis-aligned grid.
#'
#' @param mat Numeric matrix in native slice-axis order.
#' @param axis_directions A 3 x 2 signed anatomical-axis matrix.
#' @param pixel_spacing Physical spacing for the two native slice axes.
#' @param alpha_map Optional matrix aligned with \code{mat}.
#' @param downsample Positive integer grid decimation factor.
#' @keywords internal
#' @noRd
orient_matrix_for_raster <- function(mat, axis_directions, pixel_spacing,
                                     alpha_map = NULL, downsample = 1L) {
  mat <- slice_to_matrix(mat)
  axis_directions <- as.matrix(axis_directions)
  if (!identical(dim(axis_directions), c(3L, 2L))) {
    cli::cli_abort("`axis_directions` must be a 3 x 2 matrix.", call = NULL)
  }
  pixel_spacing <- as.numeric(pixel_spacing)
  if (length(pixel_spacing) != 2L || any(!is.finite(pixel_spacing)) ||
      any(pixel_spacing <= 0)) {
    cli::cli_abort("`pixel_spacing` must contain two positive finite values.", call = NULL)
  }
  downsample <- as.integer(downsample)
  if (length(downsample) != 1L || is.na(downsample) || downsample < 1L) {
    cli::cli_abort("`downsample` must be a positive integer.", call = NULL)
  }
  if (!is.null(alpha_map)) {
    alpha_map <- slice_to_matrix(alpha_map)
    if (!identical(dim(alpha_map), dim(mat))) {
      cli::cli_abort("'alpha_map' must have the same dimensions as 'mat'.", call = NULL)
    }
  }

  # Each native voxel axis has one nearest anatomical direction. Using this
  # signed permutation, rather than raw affine x/y values, keeps oblique slices
  # on a regular raster without losing their anatomical display orientation.
  anatomical_axis <- apply(abs(axis_directions), 2L, which.max)
  strength <- vapply(
    seq_len(2L),
    function(j) abs(axis_directions[anatomical_axis[[j]], j]),
    numeric(1)
  )
  if (any(strength == 0) || anyDuplicated(anatomical_axis)) {
    cli::cli_abort("Slice axes must map to two distinct anatomical directions.", call = NULL)
  }

  display_axes <- sort(anatomical_axis)
  input_order <- match(display_axes, anatomical_axis)
  direction <- vapply(
    seq_len(2L),
    function(j) sign(axis_directions[display_axes[[j]], input_order[[j]]]),
    numeric(1)
  )

  canonical <- if (identical(input_order, c(1L, 2L))) mat else t(mat)
  canonical_alpha <- if (is.null(alpha_map)) NULL else {
    if (identical(input_order, c(1L, 2L))) alpha_map else t(alpha_map)
  }

  if (direction[[1L]] < 0) {
    canonical <- canonical[nrow(canonical):1L, , drop = FALSE]
    if (!is.null(canonical_alpha)) {
      canonical_alpha <- canonical_alpha[nrow(canonical_alpha):1L, , drop = FALSE]
    }
  }
  if (direction[[2L]] < 0) {
    canonical <- canonical[, ncol(canonical):1L, drop = FALSE]
    if (!is.null(canonical_alpha)) {
      canonical_alpha <- canonical_alpha[, ncol(canonical_alpha):1L, drop = FALSE]
    }
  }

  display_spacing <- pixel_spacing[input_order]
  full_display_dim <- dim(canonical)
  x_all <- (seq_len(nrow(canonical)) - 1) * display_spacing[[1L]]
  y_all <- (seq_len(ncol(canonical)) - 1) * display_spacing[[2L]]
  keep_x <- seq.int(1L, nrow(canonical), by = downsample)
  keep_y <- seq.int(1L, ncol(canonical), by = downsample)
  canonical <- canonical[keep_x, keep_y, drop = FALSE]
  if (!is.null(canonical_alpha)) {
    canonical_alpha <- canonical_alpha[keep_x, keep_y, drop = FALSE]
  }

  # rasterGrob and c(t(mat))-backed geom_raster data both expect row 1 at the
  # top. The canonical matrix stores x in rows and y in columns, so transpose it
  # and reverse y exactly once here.
  raster_mat <- t(canonical[, ncol(canonical):1L, drop = FALSE])
  raster_alpha <- if (is.null(canonical_alpha)) NULL else {
    t(canonical_alpha[, ncol(canonical_alpha):1L, drop = FALSE])
  }

  xvals <- x_all[keep_x]
  yvals <- rev(y_all[keep_y])
  xr <- raster_extent_from_centers(xvals)
  yr <- raster_extent_from_centers(yvals)
  negative_label <- c("L", "P", "I")
  positive_label <- c("R", "A", "S")
  normal_axis <- setdiff(seq_len(3L), display_axes)
  plane <- c("Sagittal", "Coronal", "Axial")[[normal_axis]]

  list(
    mat = raster_mat,
    alpha_map = raster_alpha,
    x = xvals,
    y = yvals,
    xmin = xr[[1L]],
    xmax = xr[[2L]],
    ymin = yr[[1L]],
    ymax = yr[[2L]],
    plane = plane,
    labels = c(
      left = negative_label[[display_axes[[1L]]]],
      right = positive_label[[display_axes[[1L]]]],
      bottom = negative_label[[display_axes[[2L]]]],
      top = positive_label[[display_axes[[2L]]]]
    ),
    input_order = input_order,
    flipped = direction < 0,
    display_dim = full_display_dim,
    display_spacing = display_spacing
  )
}

#' Orient a volume slice for regular raster rendering
#' @keywords internal
#' @noRd
orient_volume_slice_for_raster <- function(vol, z, along = 3L, mat = NULL,
                                           alpha_map = NULL, downsample = 1L) {
  along <- as.integer(along)
  retained <- setdiff(seq_len(3L), along)
  if (is.null(mat)) {
    mat <- volume_slice_matrix(vol, z, along)
  }
  directions <- perm_mat(axes(space(vol)))[, retained, drop = FALSE]
  orient_matrix_for_raster(
    mat,
    axis_directions = directions,
    pixel_spacing = spacing(space(vol))[retained],
    alpha_map = alpha_map,
    downsample = downsample
  )
}

#' Convert an oriented raster description to ggplot data
#' @keywords internal
#' @noRd
oriented_raster_df <- function(oriented) {
  df <- expand.grid(x = oriented$x, y = oriented$y)
  df$value <- c(t(oriented$mat))
  df
}

#' Map native two-axis grid coordinates into display coordinates
#' @keywords internal
#' @noRd
slice_grid_to_display <- function(oriented, grid) {
  grid <- as.numeric(grid)
  if (length(grid) != 2L || any(!is.finite(grid))) {
    cli::cli_abort("`grid` must contain two finite slice-grid coordinates.", call = NULL)
  }
  display_grid <- grid[oriented$input_order]
  display_grid[oriented$flipped] <-
    oriented$display_dim[oriented$flipped] + 1 - display_grid[oriented$flipped]
  (display_grid - 1) * oriented$display_spacing
}

#' Orient slice-aligned matrices for regular raster rendering
#' @param slc A NeuroSlice describing the 2D slice axes.
#' @param mat Numeric matrix aligned with \code{slc}.
#' @param alpha_map Optional numeric matrix aligned with \code{slc}.
#' @param downsample Positive integer grid decimation factor.
#' @keywords internal
#' @noRd
orient_slice_for_raster <- function(slc, mat, alpha_map = NULL, downsample = 1L) {
  orient_matrix_for_raster(
    mat,
    axis_directions = perm_mat(axes(space(slc))),
    pixel_spacing = spacing(space(slc)),
    alpha_map = alpha_map,
    downsample = downsample
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
    cli::cli_abort(sprintf("`%s` must be a 3D volume.", reference_name), call = NULL)
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
      cli::cli_abort(
        sprintf("`%s` must have the same dimensions as `%s`.", volume_name, reference_name),
        call = NULL
      )
    }
    if (!identical(space(volumes[[i]]), ref_space)) {
      cli::cli_abort(
        sprintf("`%s` must be on the same NeuroSpace grid as `%s`.", volume_name, reference_name),
        call = NULL
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
    cli::cli_abort("`along` must be one of 1, 2, or 3.", call = NULL)
  }

  if (is.numeric(zlevels) && any(is.finite(zlevels) & zlevels != round(zlevels))) {
    cli::cli_abort(c(
      "{.arg zlevels} must be whole slice indices.",
      "i" = "To give slice positions in world coordinates, use {.code unit = \"mm\"}."
    ), call = NULL)
  }
  zlevels <- as.integer(zlevels)
  if (!length(zlevels)) {
    cli::cli_abort("`zlevels` must contain at least one slice index.", call = NULL)
  }
  if (anyNA(zlevels) || any(zlevels < 1L | zlevels > dims[[along]])) {
    cli::cli_abort(c(
      sprintf("`zlevels` must be valid slice indices along axis %d.", along),
      "i" = sprintf("Valid indices are 1 to %d.", dims[[along]]),
      "i" = "To give slice positions in world coordinates, use {.code unit = \"mm\"}."
    ), call = NULL)
  }

  ncol <- as.integer(ncol)
  if (length(ncol) != 1L || is.na(ncol) || ncol < 1L) {
    cli::cli_abort("`ncol` must be a positive integer.", call = NULL)
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







#' Validate a colorbar title
#'
#' @param x Candidate title.
#' @return A length-1 character vector.
#' @keywords internal
#' @noRd
validate_cbar_title <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    cli::cli_abort("`cbar_title` must be a single non-NA character string.", call = NULL)
  }
  x
}


#' Self-tuning nonlinear opacity curve for statistical overlays
#'
#' Parameters of the curve used by \code{plot_overlay(ov_alpha_mode = "soft")}.
#' For magnitude \eqn{m}, opacity is
#' \deqn{alpha(m) = f + (1 - f)\,\mathrm{clamp}((m - lo) / (hi - lo), 0, 1)^{\gamma},}
#' where \eqn{f} is \code{alpha_floor}. Plotting then hides values below the
#' hard threshold and multiplies by \code{ov_alpha}.
#'
#' The knee \code{lo} defaults to the threshold or, when no threshold is set,
#' to the median non-zero magnitude (a robust noise-floor proxy). The cap
#' \code{hi} defaults to the largest magnitude. When \code{gamma} is not given
#' it is tuned so the median supra-knee magnitude maps to \code{alpha_mid}
#' (before the floor), and clamped to \code{[gamma_min, gamma_max]} so the
#' default curve stays convex; an explicit \code{gamma} bypasses the clamp.
#' Fix \code{knee}, \code{cap} and \code{gamma} to reuse one curve across
#' datasets.
#'
#' @param mags Numeric vector of overlay magnitudes (typically \code{abs(values)}).
#' @param thresh Hard threshold; used as the knee when \code{> 0} and
#'   \code{knee} is not given.
#' @param cap Optional upper magnitude anchor (opacity 1).
#' @param gamma Optional fixed exponent; \code{NULL} auto-tunes it.
#' @param alpha_mid Target opacity for the median supra-knee magnitude when
#'   \code{gamma} is tuned.
#' @param gamma_min,gamma_max Clamp range for the tuned exponent.
#' @param knee Optional non-negative lower magnitude anchor, overriding the
#'   threshold/median policy. Use \code{0} to ramp from zero.
#' @param alpha_floor Minimum opacity (0 to 1) above the knee, before
#'   \code{ov_alpha} and the hard threshold are applied.
#' @return A list with \code{lo}, \code{hi}, \code{gamma} and
#'   \code{alpha_floor}.
#' @examples
#' p <- soft_alpha_params(0:8, knee = 0, cap = 3, gamma = 0.7, alpha_floor = 0.15)
#' m <- c(0, 1.6, 2.1, 3, 8)
#' p$alpha_floor + (1 - p$alpha_floor) *
#'   pmin(pmax((m - p$lo) / (p$hi - p$lo), 0), 1)^p$gamma
#' @seealso \code{\link{plot_overlay}}
#' @export
soft_alpha_params <- function(mags, thresh = 0, cap = NULL, gamma = NULL,
                              alpha_mid = 0.2, gamma_min = 1.5, gamma_max = 5,
                              knee = NULL, alpha_floor = 0) {
  scalar <- function(x, name, lower, upper = Inf, strict = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        (if (strict) x <= lower else x < lower) || x > upper) {
      cli::cli_abort("{.arg {name}} must be a single finite number in the allowed range.",
                     call = NULL)
    }
  }
  scalar(thresh, "thresh", 0)
  scalar(alpha_floor, "alpha_floor", 0, 1)
  scalar(alpha_mid, "alpha_mid", 0, 1, strict = TRUE)
  if (alpha_mid >= 1) cli::cli_abort("{.arg alpha_mid} must be less than 1.", call = NULL)
  scalar(gamma_min, "gamma_min", 0, strict = TRUE)
  scalar(gamma_max, "gamma_max", gamma_min)
  if (!is.null(gamma)) scalar(gamma, "gamma", 0, strict = TRUE)
  if (!is.null(knee)) scalar(knee, "knee", 0)
  if (!is.null(cap)) scalar(cap, "cap", 0, strict = TRUE)
  explicit_knee <- !is.null(knee)
  mags <- mags[is.finite(mags) & mags > 0]
  knee <- if (explicit_knee) {
    knee
  } else if (isTRUE(thresh > 0)) {
    thresh
  } else if (length(mags)) {
    stats::median(mags)
  } else {
    0
  }
  hi <- if (!is.null(cap)) cap else if (length(mags)) max(mags) else knee + 1
  # An explicit cap below the resolved knee is an error. Equality is only an
  # error when the knee was also explicit: a data-driven cap that lands on a
  # derived knee (e.g. constant magnitudes) falls through to the rescue below.
  if (!is.null(cap) && (hi < knee || (explicit_knee && hi <= knee))) {
    cli::cli_abort("{.arg cap} must exceed {.arg knee}.", call = NULL)
  }
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
  list(lo = knee, hi = hi, gamma = gamma, alpha_floor = alpha_floor)
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
      cli::cli_abort("A numeric range must be two distinct finite values, c(lo, hi).", call = NULL)
    }
    return(range(range_arg))
  }
  mode <- match.arg(range_arg, c("robust", "data"))
  compute_limits(values, mode = mode, probs = probs)
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
    cli::cli_abort("`enhance` must be TRUE, FALSE, or a named list of enhance_stat_map() arguments.", call = NULL)
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
  labels <- switch(
    plane,
    axial = list(left = "L", right = "R", top = "A", bottom = "P"),
    coronal = list(left = "L", right = "R", top = "S", bottom = "I"),
    sagittal = list(left = "P", right = "A", top = "S", bottom = "I")
  )
  layers <- list(
    ggplot2::annotation_custom(grid::textGrob(labels$left,  gp = gp),
                               xmin = 0.5, xmax = 0.5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(labels$right, gp = gp),
                               xmin = nc + .5, xmax = nc + .5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(labels$top,   gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = .5, ymax = .5),
    ggplot2::annotation_custom(grid::textGrob(labels$bottom,gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = nr + .5, ymax = nr + .5)
  )
  layers
}

#' Otsu foreground threshold for a background image
#'
#' Two-class Otsu split of the finite, non-constant values of an image
#' (256-bin histogram). Used to separate head/brain from air so that display
#' windowing, cropping, and default slice selection ignore empty space.
#'
#' @param x Numeric vector of image values.
#' @return A single numeric threshold (or \code{NA} if undetermined).
#' @keywords internal
#' @noRd
foreground_threshold <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  if (length(x) > 2e5) x <- x[seq.int(1L, length(x), length.out = 2e5)]
  rng <- range(x)
  if (rng[1] == rng[2]) return(NA_real_)
  br <- seq(rng[1], rng[2], length.out = 257L)
  h <- graphics::hist(x, breaks = br, plot = FALSE)
  p <- h$counts / sum(h$counts)
  mids <- h$mids
  w0 <- cumsum(p)
  mu <- cumsum(p * mids)
  mu_t <- mu[length(mu)]
  between <- (mu_t * w0 - mu)^2 / (w0 * (1 - w0))
  between[!is.finite(between)] <- -Inf
  br[which.max(between) + 1L]
}

#' Display window for a structural background image
#'
#' For \code{range = "robust"} the window is computed over foreground voxels
#' only (values above the Otsu threshold) so that air does not drag the upper
#' limit down and over-expose tissue. The lower limit is the image minimum
#' (air renders black), the upper limit the \code{probs[2]} quantile of tissue.
#'
#' @keywords internal
#' @noRd
background_display_limits <- function(range_arg, values, probs = c(.02, .98)) {
  if (is.numeric(range_arg) || !identical(match.arg(range_arg[1], c("robust", "data")), "robust")) {
    return(resolve_display_limits(range_arg, values, probs = probs))
  }
  v <- values[is.finite(values)]
  if (!length(v)) return(c(0, 1))
  thr <- foreground_threshold(v)
  fg <- if (is.finite(thr)) v[v > thr] else v
  if (length(fg) < 10L) fg <- v
  lo <- stats::quantile(v, probs[1], names = FALSE)
  if (min(v) >= 0) lo <- min(v)
  # Map air (the non-foreground noise floor) to pure black, so it matches the
  # tile colour and padding never shows as a seam.
  if (is.finite(thr)) {
    air <- v[v <= thr]
    if (length(air) >= 10L) {
      lo <- max(lo, min(stats::quantile(air, 0.9, names = FALSE), lo + (thr - lo) / 4))
    }
  }
  # Upper limit a little above the brightest tissue so white matter renders
  # light grey with texture rather than a clipped white plateau.
  hi <- stats::quantile(fg, max(probs[2], 0.995), names = FALSE)
  hi <- lo + (hi - lo) * 1.08
  if (!is.finite(hi) || hi <= lo) hi <- max(v)
  if (hi <= lo) hi <- lo + 1
  c(lo, hi)
}

#' Choose informative default slices
#'
#' Picks \code{n} evenly spaced slices along \code{along} inside the extent of
#' the foreground (values above the Otsu threshold), trimmed by \code{trim} at
#' each end so the montage does not open or close on near-empty slices.
#'
#' @keywords internal
#' @noRd
default_slice_levels <- function(vol, along = 3L, n = 9L, trim = 0.06,
                                 support = NULL) {
  d <- dim(vol)[1:3]
  brain <- structural_slice_range(vol, along)
  if (!is.null(support) && any(support, na.rm = TRUE)) {
    support[is.na(support)] <- FALSE
    prof <- apply(support, along, sum)
    keep <- which(prof > 0.02 * max(prof))
    lo <- min(keep); hi <- max(keep)
    span <- hi - lo
    lo <- lo + 0.04 * span; hi <- hi - 0.04 * span
    # Stay inside the substantial-brain range when the overlay allows it, so
    # the default grid does not open on orbits and skull base; clusters
    # outside it remain reachable through zlevels.
    if (!is.null(brain)) {
      ilo <- max(lo, brain[1]); ihi <- min(hi, brain[2])
      if (ihi - ilo >= 0.5 * (hi - lo)) { lo <- ilo; hi <- ihi }
    }
  } else if (!is.null(brain)) {
    lo <- brain[1]; hi <- brain[2]
  } else {
    lo <- 1 + trim * (d[[along]] - 1); hi <- d[[along]] - trim * (d[[along]] - 1)
  }
  unique(as.integer(round(seq(lo, hi, length.out = n))))
}

#' Slice range covering substantial brain tissue along an axis
#'
#' Two-level Otsu: the first split separates head from air, the second
#' separates bright tissue (brain parenchyma on T1) from scalp and CSF. The
#' range starts where bright tissue is substantial (skipping neck and skull
#' base) and runs to near the vertex, trimmed more at the inferior end.
#' @return c(lo, hi) slice positions, or NULL.
#' @keywords internal
#' @noRd
structural_slice_range <- function(vol, along) {
  arr <- as.array(vol)
  v <- as.numeric(arr)
  thr <- foreground_threshold(v)
  if (!is.finite(thr)) return(NULL)
  thr2 <- foreground_threshold(v[is.finite(v) & v > thr])
  use <- if (is.finite(thr2)) thr2 else thr
  fg <- arr > use
  fg[is.na(fg)] <- FALSE
  prof <- apply(fg, along, sum)
  if (!any(prof > 0)) return(NULL)
  big <- which(prof > 0.35 * max(prof))
  any_t <- which(prof > 0.12 * max(prof))
  lo <- min(big); hi <- max(any_t)
  span <- hi - lo
  c(lo + 0.2 * span, hi - 0.03 * span)
}

#' World coordinate (mm) of a slice along a native axis
#'
#' Returns the world coordinate of the slice centre along the anatomical axis
#' nearest to the native slicing axis, with a label such as \code{"z = 24"}.
#'
#' @keywords internal
#' @noRd
slice_world_label <- function(vol, z, along = 3L) {
  sp <- space(vol)
  d <- dim(vol)[1:3]
  g <- (d + 1) / 2
  g[[along]] <- z
  w <- as.numeric(grid_to_coord(sp, matrix(g, nrow = 1L)))
  directions <- perm_mat(axes(sp))
  world_axis <- which.max(abs(directions[, along]))
  val <- w[[world_axis]]
  paste0(c("x", "y", "z")[[world_axis]], " = ", format(round(val), trim = TRUE), " mm")
}

#' Display limits for a statistical overlay
#'
#' Unlike a structural image, a statistical map is mostly zeros (outside the
#' mask) or near-zero noise. For \code{range = "robust"} the limits are therefore
#' computed over the non-zero finite values only, with the upper end at the
#' \code{max(probs[2], 0.99)} quantile of the magnitude, and they are widened if
#' necessary so that the threshold always lies inside the scale.
#'
#' @keywords internal
#' @noRd
overlay_display_limits <- function(range_arg, values, probs = c(.02, .98),
                                   thresh = 0) {
  if (is.numeric(range_arg)) {
    return(resolve_display_limits(range_arg, values, probs = probs))
  }
  mode <- match.arg(range_arg[1], c("robust", "data"))
  v <- values[is.finite(values) & values != 0]
  if (!length(v)) return(c(0, 1))
  if (mode != "data" && length(v) > 5e5) {
    # Quantiles of a large map are stable on a regular subsample.
    v <- v[seq.int(1L, length(v), length.out = 5e5)]
  }
  if (mode == "data") {
    lim <- range(v)
  } else {
    hi_p <- max(probs[2], 0.99)
    # With a threshold, the scale is set by the supra-threshold values (the
    # ones actually drawn); otherwise by all non-zero values.
    vv <- if (isTRUE(thresh > 0) && sum(abs(v) >= thresh) >= 10L) v[abs(v) >= thresh] else v
    lim <- c(stats::quantile(vv, 1 - hi_p, names = FALSE),
             stats::quantile(vv, hi_p, names = FALSE))
    if (min(vv) >= 0) lim[1] <- min(0, lim[1])
    if (lim[1] == lim[2]) lim <- range(v)
  }
  if (isTRUE(thresh > 0)) {
    if (lim[2] > 0 && lim[2] < thresh) lim[2] <- max(max(v), thresh * 1.05)
    if (lim[1] < 0 && lim[1] > -thresh) lim[1] <- min(min(v), -thresh * 1.05)
  }
  if (lim[1] == lim[2]) lim <- lim + c(-0.5, 0.5)
  lim
}

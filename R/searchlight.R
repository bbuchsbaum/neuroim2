#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Searchlight Analysis Methods
#'
#' @name searchlight-methods
#' @description Methods for performing searchlight analyses on neuroimaging data
NULL

#' @importFrom assertthat assert_that
#' @importFrom purrr map map_dbl map_int
#' @importFrom stats kmeans
#' @importFrom dbscan frNN
#' @importFrom utils object.size
NULL


#' Create a spherical random searchlight iterator
#'
#' @description
#' This function generates a spherical random searchlight iterator for analyzing
#' local neighborhoods of voxels within a given radius in a brain mask.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius of the searchlight sphere in voxel units.
#' @param nonzero Logical; if \code{TRUE} (default) discard zero-valued voxels in
#'   the mask when forming each searchlight.
#'
#' @return A list of \code{\linkS4class{ROIVolWindow}} objects, each representing
#'   a spherical searchlight region.
#'
#' @examples
#' # Create a simple brain mask
#' mask_data <- array(TRUE, c(10, 10, 10))
#' mask_data[1, 1, 1] <- FALSE
#' mask <- LogicalNeuroVol(mask_data, NeuroSpace(c(10,10,10)))
#'
#' # Generate random searchlight iterator with a radius of 2 voxels
#' \donttest{
#' searchlights <- random_searchlight(mask, radius = 6)
#' }
#'
#' @export
random_searchlight <- function(mask, radius, nonzero = TRUE) {
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  assert_that(radius > 0,
              msg = "radius must be positive")
  assert_that(is.logical(nonzero), length(nonzero) == 1,
              msg = "nonzero must be TRUE or FALSE")

  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  # Logical vector tracking remaining voxels
  remaining <- rep(TRUE, length(mask.idx))

  # Lookup array: maps voxel coords to index in mask.idx
  lookup <- array(0, dim(mask))
  lookup[mask.idx] <- seq_along(mask.idx)

  # Preallocate a result list
  slist <- vector("list", length(mask.idx))
  counter <- 1

  # Vector of voxel indices that remain
  remain_indices <- seq_along(mask.idx)

  while (length(remain_indices) > 0) {
    # sample a center index from remain_indices
    sel <- sample.int(length(remain_indices), 1)
    center_idx <- remain_indices[sel]
    center_coord <- grid[center_idx, , drop=FALSE]

    # Compute spherical ROI
    search <- spherical_roi(mask, center_coord, radius, nonzero = nonzero)
    vox <- coords(search)

    # If no voxels in ROI, remove center_idx to avoid infinite loop
    if (nrow(vox) == 0) {
      # Mark center voxel as used to progress
      remaining[center_idx] <- FALSE
      remain_indices <- remain_indices[remaining[remain_indices]]
      # continue to next iteration without adding to slist
      next
    }

    # Lookup voxel indices
    idx_lookup <- lookup[vox]
    mask_hits <- idx_lookup > 0

    active_mask <- rep(FALSE, nrow(vox))
    active_mask[mask_hits] <- remaining[idx_lookup[mask_hits]]

    # If none of the masked voxels are available, drop current center and move on
    if (!any(active_mask)) {
      remaining[center_idx] <- FALSE
      remain_indices <- remain_indices[remaining[remain_indices]]
      next
    }

    # Keep masked voxels that are still available;
    # optionally keep out-of-mask voxels when nonzero = FALSE
    kept <- active_mask | (!mask_hits & !nonzero)
    kept_vox <- vox[kept, , drop = FALSE]

    # Row position of the center voxel inside the kept ROI coordinates
    center_row <- which(rowSums(kept_vox == matrix(center_coord, nrow(kept_vox), 3, byrow = TRUE)) == 3)
    center_row <- if (length(center_row) == 0) NA_integer_ else center_row[1]

    parent_idx <- grid_to_index(mask, center_coord)

    search2 <- new("ROIVolWindow",
                   rep(1, nrow(kept_vox)),
                   space=space(mask),
                   coords=kept_vox,
                   center_index=as.integer(center_row),
                   parent_index=as.integer(parent_idx))

    # Expose row index within the ROI coordinates for downstream consumers
    attr(search2, "center_row_index") <- as.integer(center_row)
    # Index of the center voxel within the mask's nonzero ordering
    attr(search2, "mask_index") <- as.integer(center_idx)

    # Mark chosen voxels (that are in the mask) as used
    idx_keep <- idx_lookup[mask_hits & active_mask]
    remaining[idx_keep] <- FALSE

    # Update remain_indices to reflect removed voxels
    remain_indices <- remain_indices[remaining[remain_indices]]

    slist[[counter]] <- search2
    counter <- counter + 1
  }

  # Trim slist to the actual number of used slots
  slist[seq_len(counter - 1)]
}

#
# random_searchlight2 <- function(mask, radius) {
#   stopifnot(inherits(mask, "NeuroVol"))
#
#   mask.idx <- which(mask != 0)
#
#   n_voxels <- length(mask.idx)
#
#   remaining <- rep(TRUE, n_voxels)
#
#   grid <- index_to_grid(mask, mask.idx)
#
#   lookup <- integer(prod(dim(mask)))
#   lookup[mask.idx] <- seq_len(n_voxels)
#
#   slist <- list()
#   counter <- 1
#
#   while (any(remaining)) {
#     # Select a random index among the remaining indices
#     remaining_indices <- which(remaining)
#     center_idx_in_remaining <- sample(remaining_indices, 1)
#
#     center_idx <- mask.idx[center_idx_in_remaining]
#     center_coord <- grid[center_idx_in_remaining, ]
#
#     # Get a searchlight surrounding the center
#     search <- spherical_roi(mask, center_coord, radius, nonzero = TRUE)
#
#     # Get the indices of the voxels in the searchlight
#     vox_coords <- coords(search)
#     vox_indices <- lookup[grid_to_index(mask, vox_coords)]
#     vox_indices <- vox_indices[vox_indices > 0]
#
#     # Keep only the voxels that are still remaining
#     in_remaining <- remaining[vox_indices]
#     vox_indices <- vox_indices[in_remaining]
#     vox_coords <- vox_coords[in_remaining, , drop = FALSE]
#
#     if (length(vox_indices) > 0) {
#       # Mark these voxels as no longer remaining
#       remaining[vox_indices] <- FALSE
#
#       # Create a new ROIVolWindow
#       search2 <- new("ROIVolWindow", rep(1, length(vox_indices)), space=space(mask),
#                               coords=vox_coords, parent_index=as.integer(1), center_index=as.integer(center_idx))
#
#       slist[[counter]] <- search2
#       counter <- counter + 1
#     }
#   }
#
#   slist
# }

#' Create a resampled searchlight iterator
#'
#' @description
#' This function generates a resampled searchlight iterator by sampling regions
#' from within a brain mask. By default it builds spherical searchlights, but
#' users can provide a custom \code{shape_fun} to return ellipsoids, cubes, or
#' arbitrary irregular searchlight shapes. Centers are drawn with replacement,
#' so the same voxel (and its neighborhood) may appear multiple times. Each
#' searchlight can also draw its radius from a user-specified set of radii.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric scalar or vector specifying candidate radii (in voxel
#'   units) for the searchlight sphere. If a vector is supplied, a radius is
#'   sampled uniformly (with replacement) for each searchlight. All radii must
#'   be positive. Default is 8.
#' @param iter An integer specifying the total number of searchlights to sample
#'   (with replacement). Default is 100.
#' @param shape_fun Either \code{NULL} (default spherical kernel), a character
#'   keyword (\code{"sphere"}, \code{"ellipsoid"}, \code{"cube"}, \code{"blobby"}),
#'   or a custom function. Custom functions are called as
#'   \code{shape_fun(mask, center, radius, iter, nonzero)} and must return either
#'   a \code{\linkS4class{ROIVolWindow}} or an \code{n x 3} integer matrix of
#'   voxel coordinates. This enables anisotropic or irregular searchlights.
#' @param nonzero Logical; if \code{TRUE} (default), the generated searchlight is
#'   intersected with the non-zero voxels of \code{mask}. Applies to both the
#'   default sphere and any \code{shape_fun} that returns coordinates.
#'
#' @return A \code{deferred_list} object containing \code{\linkS4class{ROIVolWindow}}
#'   objects, each representing a sampled searchlight region drawn from within the mask.
#'
#' @details
#' Searchlight centers are sampled with replacement, so the same center (and its
#' surrounding voxels) can be selected multiple times. When multiple radii are
#' provided, each searchlight independently samples one radius from the supplied
#' values. Supplying \code{shape_fun} lets you draw non-spherical searchlights
#' (e.g., ellipsoids, cubes, blobby deformations, or task-specific kernels).
#' Built-in shortcuts are available via \code{shape_fun = "ellipsoid"},
#' \code{"cube"}, and \code{"blobby"}; \code{"sphere"} or \code{NULL} uses the
#' default spherical kernel.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate a resampled searchlight iterator with radii drawn from {4,6,8}
#' searchlights <- resampled_searchlight(mask, radius = c(4, 6, 8))
#'
#' # Use a custom shape: random ellipsoid scaled along each axis
#' ellipsoid_fun <- function(mask, center, radius, iter, nonzero) {
#'   scales <- runif(3, 0.5, 1.5)        # axis-wise stretch/compress
#'   vox <- spherical_roi(mask, center, radius, nonzero = FALSE)@coords
#'   ctr_mat <- matrix(center, nrow(vox), 3, byrow = TRUE)
#'   keep <- rowSums(((vox - ctr_mat) * scales)^2) <= radius^2
#'   vox[keep, , drop = FALSE]
#' }
#' ellip_searchlights <- resampled_searchlight(mask, radius = c(4, 6),
#'                                             iter = 50, shape_fun = ellipsoid_fun)
#'
#' # Or use built-in named shapes
#' ellip_builtin <- resampled_searchlight(mask, radius = 6, shape_fun = "ellipsoid")
#' cube_builtin  <- resampled_searchlight(mask, radius = 6, shape_fun = "cube")
#' 
#'
#' @name resampled_searchlight
#' @aliases resampled_searchlight bootstrap_searchlight
#' @rdname resampled_searchlight
#' @export
resampled_searchlight <- function(mask,
                                  radius = 8,
                                  iter = 100,
                                  shape_fun = NULL,
                                  nonzero = TRUE) {
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  assert_that(is.numeric(radius), length(radius) > 0, all(radius > 0),
              msg = "radius must be positive numeric")
  assert_that(length(iter) == 1, is.finite(iter), iter > 0,
              msg = "iter must be a positive number")
  assert_that(
    is.null(shape_fun) ||
      is.function(shape_fun) ||
      (is.character(shape_fun) && length(shape_fun) == 1),
    msg = "shape_fun must be NULL, a function, or one of the built-in shape names"
  )
  assert_that(is.logical(nonzero), length(nonzero) == 1,
              msg = "nonzero must be TRUE or FALSE")

  iter <- as.integer(iter)

  mask.idx <- which(mask != 0)
  assert_that(length(mask.idx) > 0,
              msg = "mask contains no nonzero voxels to sample")

  grid <- index_to_grid(mask, mask.idx)

  # Sample centers with replacement; allows iter > nrow(grid)
  center_idx <- sample.int(nrow(grid), iter, replace = TRUE)

  # Sample radius per iteration when a vector is supplied
  radii <- if (length(radius) == 1L) rep(radius, iter) else sample(radius, iter, replace = TRUE)

  # Resolve built-in shape keywords to concrete functions
  if (is.character(shape_fun)) {
    shape_fun <- match.arg(shape_fun, c("sphere", "ellipsoid", "cube", "blobby"))
    shape_fun <- switch(shape_fun,
                        sphere    = NULL,
                        ellipsoid = ellipsoid_shape(),
                        cube      = cube_shape(),
                        blobby    = blobby_shape())
  }

  force(mask)
  # helper to coerce arbitrary shape output into a ROIVolWindow
  to_roi_window <- function(obj, center_coord) {
    center_vec <- drop(center_coord)
    parent_index <- grid_to_index(mask, center_vec)

    # User supplied a ROIVolWindow already; optionally filter nonzero voxels
    if (inherits(obj, "ROIVolWindow")) {
      coords <- obj@coords
      vals   <- obj@.Data

      if (nonzero) {
        keep <- mask[coords] != 0
        coords <- coords[keep, , drop = FALSE]
        vals   <- vals[keep]
      }

      # recompute center row in case filtering changed indexing
      center_row <- if (nrow(coords) == 0) {
        NA_integer_
      } else {
        which(rowSums(coords == matrix(center_vec, nrow(coords), 3, byrow = TRUE)) == 3)
      }
      center_row <- if (length(center_row) == 0) NA_integer_ else center_row[1]

      return(new("ROIVolWindow",
                 vals,
                 space = space(mask),
                 coords = coords,
                 center_index = as.integer(center_row),
                 parent_index = as.integer(parent_index)))
    }

    # Allow a bare matrix of voxel coordinates (integer, 3 cols)
    if (is.matrix(obj) && ncol(obj) == 3) {
      coords <- obj

      # keep in-bounds coordinates
      vdim <- dim(mask)
      in_bounds <- coords[,1] >= 1 & coords[,1] <= vdim[1] &
                   coords[,2] >= 1 & coords[,2] <= vdim[2] &
                   coords[,3] >= 1 & coords[,3] <= vdim[3]
      coords <- coords[in_bounds, , drop = FALSE]

      # apply mask filtering if requested
      if (nonzero && nrow(coords) > 0) {
        keep <- mask[coords] != 0
        coords <- coords[keep, , drop = FALSE]
      }

      if (nrow(coords) == 0) {
        return(new("ROIVolWindow",
                   numeric(0),
                   space = space(mask),
                   coords = matrix(ncol = 3, nrow = 0),
                   center_index = as.integer(NA),
                   parent_index = as.integer(parent_index)))
      }

      # identify where (if anywhere) the sampled center lives in coords
      center_row <- which(rowSums(coords == matrix(center_coord, nrow(coords), 3, byrow = TRUE)) == 3)
      center_row <- if (length(center_row) == 0) NA_integer_ else center_row[1]

      return(new("ROIVolWindow",
                 rep(1, nrow(coords)),
                 space = space(mask),
                 coords = coords,
                 center_index = as.integer(center_row),
                 parent_index = as.integer(parent_index)))
    }

    stop("shape_fun must return a ROIVolWindow or an n x 3 matrix of voxel coordinates")
  }

  f <- function(i) {
    ctr <- grid[center_idx[i], , drop = FALSE]
    rad <- radii[i]

    roi_obj <- if (is.null(shape_fun)) {
      spherical_roi(mask, ctr, rad, nonzero = nonzero)
    } else {
      shape_fun(mask = mask, center = ctr, radius = rad, iter = i, nonzero = nonzero)
    }

    out <- to_roi_window(roi_obj, center_coord = ctr)
    attr(out, "mask_index") <- as.integer(center_idx[i])
    out
  }

  deflist::deflist(f, iter)
}

#' @rdname resampled_searchlight
#' @export
bootstrap_searchlight <- function(mask, radius=8, iter=100) {
  .Deprecated("resampled_searchlight", package = "neuroim2")
  resampled_searchlight(mask, radius=radius, iter=iter)
}

#' Convenience shape generators for \code{resampled_searchlight()}
#'
#' Helpers that return ready-to-use \code{shape_fun} callbacks for
#' \code{resampled_searchlight()}, covering a few sensible non-spherical
#' defaults.
#'
#' Each returned function has signature \code{function(mask, center, radius, iter, nonzero)}
#' and should return an \eqn{n \times 3} integer coordinate matrix. The
#' coordinates are later converted to a \code{ROIVolWindow} internally.
#'
#' @param scales Length-3 positive numeric vector scaling the x/y/z axes relative
#'   to a sphere (for \code{ellipsoid_shape}). Values >1 stretch; <1 compress.
#' @param jitter Non-negative numeric; standard deviation of multiplicative
#'   Gaussian noise applied to \code{scales} each draw (ellipsoid).
#' @param drop Numeric in [0,1]; probability of dropping a voxel (blobby).
#' @param edge_fraction Numeric in (0,1]; fraction of farthest voxels (by
#'   Euclidean distance from the center, in voxel units) considered "edge" and
#'   eligible for random dropping (blobby).
#'
#' @return A function suitable for the \code{shape_fun} argument of
#'   \code{resampled_searchlight()}.
#'
#' @examples
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Ellipsoid stretched along z with modest per-iteration jitter
#' sl_ellip <- resampled_searchlight(mask, radius = 6,
#'                                    shape_fun = ellipsoid_shape(scales = c(1, 1, 1.4),
#'                                                               jitter = 0.1))
#'
#' # Simple axis-aligned cube (Chebyshev ball)
#' sl_cube <- resampled_searchlight(mask, radius = 5, shape_fun = "cube")
#'
#' # Blobby sphere with 40% dropout on boundary voxels
#' sl_blob <- resampled_searchlight(mask, radius = 6,
#'                                  shape_fun = blobby_shape(drop = 0.4, edge_fraction = 0.6))
#'
#' @name searchlight_shape_functions
NULL

#' @rdname searchlight_shape_functions
#' @export
ellipsoid_shape <- function(scales = c(1, 1, 1), jitter = 0) {
  assert_that(is.numeric(scales), length(scales) == 3, all(scales > 0))
  assert_that(is.numeric(jitter), length(jitter) == 1, jitter >= 0)

  function(mask, center, radius, iter, nonzero) {
    vox <- spherical_roi(mask, center, radius, nonzero = FALSE)@coords
    if (nrow(vox) == 0) {
      return(matrix(ncol = 3, nrow = 0))
    }

    # Optionally jitter axis scales per draw
    sc <- scales
    if (jitter > 0) {
      sc <- pmax(scales * (1 + stats::rnorm(3, 0, jitter)), .Machine$double.eps)
    }

    ctr_mat <- matrix(center, nrow(vox), 3, byrow = TRUE)
    keep <- rowSums(((vox - ctr_mat) * sc)^2) <= radius^2
    vox[keep, , drop = FALSE]
  }
}

#' @rdname searchlight_shape_functions
#' @export
cube_shape <- function() {
  function(mask, center, radius, iter, nonzero) {
    spc <- spacing(mask)
    half_width <- ceiling(radius / spc)

    coords <- as.matrix(expand.grid(
      seq.int(center[1] - half_width[1], center[1] + half_width[1]),
      seq.int(center[2] - half_width[2], center[2] + half_width[2]),
      seq.int(center[3] - half_width[3], center[3] + half_width[3])
    ))

    coords
  }
}

#' @rdname searchlight_shape_functions
#' @export
blobby_shape <- function(drop = 0.3, edge_fraction = 0.7) {
  assert_that(is.numeric(drop), length(drop) == 1, drop >= 0, drop <= 1)
  assert_that(is.numeric(edge_fraction), length(edge_fraction) == 1,
              edge_fraction > 0, edge_fraction <= 1)

  function(mask, center, radius, iter, nonzero) {
    roi <- spherical_roi(mask, center, radius, nonzero = FALSE)
    coords <- roi@coords
    if (nrow(coords) == 0) {
      return(coords)
    }

    ctr_mat <- matrix(center, nrow(coords), 3, byrow = TRUE)
    dist <- sqrt(rowSums((coords - ctr_mat)^2))

    # Identify edge voxels; drop a fraction of them at random
    threshold <- stats::quantile(dist, probs = edge_fraction)
    edge_idx <- dist >= threshold
    drop_mask <- edge_idx & stats::runif(length(dist)) < drop
    keep <- !drop_mask

    coords[keep, , drop = FALSE]
  }
}


#' Create an exhaustive searchlight iterator for voxel coordinates using spherical_roi
#'
#' @description
#' This function generates an exhaustive searchlight iterator that returns voxel
#' coordinates for each searchlight sphere within the provided mask, using
#' `spherical_roi` for neighborhood computation. The iterator
#' visits every non-zero voxel in the mask as a potential center voxel.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius (in mm) of the spherical searchlight.
#' @param nonzero A logical value indicating whether to include only coordinates
#'   with nonzero values in the supplied mask. Default is FALSE.
#' @param cores An integer specifying the number of cores to use for parallel
#'   computation. Default is 0, which uses a single core.
#'
#' @return A \code{deferred_list} object containing matrices of integer-valued
#'   voxel coordinates, each representing a searchlight region.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate an exhaustive searchlight iterator with a radius of 6 mm
#'
#' searchlights <- searchlight_coords(mask, radius = 6)
#' 
#'
#' @export
searchlight_coords <- function(mask, radius, nonzero=FALSE, cores=0) {
  # Decide which voxels to consider
  if (nonzero) {
    mask.idx <- which(mask != 0)
  } else {
    mask.idx <- seq_len(prod(dim(mask)))
  }

  # Convert voxel indices to coordinates
  grid <- index_to_grid(mask, mask.idx) # Nx3 integer voxel coords

  # Define a function to get the spherical neighborhood for a single voxel
  f <- function(i) {
    centroid <- grid[i, , drop=FALSE]
    roi <- spherical_roi(mask, centroid, radius=radius, nonzero=nonzero)
    coords(roi) # returns an Nx3 matrix of voxel coordinates
  }

  # Create a deferred_list for lazy evaluation
  deflist::deflist(f, length(mask.idx))
}


#' Create an exhaustive searchlight iterator
#'
#' @description
#' This function generates an exhaustive searchlight iterator that returns either
#' voxel coordinates or ROIVolWindow objects for each searchlight sphere within
#' the provided mask. The iterator visits every non-zero voxel in the mask as a
#' potential center voxel.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius (in mm) of the spherical searchlight.
#' @param eager A logical value specifying whether to eagerly compute the
#'   searchlight ROIs. Default is FALSE, which uses lazy evaluation.
#' @param nonzero A logical value indicating whether to include only coordinates
#'   with nonzero values in the supplied mask. Default is FALSE.
#' @param cores An integer specifying the number of cores to use for parallel
#'   computation. Default is 0, which uses a single core.
#'
#' @return A \code{deferred_list} object containing either matrices of integer-valued
#'   voxel coordinates or \code{\linkS4class{ROIVolWindow}} objects, each representing
#'   a searchlight region.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate an exhaustive searchlight iterator with a radius of 6 mm
#' 
#' searchlights <- searchlight(mask, radius = 6, eager = FALSE)
#' 
#'
#' @export
#' @rdname searchlight
searchlight <- function(mask, radius, eager=FALSE, nonzero=FALSE, cores=0) {
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  assert_that(radius > 0,
              msg = "radius must be positive")
  assert_that(is.logical(eager),
              msg = "eager must be TRUE or FALSE")
  assert_that(is.logical(nonzero),
              msg = "nonzero must be TRUE or FALSE")
  assert_that(cores >= 0,
              msg = "cores must be non-negative")

  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  if (!eager) {
    force(mask)
    force(radius)
    f <- function(i) { 
      roi <- spherical_roi(mask, grid[i,], radius, nonzero=nonzero)
      attr(roi, "mask_index") <- as.integer(i)
      roi
    }
    deflist::deflist(f, nrow(grid))
  } else {
    # Use spherical_roi_set to get all ROIs at once
    result_list <- spherical_roi_set(
      bvol = mask,
      centroids = grid,
      radius = radius,
      nonzero = nonzero
    )

    for (i in seq_along(result_list)) {
      attr(result_list[[i]], "mask_index") <- as.integer(i)
    }

    result_list
  }
}


#' Create a clustered searchlight iterator
#'
#' @description
#' This function generates a searchlight iterator that iterates over successive
#' spatial clusters in an image volume. It allows for the exploration of spatially
#' clustered regions within the provided mask by using either a pre-defined
#' clustered volume or performing k-means clustering to generate the clusters.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param cvol An optional \code{ClusteredNeuroVol} instance representing pre-defined
#'   clusters within the mask. If provided, the 'csize' parameter is ignored.
#' @param csize An optional integer specifying the number of clusters to be
#'   generated using k-means clustering (ignored if \code{cvol} is provided).
#'
#' @return A \code{deferred_list} object containing \code{ROIVol} objects, each
#'   representing a clustered region within the image volume.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate a clustered searchlight iterator with 5 clusters
#' clust_searchlight <- clustered_searchlight(mask, csize = 5)
#' 
#'
#' @rdname clustered_searchlight
#' @export
clustered_searchlight <- function(mask, cvol=NULL, csize=NULL) {
  assert_that(!is.null(csize) || !is.null(cvol),
              msg = "must provide either 'cvol' or 'csize' argument")
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  if (!is.null(csize)) {
    assert_that(csize > 0,
                msg = "csize must be positive")
  }

  mask.idx <- which(mask != 0)
  grid <- index_to_coord(mask, mask.idx)

  if (is.null(cvol)) {
    kres <- kmeans(grid, centers=csize, iter.max=500)
    cvol <- ClusteredNeuroVol(mask, clusters=kres$cluster)
  }

  index_list <- as.list(cvol@cluster_map)
  csize <- num_clusters(cvol)

  sp <- space(mask)
  f <- function(i) {
    ind <- index_list[[as.character(i)]]
    ROIVol(sp, index_to_grid(sp,ind), data=rep(1, length(ind)))
  }


  deflist::deflist(f, csize)

}

#' Cluster-centroid searchlight over cluster time-series
#'
#' @description
#' Iterate over clusters by their centroids and, for each seed cluster, return the
#' time-series of the K nearest clusters (or those within a radius). This enables
#' searchlight analysis at the cluster level rather than individual voxels.
#'
#' @param x A `ClusteredNeuroVec` object or a `NeuroVec` plus `cvol`
#' @param cvol A `ClusteredNeuroVol` (required if `x` is a `NeuroVec`)
#' @param k Integer, number of nearest clusters including the seed (default: 10).
#'   Will be capped at the total number of clusters if specified value exceeds it.
#' @param radius Numeric distance in mm. If given, use all clusters within this radius
#'   instead of k-nearest neighbors. Cannot be used together with \code{k}.
#' @param label Optional character label for the returned windows
#'
#' @return A list of \code{ROIVec} objects, one per cluster, where each ROIVec contains:
#' \describe{
#'   \item{values}{A T×N matrix where T is the number of timepoints and N is the number
#'     of neighboring clusters (including the seed itself)}
#'   \item{coords}{The centroid coordinates of the neighboring clusters}
#' }
#' The seed cluster's time-series is always the first column in each ROIVec.
#'
#' @details
#' The function creates a searchlight around each cluster's centroid, selecting either:
#' \itemize{
#'   \item The k nearest clusters (when \code{k} is specified)
#'   \item All clusters within a given radius (when \code{radius} is specified)
#' }
#' 
#' This is particularly useful for cluster-level connectivity analyses or when
#' working with parcellated data where voxel-level searchlights would be redundant.
#'
#' @seealso 
#' \code{\link{ClusteredNeuroVec}} for creating clustered neuroimaging vectors,
#' \code{\link{searchlight}} for voxel-level searchlight analysis,
#' \code{\link{ROIVec}} for the structure of returned windows
#'
#' @importFrom stats dist
#' @export
#' @examples
#' # Create synthetic 4D data (8x8x8 volume, 10 timepoints)
#' sp4 <- NeuroSpace(c(8,8,8,10), c(1,1,1))
#' arr <- array(rnorm(8*8*8*10), dim=c(8,8,8,10))
#' vec <- NeuroVec(arr, sp4)
#' 
#' # Create a mask covering most of the volume
#' sp3 <- NeuroSpace(c(8,8,8), c(1,1,1))
#' mask_arr <- array(FALSE, dim=c(8,8,8))
#' mask_arr[2:7, 2:7, 2:7] <- TRUE
#' mask <- LogicalNeuroVol(mask_arr, sp3)
#' 
#' # Assign voxels to 10 clusters
#' n_voxels <- sum(mask_arr)
#' clusters <- sample(1:10, n_voxels, replace=TRUE)
#' cvol <- ClusteredNeuroVol(mask, clusters)
#' 
#' # Create clustered representation
#' cv <- ClusteredNeuroVec(vec, cvol)
#' 
#' # Get cluster searchlight with 3 nearest neighbors
#' windows <- cluster_searchlight_series(cv, k = 3)
#' length(windows)  # 10 windows (one per cluster)
#' 
#' # Check first window
#' roi1 <- windows[[1]]
#' dim(values(roi1))  # 10 x 3 (timepoints x neighbors)
#' 
#' # Use radius-based neighborhoods (5mm radius)
#' windows_radius <- cluster_searchlight_series(cv, radius = 5)
#' # Each window may have different number of neighbors
cluster_searchlight_series <- function(x, cvol = NULL, k = 10L, radius = NULL, label = "") {
  if (inherits(x, "NeuroVec")) {
    stopifnot(!is.null(cvol), inherits(cvol, "ClusteredNeuroVol"))
    x <- ClusteredNeuroVec(x, cvol, label = label)
  } else {
    stopifnot(inherits(x, "ClusteredNeuroVec"))
    cvol <- x@cvol
  }
  
  # Get cluster centroids
  ctr <- centroids(cvol)  # K x 3 matrix
  K <- nrow(ctr)
  
  # Default to k=10 if neither k nor radius specified
  if (is.null(k) && is.null(radius)) {
    k <- min(10L, K)
  }
  
  # Get time-series matrix and add column names
  TS <- x@ts  # T x K
  colnames(TS) <- if (!is.null(cvol@label_map) && length(cvol@label_map) == ncol(TS)) {
    names(cvol@label_map)
  } else {
    paste0("Cluster_", seq_len(ncol(TS)))
  }
  
  # Calculate pairwise distances between centroids
  dmat <- as.matrix(stats::dist(ctr))  # K x K
  
  # Function to create one ROIVec window
  make_one <- function(seed) {
    # Select neighbors based on radius or k-NN
    neigh <- if (!is.null(radius)) {
      which(dmat[seed, ] <= radius)
    } else {
      head(order(dmat[seed, ]), k)
    }
    
    # ROIVec expects: coords = N x 3, data = T x N
    coords <- ctr[neigh, , drop = FALSE]  # N x 3
    mat <- TS[, neigh, drop = FALSE]      # T x N
    
    # Create ROIVec using constructor function
    sp4 <- space(x)
    ROIVec(sp4, coords, mat)
  }
  
  # Create list of ROIVec windows
  out <- lapply(seq_len(K), make_one)
  names(out) <- colnames(TS)
  
  # Return as regular list (deflist can be added later if needed)
  out
}

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
random_searchlight <- function(mask, radius) {
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  assert_that(radius > 0,
              msg = "radius must be positive")

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
    search <- spherical_roi(mask, center_coord, radius, nonzero=TRUE)
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
    idx <- lookup[vox]
    keep_mask <- remaining[idx]

    # If no voxels kept, remove at least the center_idx
    # to ensure progress
    if (!any(keep_mask)) {
      remaining[center_idx] <- FALSE
      remain_indices <- remain_indices[remaining[remain_indices]]
      next
    }

    kept_vox <- vox[keep_mask, , drop=FALSE]

    # Row position of the center voxel inside the kept ROI coordinates
    center_row <- which(idx[keep_mask] == center_idx)
    if (length(center_row) == 0) {
      center_row <- NA_integer_  # fallback; should not happen
    }

    search2 <- new("ROIVolWindow",
                   rep(1, nrow(kept_vox)),
                   space=space(mask),
                   coords=kept_vox,
                   center_index=as.integer(center_idx),
                   parent_index=as.integer(search@parent_index))

    # Expose row index within the ROI coordinates for downstream consumers
    attr(search2, "center_row_index") <- as.integer(center_row[1])

    # Mark chosen voxels as used
    remaining[idx[keep_mask]] <- FALSE

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

#' Create a resampled spherical searchlight iterator
#'
#' @description
#' This function generates a spherical searchlight iterator by sampling regions
#' from within a brain mask. It creates searchlight spheres around randomly
#' selected center voxels (with replacement), allowing the same center and
#' surrounding voxels to appear in multiple samples. Each searchlight can also
#' draw its radius from a user-specified set of radii.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric scalar or vector specifying candidate radii (in voxel
#'   units) for the searchlight sphere. If a vector is supplied, a radius is
#'   sampled uniformly (with replacement) for each searchlight. All radii must
#'   be positive. Default is 8.
#' @param iter An integer specifying the total number of searchlights to sample.
#'   Default is 100.
#'
#' @return A \code{deferred_list} object containing \code{\linkS4class{ROIVolWindow}}
#'   objects, each representing a spherical searchlight region sampled from within the mask.
#'
#' @details
#' Searchlight centers are sampled with replacement, so the same center (and its
#' surrounding voxels) can be selected multiple times. When multiple radii are
#' provided, each searchlight independently samples one radius from the supplied
#' values.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate a resampled searchlight iterator with radii drawn from {4,6,8}
#' searchlights <- resampled_searchlight(mask, radius = c(4, 6, 8))
#' 
#'
#' @name resampled_searchlight
#' @aliases resampled_searchlight bootstrap_searchlight
#' @rdname resampled_searchlight
#' @export
resampled_searchlight <- function(mask, radius=8, iter=100) {
  assert_that(inherits(mask, "NeuroVol"),
              msg = "mask must be a NeuroVol object")
  assert_that(is.numeric(radius), length(radius) > 0, all(radius > 0),
              msg = "radius must be positive numeric")
  assert_that(length(iter) == 1, is.finite(iter), iter > 0,
              msg = "iter must be a positive number")

  iter <- as.integer(iter)

  mask.idx <- which(mask != 0)
  assert_that(length(mask.idx) > 0,
              msg = "mask contains no nonzero voxels to sample")

  grid <- index_to_grid(mask, mask.idx)

  # Sample centers with replacement; allows iter > nrow(grid)
  center_idx <- sample.int(nrow(grid), iter, replace = TRUE)

  # Sample radius per iteration when a vector is supplied
  radii <- if (length(radius) == 1L) rep(radius, iter) else sample(radius, iter, replace = TRUE)

  force(mask)
  f <- function(i) spherical_roi(mask, grid[center_idx[i],], radii[i], nonzero=TRUE)

  deflist::deflist(f, iter)
}

#' @rdname resampled_searchlight
#' @export
bootstrap_searchlight <- function(mask, radius=8, iter=100) {
  .Deprecated("resampled_searchlight", package = "neuroim2")
  resampled_searchlight(mask, radius=radius, iter=iter)
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
    f <- function(i) { spherical_roi(mask, grid[i,], radius, nonzero=nonzero) }
    deflist::deflist(f, nrow(grid))
  } else {
    # Use spherical_roi_set to get all ROIs at once
    result_list <- spherical_roi_set(
      bvol = mask,
      centroids = grid,
      radius = radius,
      nonzero = nonzero
    )

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

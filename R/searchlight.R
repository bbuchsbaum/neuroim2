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
#' \dontrun{
#' searchlights <- random_searchlight(mask, radius = 2)
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

    search2 <- new("ROIVolWindow",
                   rep(1, sum(keep_mask)),
                   space=space(mask),
                   coords=kept_vox,
                   center_index=as.integer(center_idx),
                   parent_index=as.integer(search@parent_index))

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

#' Create a bootstrap spherical searchlight iterator
#'
#' @description
#' This function generates a spherical searchlight iterator by sampling regions
#' from within a brain mask. It creates searchlight spheres around random center
#' voxels, allowing the same surrounding voxel to belong to multiple searchlight samples.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius of the searchlight sphere
#'   in voxel units. Default is 8.
#' @param iter An integer specifying the total number of searchlights to sample.
#'   Default is 100.
#'
#' @return A \code{deferred_list} object containing \code{\linkS4class{ROIVolWindow}}
#'   objects, each representing a spherical searchlight region sampled from within the mask.
#'
#' @details
#' Searchlight centers are sampled without replacement, but the same surrounding
#' voxel can belong to multiple searchlight samples.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Generate a bootstrap searchlight iterator with a radius of 6 voxels
#' \dontrun{
#' searchlights <- bootstrap_searchlight(mask, radius = 6)
#' }
#'
#' @export
#' @rdname bootstrap_searchlight
bootstrap_searchlight <- function(mask, radius=8, iter=100) {
  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  sample.idx <- sample(1:nrow(grid), iter)

  force(mask)
  f <- function(i) spherical_roi(mask, grid[sample.idx[i],], radius, nonzero=TRUE)

  #dlis <- deferred_list(lapply(1:iter, function(i) f))
  deflist::deflist(f, iter)
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
#' \dontrun{
#' searchlights <- searchlight_coords(mask, radius = 6)
#' }
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
#' \dontrun{
#' searchlights <- searchlight(mask, radius = 6, eager = TRUE)
#' }
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
#' \dontrun{
#' clust_searchlight <- clustered_searchlight(mask, csize = 5)
#' }
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

  #dlis <- deferred_list(lapply(1:csize, function(i) f))
  #dlis

  deflist::deflist(f, csize)

}

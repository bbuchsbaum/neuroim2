
#searchlight_table <- function(x, mask, radius, type=c("standard", "random")) {
#
#}


#' Create a spherical random searchlight iterator
#'
#' This function generates a spherical random searchlight iterator, which can be used to
#' analyze the local neighborhood of voxels within a given radius in a brain mask.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius of the searchlight sphere in voxel units.
#'
#' @return A list of \code{\linkS4class{ROIVolWindow}} objects, each representing a spherical searchlight region.
#'
#' @examples
#' # Create a simple brain mask
#' mask <- array(TRUE, c(10, 10, 10))
#' mask[1, 1, 1] <- FALSE
#' mask <- LogicalNeuroVol(mask, NeuroSpace(c(10,10,10)))
#' # Generate random searchlight iterator with a radius of 2 voxels
#'
#' \dontrun{searchlights <- random_searchlight(mask, radius = 2)
#' }
#'
#' @export
#' @rdname random_searchlight
random_searchlight <- function(mask, radius) {
  assert_that(inherits(mask, "NeuroVol"))

  done <- array(FALSE, dim(mask))

  mask.idx <- which(mask != 0)

  grid <- index_to_grid(mask, as.numeric(mask.idx))
  hmap <- as.list(mask.idx)
  names(hmap) <- 1:length(hmap)
  hmap <- list2env(hmap)

  lookup <- array(0, dim(mask))
  lookup[mask.idx] <- 1:length(mask.idx)

  slist <- list()
  counter <- 1
  len <- length(mask.idx)
  keys <- ls(hmap,sorted=FALSE)

  while (len > 0) {
    ## select a center voxel from remaining indices
    center <- as.integer(sample(keys,1))

    ## get a searchlight surrounding the key
    search <- spherical_roi(mask, grid[center,], radius, nonzero=TRUE)


    vox <- coords(search)
    idx <- lookup[vox]
    ret <- mget(format(idx, scientific=FALSE, trim=TRUE), envir=hmap, ifnotfound=NA)
    keep <- !is.na(unlist(ret))

    search2 <- new("ROIVolWindow", rep(1,sum(keep)), space=space(mask), coords=coords(search)[keep,,drop=FALSE],
                   center_index=as.integer(1), parent_index=as.integer(search@parent_index))

    rm(list=format(idx[keep], scientific=FALSE, trim=TRUE),envir=hmap)
    slist[[counter]] <- search2
    counter <- counter+1

    keys <- ls(hmap,sorted=FALSE)
    len <- length(keys)
  }

  slist
}


#' Create a spherical searchlight iterator that samples regions from within a mask
#'
#' This function generates a spherical searchlight iterator by sampling regions from within a brain mask.
#' It creates searchlight spheres around random center voxels, allowing the same surround voxel to belong
#' to multiple searchlight samples.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask.
#' @param radius A numeric value specifying the radius of the searchlight sphere in voxel units (default is 8).
#' @param iter An integer specifying the total number of searchlights to sample (default is 100).
#'
#' @return A \code{deferred_list} object containing \code{\linkS4class{ROIVolWindow}} objects,
#'         each representing a spherical searchlight region sampled from within the mask.
#'
#' @details Searchlight centers are sampled without replacement, but the same surround voxel can belong to multiple searchlight samples.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
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
  deferred_list2(f, iter)
}

#' Create an exhaustive searchlight iterator that only returns voxel coordinates
#'
#' This function generates an exhaustive searchlight iterator that returns voxel coordinates for each searchlight
#' sphere within the provided mask. The searchlight iterator visits every non-zero voxel in the mask as a potential center voxel.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask, containing valid central voxels for the roving searchlight.
#' @param radius A numeric value specifying the radius (in mm) of the spherical searchlight.
#' @param nonzero A logical value indicating whether to include only coordinates with nonzero values in the supplied mask (default is FALSE).
#' @param cores An integer specifying the number of cores to use for parallel computation (default is 0, which uses a single core).
#'
#' @return A \code{deferred_list} object containing matrices of integer-valued voxel coordinates, each representing a searchlight region.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#'
#' # Generate an exhaustive searchlight iterator with a radius of 6 mm
#' \dontrun{ searchlights <- searchlight_coords(mask, radius = 6)
#' }
#'
#' @export
#' @rdname searchlight_coords
#' @importFrom dbscan frNN
searchlight_coords <- function(mask, radius, nonzero=FALSE, cores=0) {
  mask.idx <- which(mask != 0)

  grid <- index_to_grid(mask, mask.idx)
  cds <- index_to_coord(mask, mask.idx)

  #rad <- rflann::RadiusSearch(cds, cds, radius=radius^2,
  #                            max_neighbour=as.integer((radius+1))^3,
  #                            build="kdtree", cores=cores, checks=1)

  rad <- dbscan::frNN(cds, eps=radius, cds)

  spmask <- space(mask)

  f <- function(i) {
    #ind <- rad$indices[[i]]
    ind <- rad$id[[i]]
    grid[ind,,drop=FALSE]
  }

  #deferred_list(map(seq_along(rad$indices), ~ f))
  len <- nrow(cds)
  deferred_list2(f, len)

  #purrr::map(seq_along(rad$indices), function(i) {
  #    ind <- rad$indices[[i]]
  #    grid[ind,,drop=FALSE]
  #})
}



#' Create an exhaustive searchlight iterator that only returns voxel coordinates
#'
#' This function generates an exhaustive searchlight iterator that returns voxel coordinates for each searchlight
#' sphere within the provided mask. The searchlight iterator visits every non-zero voxel in the mask as a potential center voxel.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask, containing valid central voxels for the roving searchlight.
#' @param radius A numeric value specifying the radius (in mm) of the spherical searchlight.
#' @param eager A logical value specifying whether to eagerly compute the searchlight ROIs (default is FALSE, which uses lazy evaluation).
#' @param nonzero A logical value indicating whether to include only coordinates with nonzero values in the supplied mask (default is FALSE).
#' @param cores An integer specifying the number of cores to use for parallel computation (default is 0, which uses a single core).
#'
#' @return A \code{deferred_list} object containing matrices of integer-valued voxel coordinates, each representing a searchlight region.
#'
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#'
#' # Generate an exhaustive searchlight iterator with a radius of 6 mm
#' \dontrun{
#' searchlights <- searchlight(mask, radius = 6, eager = TRUE)
#' }
#'
#' @export
#' @rdname searchlight_coords
#' @importFrom dbscan frNN
searchlight <- function(mask, radius, eager=FALSE, nonzero=FALSE, cores=0) {
  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  if (!eager) {
    force(mask)
    force(radius)
    f <- function(i) { spherical_roi(mask, grid[i,], radius, nonzero=nonzero) }
    #deferred_list(lapply(1:nrow(grid), function(i) f))
    deferred_list2(f, nrow(grid))
  } else {
    cds <- index_to_coord(mask, mask.idx)

    ocds <- if (nonzero) {
      grid <- index_to_grid(mask, mask.idx)
      cds
    } else {
      tmp <- index_to_coord(mask, 1:prod(dim(mask)))
      grid <- index_to_grid(mask, 1:prod(dim(mask)))
      tmp
    }

    #rad <- rflann::RadiusSearch(cds, cds, radius=radius^2, max_neighbour=as.integer((radius+1))^3,
    #                            build="kdtree", cores=cores, checks=1)

    rad <- dbscan::frNN(ocds, eps=radius, cds )
    spmask <- space(mask)

    purrr::map(seq_along(rad$id), function(i) {
      #ind <- rad$indices[[i]]
      ind <- rad$id[[i]]
      search <- new("ROIVolWindow", mask[mask.idx[ind]], space=spmask, coords=grid[ind,,drop=FALSE],
                    center_index=as.integer(1), parent_index=as.integer(mask.idx[ind[1]]))
      search
    })

  }
}


#' Create a clustered searchlight iterator
#'
#' This function generates a searchlight iterator that iterates over successive spatial clusters in an image volume.
#' It allows for the exploration of spatially clustered regions within the provided mask by using either a pre-defined
#' clustered volume or performing k-means clustering to generate the clusters.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object representing the brain mask, containing valid central voxels for the roving searchlight.
#' @param cvol An optional \code{ClusteredNeuroVol} instance representing pre-defined clusters within the mask. If provided, the 'csize' parameter is ignored.
#' @param csize An optional integer specifying the number of clusters to be generated using k-means clustering (ignored if \code{cvol} is provided).
#'
#' @return A \code{deferred_list} object containing \code{ROIVol} objects, each representing a clustered region within the image volume.
#'
#' @importFrom stats kmeans
#' @examples
#' # Load an example brain mask
#' mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#'
#' # Generate a clustered searchlight iterator with 5 clusters
#' \dontrun{
#' clust_searchlight <- clustered_searchlight(mask, csize = 5)
#' }
#'
#' @rdname searchlight
#' @export
clustered_searchlight <- function(mask, cvol=NULL, csize=NULL) {
  if (is.null(csize) && is.null(cvol)) {
    stop(paste("must provide either 'cvol' or 'csize' argument"))
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

  deferred_list2(f, csize)

}

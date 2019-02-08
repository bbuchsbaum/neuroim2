
#' random_searchlight
#'
#' Create an spherical random searchlight iterator
#'
#' @inheritParams searchlight_coords
#' @export
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

    center <- as.integer(keys[1])

    search <- spherical_roi(mask, grid[center,], radius, nonzero=TRUE)
    vox <- coords(search)
    idx <- lookup[vox]
    ret <- mget(as.character(idx), envir=hmap, ifnotfound=NA)
    keep <- !is.na(unlist(ret))

    search2 <- new("ROIVolWindow", rep(1,sum(keep)), space=space(mask), coords=coords(search)[keep,,drop=FALSE],
                   center_index=as.integer(1), parent_index=as.integer(search@parent_index))

    rm(list=as.character(idx[keep]),envir=hmap)
    slist[[counter]] <- search2
    counter <- counter+1

    keys <- ls(hmap,sorted=FALSE)
    len <- length(keys)
  }

  slist
}


#' bootstrap_searchlight
#'
#' Create a spherical searchlight iterator that samples regions from within a mask.
#'
#' @inheritParams searchlight_coords
#' @param iter the total number of searchlights to sample (default is 100).
#' @export
#' @rdname searchlight
#' @details searchlight centers are sampled without replacement, but the same surround voxel can belong to multiple searchlight samples.
#' @examples
#'
#' mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#' slight <- bootstrap_searchlight(mask, 6)
#'
bootstrap_searchlight <- function(mask, radius=8, iter=100) {
  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  sample.idx <- sample(1:nrow(grid), iter)

  force(mask)
  f <- function(i) spherical_roi(mask, grid[sample.idx[i],], radius, nonzero=TRUE)

  #dlis <- deferred_list(lapply(1:iter, function(i) f))
  deferred_list2(f, iter)
}

#' searchlight_coords
#'
#' Create an exhaustive searchlight iterator that only returns voxel coordinates
#'
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight
#' @param nonzero only include nonzero coordinates
#' @return a list ofmatrices containing of integer-valued voxel coordinates
#' @rdname searchlight
#' @importFrom rflann RadiusSearch
#' @export
searchlight_coords <- function(mask, radius, nonzero=FALSE) {
  mask.idx <- which(mask != 0)

  grid <- index_to_grid(mask, mask.idx)
  cds <- index_to_coord(mask, mask.idx)

  rad <- rflann::RadiusSearch(cds, cds, radius=radius^2, max_neighbour=as.integer((radius+1))^3, build="kdtree", cores=0, checks=1)

  spmask <- space(mask)

  f <- function(i) {
    ind <- rad$indices[[i]]
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



#' Create an exhaustive searchlight iterator
#'
#' Generates an \code{ROIVol} around each non-zero voxel in a 3D image mask
#'
#' @inheritParams searchlight_coords
#' @param eager if TRUE, then all searchlight coordinates set are generated up front. This is faster but requires more memory to store all coordinates.
#' @return a list of \code{ROIVolWindow} objects
#' @rdname searchlight
#' @importFrom rflann RadiusSearch
#' @export
searchlight <- function(mask, radius, eager=FALSE, nonzero=FALSE) {
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
    rad <- rflann::RadiusSearch(cds, cds, radius=radius^2, max_neighbour=as.integer((radius+1))^3, build="kdtree", cores=0, checks=1)

    spmask <- space(mask)
    purrr::map(seq_along(rad$indices), function(i) {
      ind <- rad$indices[[i]]
      search <- new("ROIVolWindow", mask[mask.idx[ind]], space=spmask, coords=grid[ind,,drop=FALSE],
                    center_index=as.integer(1), parent_index=as.integer(mask.idx[ind[1]]))
      search
    })

  }
}


#' Create a clustered searchlight iterator
#'
#' @inheritParams searchlight_coords
#' @param cvol a \code{ClusteredNeuroVol} instance
#' @param csize the number of clusters (ignored if \code{cvol} is provided)
#' @return an \code{iter} class
#' @importFrom stats kmeans
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

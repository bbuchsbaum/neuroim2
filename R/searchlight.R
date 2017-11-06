
#' Create an spherical random searchlight iterator

#' @param mask an volumetric image mask of type \code{\linkS4class{NeuroVol}}
#'       containing valid searchlight voxel set.
#' @param radius width in mm of spherical searchlight
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
    search2 <- ROIVol(space(mask), coords=coords(search)[keep,,drop=FALSE], data=rep(1,sum(keep)))
    attr(search2, "center") <- grid[center,]
    attr(search2, "center.index") <- mask.idx[center]
    attr(search2, "indices") <- grid_to_index(mask, vox)
    attr(search2, "length") <- nrow(vox)

    rm(list=as.character(idx[keep]),envir=hmap)
    slist[[counter]] <- search2
    counter <- counter+1

    keys <- ls(hmap,sorted=FALSE)
    len <- length(keys)
  }

  slist
}

#'
#' #' Create an spherical random searchlight iterator
#' #'
#' #' @param mask an volumetric image mask of type \code{\linkS4class{NeuroVol}}
#' #'        containing valid searchlight voxel set.
#' #' @param radius width in mm of spherical searchlight
#' #' @export
#' random_searchlight <- function(mask, radius) {
#'   assert_that(inherits(mask, "NeuroVol"))
#'
#'   done <- array(FALSE, dim(mask))
#'
#'   mask.idx <- which(mask != 0)
#'
#'   grid <- index_to_grid(mask, as.numeric(mask.idx))
#'
#'   prog <- function() { sum(done)/length(mask.idx) }
#'
#'   nextEl <- function() {
#'     if (!all(done[mask.idx])) {
#'       center <- .resample(which(!done[mask.idx]), 1)
#'       search <- spherical_roi(mask, grid[center,], radius, nonzero=TRUE)
#'       vox <- coords(search)
#'       keep <- !done[vox]
#'       vox <- vox[keep,,drop=FALSE]
#'       search2 <- ROIVol(space(mask), coords=vox, data=rep(1,sum(keep)))
#'       done[vox] <<- TRUE
#'       attr(search2, "center") <- grid[center,]
#'       attr(search2, "center.index") <- mask.idx[center]
#'       attr(search2, "indices") <- grid_to_index(mask, vox)
#'       attr(search2, "length") <- nrow(vox)
#'       search2
#'
#'     } else {
#'       stop('StopIteration')
#'     }
#'   }
#'   obj <- list(nextElem=nextEl, progress=prog)
#'   class(obj) <- c("random_searchlight", "searchlight", 'abstractiter', 'iter')
#'   obj
#' }



#' Create a spherical searchlight iterator that samples regions from within a mask.
#'
#' searchlight centers are sampled without replacement, but the same surround voxel can belong to multiple searchlight samples.
#'
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight (can be a vector which is randomly sampled)
#' @param iter the total number of searchlights to sample (default is 100).
#' @export
bootstrap_searchlight <- function(mask, radius=8, iter=100) {
  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  sample.idx <- sample(1:nrow(grid), iter)

  f <- function(i) {
    search <- spherical_roi(mask, grid[sample.idx[i],], radius, nonzero=TRUE)
    attr(search, "center") <- grid[i,]
    attr(search, "center.index") <- mask.idx[i]
    attr(search, "length") <- nrow(search@coords)
    search
  }

  dlis <- deferred_list(lapply(1:iter, function(i) f))
}


#' Create an exhaustive searchlight iterator
#'
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight
#' @return an \code{iter} class
#' @importFrom rflann RadiusSearch
#' @export
searchlight <- function(mask, radius) {
  mask.idx <- which(mask != 0)
  grid <- index_to_grid(mask, mask.idx)

  f <- function(i) {
    search <- spherical_roi(mask, grid[i,], radius, nonzero=TRUE)
    attr(search, "center") <- grid[i,]
    attr(search, "center.index") <- mask.idx[i]
    attr(search, "length") <- nrow(search@coords)
    search
  }

  deferred_list(lapply(1:nrow(grid), function(i) f))
}


#' Create a clustered searchlight iterator
#'
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param cvol a \code{ClusteredNeuroVol} instance
#' @param csize the number of clusters (ignored if \code{cvol} is provided)
#' @return an \code{iter} class
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
    ret <- ROIVol(sp, index_to_grid(sp,ind), data=rep(1, length(ind)))
    attr(ret, "mask_indices") <- mask.idx[ind]
    attr(ret, "indices") <- ind
    attr(ret, "length") <- length(ind)
    ret
  }

  dlis <- deferred_list(lapply(1:csize, function(i) f))
  dlis

}

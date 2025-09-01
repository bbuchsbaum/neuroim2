#' ClusteredNeuroVec: Cluster-aware 4D neuroimaging data
#'
#' @description
#' `ClusteredNeuroVec` creates a 4D array-like object where voxels are grouped into clusters,
#' with one time-series per cluster. All voxels within a cluster share the same time-series,
#' making it ideal for parcellated analyses (e.g., Schaefer-Yeo parcellations).
#'
#' @param x Either a `NeuroVec` object to be reduced by clusters, or a pre-computed 
#'   numeric matrix of cluster time-series (T x K, where T=time points, K=clusters)
#' @param cvol A `ClusteredNeuroVol` object defining the cluster assignments
#' @param FUN Reduction function to aggregate voxels within clusters (default: mean).
#'   Common choices include \code{mean}, \code{median}, or custom functions.
#' @param weights Optional numeric vector of per-voxel weights for weighted aggregation.
#'   Must have length equal to the number of non-zero voxels in the mask.
#' @param label Optional character label for the object (default: "")
#'
#' @return A \code{ClusteredNeuroVec} object containing:
#' \describe{
#'   \item{cvol}{The input ClusteredNeuroVol defining cluster structure}
#'   \item{ts}{A T×K matrix of cluster time-series (T=timepoints, K=clusters)}
#'   \item{cl_map}{Integer vector mapping linear voxel indices to cluster IDs}
#'   \item{label}{Character label for the object}
#' }
#'
#' @details
#' This class implements array-like 4D access while storing data efficiently as a T×K matrix
#' instead of the full voxel×time representation. Each cluster's time-series is computed by
#' applying the aggregation function (\code{FUN}) to all voxels within that cluster.
#' 
#' The object supports standard NeuroVec operations:
#' \itemize{
#'   \item Indexing: \code{x[,,,t]} to extract 3D volumes at time t
#'   \item Series extraction: \code{series(x, i, j, k)} for time-series at voxel (i,j,k)
#'   \item Matrix conversion: \code{as.matrix(x)} to get the T×K cluster matrix
#' }
#' 
#' Single-voxel clusters are handled efficiently without aggregation overhead.
#'
#' @seealso 
#' \code{\link{ClusteredNeuroVol}} for creating cluster assignments,
#' \code{\link{cluster_searchlight_series}} for cluster-based searchlight analysis,
#' \code{\link{series}} for extracting time-series
#'
#' @export
#' @examples
#' # Create synthetic 4D data (10x10x10 volume, 20 timepoints)
#' sp4 <- NeuroSpace(c(10,10,10,20), c(1,1,1))
#' arr <- array(rnorm(10*10*10*20), dim=c(10,10,10,20))
#' vec <- NeuroVec(arr, sp4)
#' 
#' # Create a mask covering the central region
#' sp3 <- NeuroSpace(c(10,10,10), c(1,1,1))
#' mask_arr <- array(FALSE, dim=c(10,10,10))
#' mask_arr[3:8, 3:8, 3:8] <- TRUE
#' mask <- LogicalNeuroVol(mask_arr, sp3)
#' 
#' # Assign voxels to 5 random clusters
#' n_voxels <- sum(mask_arr)
#' clusters <- sample(1:5, n_voxels, replace=TRUE)
#' cvol <- ClusteredNeuroVol(mask, clusters)
#' 
#' # Create clustered representation
#' cv <- ClusteredNeuroVec(vec, cvol)
#' 
#' # Access like a regular NeuroVec
#' vol_t1 <- cv[,,,1]  # 3D volume at time 1
#' ts <- series(cv, 5, 5, 5)  # time-series at voxel (5,5,5)
#' 
#' # Get cluster time-series matrix
#' cluster_matrix <- as.matrix(cv)  # T x K matrix
#' dim(cluster_matrix)  # 20 x 5
ClusteredNeuroVec <- function(x, cvol, FUN = mean, weights = NULL, label = "") {
  stopifnot(inherits(cvol, "ClusteredNeuroVol"))
  
  # Get 3D space from mask
  sp3 <- space(cvol)
  
  if (inherits(x, "NeuroVec")) {
    # Build factor for split_reduce
    nsp <- prod(dim(sp3))
    cl_map <- integer(nsp)
    nz <- which(values(cvol@mask) != 0L)
    cl_map[nz] <- as.integer(cvol@clusters)
    
    # Create factor, excluding zeros
    fac <- factor(ifelse(cl_map == 0L, NA_integer_, cl_map))
    
    # Special handling for single-voxel clusters
    # Get unique cluster IDs
    unique_clusters <- sort(unique(cvol@clusters))
    K <- length(unique_clusters)
    T <- dim(x)[4]
    
    # Initialize result matrix
    ts <- matrix(NA_real_, T, K)
    
    # Process each cluster
    for (k in seq_along(unique_clusters)) {
      cluster_id <- unique_clusters[k]
      voxel_indices <- which(cl_map == cluster_id)
      
      if (length(voxel_indices) == 1) {
        # Single voxel - extract directly
        ts[, k] <- series(x, voxel_indices)
      } else if (length(voxel_indices) > 1) {
        # Multiple voxels - use FUN to aggregate
        voxel_series <- series(x, voxel_indices, drop = FALSE)
        ts[, k] <- apply(voxel_series, 1, FUN)
      }
    }
    
  } else if (is.matrix(x)) {
    # Pre-computed matrix provided
    stopifnot(is.numeric(x))
    ts <- x
    
    # Build cluster map
    nsp <- prod(dim(sp3))
    cl_map <- integer(nsp)
    nz <- which(values(cvol@mask) != 0L)
    cl_map[nz] <- as.integer(cvol@clusters)
    
  } else {
    stop("x must be a NeuroVec or numeric matrix")
  }
  
  # Create a 4D space with time dimension
  # Use the 3D space from cvol and add the time dimension
  dims_3d <- dim(sp3)
  dims_4d <- c(dims_3d, nrow(ts))
  
  # Create a 4D NeuroSpace with the same spatial properties as sp3
  space_4d <- NeuroSpace(dims_4d, 
                        spacing = spacing(sp3),
                        origin = origin(sp3),
                        axes = sp3@axes,
                        trans = sp3@trans)
  
  new("ClusteredNeuroVec",
      cvol = cvol,
      ts = ts,
      cl_map = as.integer(cl_map),
      label = as.character(label),
      space = space_4d)
}

#' @rdname series-methods
#' @export
setMethod("series", 
          signature(x="ClusteredNeuroVec", i="numeric"),
          function(x, i, j, k, ...) {
            sp3 <- space(x@cvol)
            dims <- dim(sp3)
            
            # Convert grid coordinates to linear index
            idx <- ((k-1) * dims[1] * dims[2]) + ((j-1) * dims[1]) + i
            
            # Get cluster ID(s)
            cid <- x@cl_map[idx]
            
            if (any(cid == 0L)) {
              # Outside mask - return NA time-series
              out <- matrix(NA_real_, nrow(x@ts), length(cid))
              if (any(cid > 0L)) {
                out[, cid > 0L] <- x@ts[, cid[cid > 0L], drop = FALSE]
              }
              return(drop(out))
            }
            
            x@ts[, cid, drop = TRUE]
          })

#' @rdname space-methods
#' @export
setMethod("space", "ClusteredNeuroVec", function(x) {
  # Add time dimension to 3D space
  sp3 <- space(x@cvol)
  NeuroSpace(c(dim(sp3), nrow(x@ts)), spacing = spacing(sp3))
})

#' @rdname dim-methods
#' @export
setMethod("dim", "ClusteredNeuroVec", function(x) {
  c(dim(space(x@cvol)), nrow(x@ts))
})

#' @rdname ndim-methods
#' @export
setMethod("ndim", "ClusteredNeuroVec", function(x) 4L)

#' @rdname length-methods
#' @export
setMethod("length", signature(x = "ClusteredNeuroVec"), function(x) {
  nrow(x@ts)
})

#' Get Labels from ClusteredNeuroVec
#' 
#' @param object A ClusteredNeuroVec object
#' @rdname labels-methods
#' @export
setMethod("labels", "ClusteredNeuroVec", function(object) object@label)

#' @rdname num_clusters-methods
#' @export
setMethod("num_clusters", "ClusteredNeuroVec", function(x) ncol(x@ts))

#' @rdname as.matrix-methods
#' @param by For ClusteredNeuroVec: controls the conversion target.
#' Defaults to "cluster" to return a T×K matrix of cluster time-series.
#' "voxel" is reserved for future use.
#' @export
setMethod("as.matrix", signature(x = "ClusteredNeuroVec"),
          function(x, by = c("cluster", "voxel")) {
            by <- match.arg(by)
            if (by == "cluster") {
              # Return T x K matrix with cluster labels as column names
              m <- x@ts
              if (!is.null(x@cvol@label_map) && length(x@cvol@label_map) == ncol(m)) {
                colnames(m) <- names(x@cvol@label_map)
              } else {
                colnames(m) <- paste0("Cluster_", seq_len(ncol(m)))
              }
              m
            } else {
              stop("by='voxel' not yet implemented for ClusteredNeuroVec")
            }
          })

#' @rdname centroids-methods
#' @export
setMethod("centroids", "ClusteredNeuroVec", 
          function(x, type = c("center_of_mass", "medoid")) {
            centroids(x@cvol, type = type)
          })

#' @rdname values-methods
#' @export
setMethod("values", "ClusteredNeuroVec", function(x) {
  x@ts
})

#' @rdname extractor4d
#' @export
setMethod("[", signature(x = "ClusteredNeuroVec", i = "missing", j = "missing"),
          function(x, i, j, k, m, ..., drop = TRUE) {
            # Extract 3D volume at time m by broadcasting cluster values
            stopifnot(length(m) == 1, m >= 1, m <= nrow(x@ts))
            
            sp3 <- space(x@cvol)
            vol_data <- array(NA_real_, dim = dim(sp3))
            
            # Fill in cluster values
            nz <- which(x@cl_map > 0)
            vol_data[nz] <- x@ts[m, x@cl_map[nz]]
            
            NeuroVol(vol_data, sp3)
          })

#' @rdname extractor4d
#' @export
setMethod("[", signature(x = "ClusteredNeuroVec", i = "numeric", j = "numeric"),
          function(x, i, j, k, m, ..., drop = TRUE) {
            # Extract specific voxel time-series value(s)
            stopifnot(length(i) == length(j), length(i) == length(k))
            
            sp3 <- space(x@cvol)
            dims <- dim(sp3)
            
            # Convert grid coordinates to linear indices
            idx <- ((k-1) * dims[1] * dims[2]) + ((j-1) * dims[1]) + i
            
            # Get cluster IDs
            cids <- x@cl_map[idx]
            
            # Extract values for each time point
            result <- matrix(NA_real_, length(m), length(idx))
            for (t_idx in seq_along(m)) {
              t <- m[t_idx]
              for (v_idx in seq_along(idx)) {
                if (cids[v_idx] > 0) {
                  result[t_idx, v_idx] <- x@ts[t, cids[v_idx]]
                }
              }
            }
            
            if (drop) {
              drop(result)
            } else {
              result
            }
          })

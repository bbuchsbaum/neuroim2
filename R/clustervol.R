#' ClusteredNeuroVol
#'
#' Construct a \code{\linkS4class{ClusteredNeuroVol}} instance
#'
#' @param mask an instance of class \code{\linkS4class{LogicalNeuroVol}}
#' @param clusters a vector of clusters ids with length equal to number of nonzero
#' voxels in mask \code{mask}
#' @param label_map an optional \code{list} that maps from cluster id to a cluster
#' label, e.g. (1 -> "FFA", 2 -> "PPA")
#' @param label an optional \code{character} string used to label of the volume
#' @return \code{\linkS4class{ClusteredNeuroVol}} instance
#'
#' @details
#'
#' The use case of \code{ClusteredNeuroVol} is to store volumetric data that has been clustered into discrete sets of voxels,
#' each of which has an associated id. For example, this class can be used to represent parcellated neuroimaging volumes.
#'
#' @export ClusteredNeuroVol
#' @examples
#'
#' # Create a simple space and volume
#' space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
#' vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
#' vol <- NeuroVol(vol_data, space)
#'
#' # Create a binary mask (e.g., values > 0)
#' mask_data <- vol_data > 0
#' mask_vol <- LogicalNeuroVol(mask_data, space)
#'
#' # Get coordinates of masked voxels
#' mask_idx <- which(mask_data)
#' coords <- index_to_coord(mask_vol, mask_idx)
#'
#' # Cluster the coordinates into 10 groups
#' set.seed(123)  # for reproducibility
#' kmeans_result <- kmeans(coords, centers = 10)
#'
#' # Create the clustered volume
#' clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)
#'
#' # Print information about the clusters
#' print(clustered_vol)
#' @rdname ClusteredNeuroVol-class
ClusteredNeuroVol <- function(mask, clusters, label_map=NULL, label="") {
  mask <- as(mask, "LogicalNeuroVol")
  space <- space(mask)
  ids <- sort(unique(clusters))

  stopifnot(length(clusters) == sum(mask))

  if (length(ids) == 1) {
    warning("clustered volume only contains 1 partition")
  }

  if (is.null(label_map)) {
    labs <- paste("Clus_", ids, sep="")
    label_map <- as.list(ids)
    names(label_map) <- labs
  } else {
    stopifnot(length(label_map) == length(ids))
    stopifnot(all(unlist(label_map) %in% ids))
  }


  clus_idx <- which(mask == TRUE)
  #cds <- index_to_coords(mask, clus_idx)

  clus_split <- split(clus_idx, clusters)
  clus_names <- names(clus_split)
  cluster_map <- new.env()

  for (i in 1:length(clus_split)) {
    cluster_map[[clus_names[[i]]]] <- clus_split[[clus_names[[i]]]]
  }

  sv <- Matrix::sparseVector(x=clusters, i=clus_idx, length=prod(dim(space)))
  #svol <- SparseNeuroVol(clusters, space(mask), indices=which(mask != 0))
  new("ClusteredNeuroVol", data=sv, mask=mask, clusters=as.integer(clusters),
      label_map=label_map, cluster_map=cluster_map, space=space)
}

#' Convert a ClusteredNeuroVol Object to a DenseNeuroVol Object
#' @name as-ClusteredNeuroVol-DenseNeuroVol
#' @aliases coerce,ClusteredNeuroVol,DenseNeuroVol-method
#' @title Convert ClusteredNeuroVol to DenseNeuroVol
#' @description This method converts a ClusteredNeuroVol into an equivalent DenseNeuroVol object.
#' @param from A \code{\linkS4class{ClusteredNeuroVol}} object to be converted
#' @return A \code{\linkS4class{DenseNeuroVol}} object
#' @examples
#'
#' # Create a clustered volume
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' clusters <- rep(1:5, length.out=sum(mask))
#' cvol <- ClusteredNeuroVol(mask, clusters)
#'
#' # Convert to DenseNeuroVol
#' dvol <- as(cvol, "DenseNeuroVol")
#'
#' @seealso \code{\linkS4class{ClusteredNeuroVol}}, \code{\linkS4class{DenseNeuroVol}}
setAs(from="ClusteredNeuroVol", to="DenseNeuroVol",
    def=function(from) {
        if (length(from@clusters) == 0) {
            stop("Cannot coerce empty ClusteredNeuroVol to DenseNeuroVol")
        }
        data <- from@clusters
        indices <- which(from@mask == TRUE)
        if (length(indices) != length(data)) {
            stop("Mismatch between mask indices and cluster data length")
        }
        DenseNeuroVol(data, space(from), indices=indices)
    })


#' Convert ClusteredNeuroVol to a base array
#'
#' Ensures that clustered volumes dispatch through the `as.array` S4 generic and
#' return dense arrays of cluster labels aligned to the underlying space.
#'
#' @param x A `ClusteredNeuroVol` instance.
#' @param ... Additional arguments (currently ignored).
#' @return A dense array of cluster ids.
#' @export
setMethod("as.array", signature(x="ClusteredNeuroVol"), function(x, ...) {
  as(x, "array")
})

#' @export
#' @rdname show-methods
setMethod(f="show", signature=signature("ClusteredNeuroVol"),
    def=function(object) {
      sp <- space(object)
      dims <- dim(object)
      spacing <- sp@spacing
      origin <- sp@origin
      n_clusters <- num_clusters(object)
      n_voxels <- sum(object@mask)

      # Header
      cat("\n")
      cat(crayon::bold(crayon::blue("ClusteredNeuroVol")), "\n")
      cat(crayon::silver(paste0(rep("=", 60), collapse="")), "\n\n")

      # Type and basic info
      cat(crayon::yellow(" > Type:          "),
          crayon::white("Clustered Volume"), "\n")

      # Dimensions section
      cat(crayon::yellow(" > Dimensions:    "),
          crayon::white(paste(dims, collapse=" x ")), "\n")

      # Space information
      cat(crayon::yellow(" > Spacing:       "),
          crayon::white(paste(signif(spacing, 3), collapse=" x ")), " ",
          crayon::silver("mm"), "\n")

      cat(crayon::yellow(" > Origin:        "),
          crayon::white(paste(signif(origin, 3), collapse=" x ")), " ",
          crayon::silver("mm"), "\n")

      # Anatomical orientation
      cat(crayon::yellow(" > Orientation:   "),
          crayon::white(paste(sp@axes@i@axis, sp@axes@j@axis, sp@axes@k@axis,
                            collapse=" x ")), "\n")

      # Cluster information
      cat("\n", crayon::bold(crayon::cyan("Cluster Information")), "\n")
      cat(crayon::silver(paste0(rep("-", 40), collapse="")), "\n")

      cat(crayon::yellow(" > Total Clusters:"),
          crayon::white(sprintf("%d", n_clusters)), "\n")

      cat(crayon::yellow(" > Active Voxels: "),
          crayon::white(sprintf("%d", n_voxels)), " ",
          crayon::silver(sprintf("(%.1f%% of volume)",
                               100 * n_voxels/prod(dims[1:3]))), "\n")

      # Label information if available
      if (!is.null(names(object@label_map))) {
        cat("\n", crayon::bold(crayon::magenta("Region Labels")), "\n")
        cat(crayon::silver(paste0(rep("-", 40), collapse="")), "\n")

        # Show first few labels with ellipsis if there are more
        n_show <- min(5, length(object@label_map))
        label_sample <- head(names(object@label_map), n_show)

        for (i in seq_along(label_sample)) {
          cat(crayon::yellow(" > "),
              crayon::white(sprintf("%-20s", label_sample[i])),
              crayon::silver(sprintf("[%d]", unlist(object@label_map[label_sample[i]]))),
              "\n")
        }

        if (length(object@label_map) > n_show) {
          cat(crayon::silver(sprintf("   ... and %d more regions\n",
                                   length(object@label_map) - n_show)))
        }
      }

      cat("\n")
    }
)


#' @export
#' @param type the type of center of mass: one of "center_of_mass" or "medoid"
#' @details For `type = "center_of_mass"`, returns arithmetic mean coordinates; for `"medoid"`, returns the most central point.
#' @return A matrix of coordinates where each row represents the centroid of a cluster.
#' @rdname centroids-methods
setMethod(f="centroids", signature=signature(x="ClusteredNeuroVol"),
          def = function(x, type=c("center_of_mass", "medoid")) {
            type <- match.arg(type)
            if (type == "center_of_mass") {
              do.call(rbind, split_clusters(x@mask, x) %>% map(~ centroid(.)) )
            } else {
              if (!requireNamespace("Gmedian", quietly = TRUE)) {
                stop("Package \"Gmedian\" needed for this function to work. Please install it.",
                     call. = FALSE)
              }
              do.call(rbind, split_clusters(x@mask, x) %>% map(~ Gmedian::Gmedian(coords(.)) ))
            }
          })

#' Split Clusters for NeuroVec Objects
#'
#' @description
#' These methods split a NeuroVec object into multiple ROIVec objects based on cluster assignments.
#'
#' @param x A NeuroVec object to be split.
#' @param clusters Either a ClusteredNeuroVol object or an integer vector of cluster assignments.
#'
#' @return A deflist (lazy-loading list) of ROIVec objects, where each element corresponds to a cluster.
#'
#' @details
#' There are two methods for splitting clusters:
#' \itemize{
#'   \item Using a ClusteredNeuroVol object: This method uses the pre-defined clusters in the ClusteredNeuroVol object.
#'   \item Using an integer vector: This method allows for custom cluster assignments.
#' }
#'
#' methods return a deflist, which is a lazy-loading list of ROIVec objects.
#'
#' @seealso
#' \code{\link{NeuroVec-class}}, \code{\link{ClusteredNeuroVol-class}}, \code{\link{ROIVec-class}}
#'
#' @examples
#' \donttest{
#'   # Create a synthetic 3D volume and its NeuroSpace
#'   space <- NeuroSpace(c(10, 10, 10,4))
#'   vol_data <- array(rnorm(10 * 10 * 10 * 4), dim = c(10, 10, 10,4))
#'   neuro_vec <- NeuroVec(vol_data, space)
#'
#'   # Create a binary mask (e.g., select voxels with values > 0)
#'   mask_data <- as.logical(neuro_vec[[1]] > .5)
#'   mask_vol <- LogicalNeuroVol(mask_data, NeuroSpace(c(10, 10, 10)))
#'
#'   # Extract indices and coordinates for the masked voxels
#'   mask_idx <- which(mask_data)
#'   coords <- index_to_coord(mask_vol, mask_idx)
#'
#'   # Perform k-means clustering on the coordinates (e.g., 3 clusters)
#'   set.seed(123)  # for reproducibility
#'   k_res <- kmeans(coords, centers = 3)
#'
#'   # Create a ClusteredNeuroVol using the mask and k-means cluster assignments
#'   clustered_vol <- ClusteredNeuroVol(mask_vol, k_res$cluster)
#'
#'   # Split the NeuroVec by clusters using the ClusteredNeuroVol method
#'   split_result_clust <- split_clusters(neuro_vec, clustered_vol)
#'
#'   # Calculate and print the mean value for each cluster
#'   means_clust <- sapply(split_result_clust, function(x) mean(values(x)))
#'   print(means_clust)
#'
#'   # Alternatively, create an integer vector of cluster assignments:
#'   cluster_assignments <- numeric(prod(dim(space)[1:3]))
#'   cluster_assignments[mask_idx] <- k_res$cluster
#'   split_result_int <- split_clusters(neuro_vec, as.integer(cluster_assignments))
#'
#'   # Verify that both splitting methods yield the same cluster means
#'   means_int <- sapply(split_result_int, function(x) mean(values(x)))
#'   print(all.equal(sort(means_clust), sort(means_int)))
#' }
#'
#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="ClusteredNeuroVol"),
          def = function(x, clusters) {
            # reuse integer path to ensure ordering matches integer split
            assert_that(prod(dim(x)[1:3]) == length(clusters@mask))
            clus_vec <- integer(length(clusters@mask))
            active_idx <- which(clusters@mask > 0)
            clus_vec[active_idx] <- clusters@clusters
            split_clusters(x, clus_vec)
          })

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="integer"),
          def = function(x, clusters) {
            assert_that(length(clusters) == prod(dim(x)[1:3]))

            unique_clusters <- sort(unique(clusters[clusters != 0]))

            f <- function(i) {
              id <- unique_clusters[i]
              idx <- which(clusters == id)
              ROIVec(space(x), index_to_grid(x, idx), x[idx])
            }

            deflist::deflist(f, length(unique_clusters))
          })



#' split_clusters
#'
#' @export
#' @rdname split_clusters-methods
#' @examples
#'
#' # Create a simple example space and data
#' space <- NeuroSpace(c(10, 10, 10,4))
#' data <- array(rnorm(1000*4), dim = c(10, 10, 10,4))
#' vec <- NeuroVec(data, space)
#'
#' # Create a mask for clustering (e.g., values > 0)
#' mask <- vec[,,,1] > 0
#' mask_vol <- LogicalNeuroVol(as.array(mask), NeuroSpace(c(10, 10, 10)))
#'
#' # Get coordinates of masked voxels for clustering
#' mask_idx <- which(mask)
#' coords <- index_to_coord(mask_vol, mask_idx)
#'
#' # Perform clustering on the coordinates (3 clusters for example)
#' set.seed(123) # for reproducibility
#' kmeans_result <- kmeans(coords, centers = 3)
#'
#' # Create a ClusteredNeuroVol
#' clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)
#'
#' # Split the NeuroVec by clusters
#' split_result <- split_clusters(vec, clustered_vol)
#'
#' # Calculate mean value for each cluster
#' cluster_means <- sapply(split_result, function(x) mean(values(x)))
#' print(cluster_means)
#'
#' # Alternative: using integer cluster assignments
#' cluster_indices <- numeric(prod(dim(space)[1:3]))
#' cluster_indices[mask_idx] <- kmeans_result$cluster
#' split_result2 <- split_clusters(vec, as.integer(cluster_indices))
#'
#' # Verify both methods give same results
#' cluster_means2 <- sapply(split_result2, function(x) mean(values(x)))
#' print(all.equal(sort(cluster_means), sort(cluster_means2)))
setMethod(f="split_clusters", signature=signature(x="NeuroVol", clusters="ClusteredNeuroVol"),
          def = function(x,clusters) {
            ids <- sort(unique(clusters@clusters))

            f <- function(i) {
              id <- ids[i]
              idx <- clusters@cluster_map[[as.character(id)]]
              if (is.null(idx)) {
                # fallback: derive from clusters slot
                idx <- which(clusters@clusters == id)
                idx <- which(clusters@mask@.Data)[idx]
              }
              dat <- linear_access(x, as.numeric(idx))
              ROIVol(space(x), index_to_grid(x,as.numeric(idx)), dat)
            }

            dlis <- deflist::deflist(f, length(ids))

          })


#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVol", clusters="integer"),
          def = function(x,clusters) {
            assert_that(length(clusters) == prod(dim(x)[1:3]))
            ind <- which(clusters > 0)
            clusters <- clusters[ind]
            clist <- split(ind, clusters)

            f <- function(i) {
              idx <- clist[[i]]
              dat <- linear_access(x, as.numeric(idx))
              ROIVol(space(x), index_to_grid(x,as.numeric(idx)), dat)
            }


            dlis <- deflist::deflist(f, length(clist))
          })

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVol", clusters="numeric"),
          def = function(x,clusters) {
            callGeneric(x,as.integer(clusters))
          })

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="ClusteredNeuroVol", clusters="missing"),
          def = function(x,clusters) {
            ids <- ls(envir = x@cluster_map)
            f <- function(i) {
              id <- ids[i]
              idx <- x@cluster_map[[id]]
              coords <- index_to_grid(x@mask, as.numeric(idx))
              ROIVol(space(x@mask), coords, rep(as.integer(id), length(idx)))
            }
            deflist::deflist(f, length(ids))
          })



#' Number of Clusters
#'
#' This function returns the number of clusters in a ClusteredNeuroVol object.
#'
#' @param x A ClusteredNeuroVol object.
#'
#' @return An integer representing the number of clusters in the input object.
#'
#' @export
#' @rdname num_clusters-methods
setMethod(f="num_clusters", signature=signature(x="ClusteredNeuroVol"),
          def=function(x) {
            length(x@cluster_map)
          })


#' @rdname as.dense-methods
#' @export
#' @return A \code{\linkS4class{NeuroVol}} object representing the dense version of the clustered volume.
setMethod("as.dense", signature(x="ClusteredNeuroVol"),
          function(x) {
            NeuroVol(as.vector(x@data), space(x@mask))
          })

#' @rdname mask-methods
#' @export
setMethod("mask", "ClusteredNeuroVol",
          function(x) {
            x@mask
          })


# ---- sub_clusters for ClusteredNeuroVol ------------------------------------

#' @rdname sub_clusters-methods
#' @export
setMethod("sub_clusters", signature(x = "ClusteredNeuroVol", ids = "integer"),
          function(x, ids, ...) {
            valid_ids <- sort(unique(x@clusters))
            bad <- setdiff(ids, valid_ids)
            if (length(bad) > 0) {
              stop("cluster ids not found: ", paste(bad, collapse = ", "))
            }
            keep <- x@clusters %in% ids
            mask_indices <- which(x@mask@.Data)
            new_mask_data <- array(FALSE, dim(x@mask))
            new_mask_data[mask_indices[keep]] <- TRUE
            new_mask <- LogicalNeuroVol(new_mask_data, space(x@mask))
            new_clusters <- x@clusters[keep]
            new_label_map <- x@label_map[vapply(x@label_map,
                                                function(v) v %in% ids, logical(1))]
            ClusteredNeuroVol(new_mask, new_clusters, label_map = new_label_map)
          })

#' @rdname sub_clusters-methods
#' @export
setMethod("sub_clusters", signature(x = "ClusteredNeuroVol", ids = "numeric"),
          function(x, ids, ...) sub_clusters(x, as.integer(ids)))

#' @rdname sub_clusters-methods
#' @export
setMethod("sub_clusters", signature(x = "ClusteredNeuroVol", ids = "character"),
          function(x, ids, ...) {
            lm <- x@label_map
            found <- ids %in% names(lm)
            if (!all(found)) {
              stop("cluster names not found: ", paste(ids[!found], collapse = ", "))
            }
            int_ids <- as.integer(unlist(lm[ids]))
            sub_clusters(x, int_ids)
          })

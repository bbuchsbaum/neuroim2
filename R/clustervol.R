#' ClusteredNeuroVol
#'
#' Construct a \code{\linkS4class{ClusteredNeuroVol}} instance
#' @param mask an instance of class \code{\linkS4class{LogicalNeuroVol}}
#' @param clusters a vector of clusters ids with length equal to number of nonzero voxels in mask \code{mask}
#' @param label_map an optional \code{list} that maps from cluster id to a cluster label, e.g. (1 -> "FFA", 2 -> "PPA")
#' @param label an optional \code{character} string used to label of the volume
#' @return \code{\linkS4class{ClusteredNeuroVol}} instance
#' @export ClusteredNeuroVol
#' @examples
#'
#' bspace <- NeuroSpace(c(16,16,16), spacing=c(1,1,1))
#' grid <- index_to_grid(bspace, 1:(16*16*16))
#' kres <- kmeans(grid, centers=10)
#' mask <- NeuroVol(rep(1, 16^3),bspace)
#' clusvol <- ClusteredNeuroVol(mask, kres$cluster)
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

#' @name
#'
#' conversion from ClusteredNeuroVol to LogicalNeuroVol
#' @rdname as-methods
setAs(from="ClusteredNeuroVol", to="DenseNeuroVol", def=function(from) {
  data = from@clusters
  indices <- which(from@mask == TRUE)
  DenseNeuroVol(data, space(from), indices=indices)
})


#' show a \code{ClusteredNeuroVol}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("ClusteredNeuroVol"),
          def=function(object) {
            sp <- space(object)
            cat("NeuroVol\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(signif(sp@spacing[1:(length(sp@spacing)-1)],2), " X ", collapse=" "),
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(signif(sp@origin[1:(length(sp@origin)-1)],2), " X ", collapse=" "),
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")
            cat("  Num Clusters   :", num_clusters(object))
          }
)

#' @export
#' @param type the type of center of mass: one of "center_of_mass" or "medoid"
#' @rdname centroids-methods
#' @importFrom Gmedian Gmedian
setMethod(f="centroids", signature=signature(x="ClusteredNeuroVol"),
          def = function(x, type=c("center_of_mass", "medoid")) {
            type <- match.arg(type)
            if (type == "center_of_mass") {
              do.call(rbind, split_clusters(x@mask, x) %>% map(~ centroid(.)) )
            } else {
              do.call(rbind, split_clusters(x@mask, x) %>% map(~ Gmedian(coords(., real=TRUE)) ))
            }
          })

## TODO add split_clusters for neurovec

#' split_clusters
#'
#' @export
#' @rdname split_clusters-methods
#' @examples
#'
#' ## split 'NeuroVol' with a 'ClusteredNeuroVol'
#' vol <- NeuroVol(array(runif(10*10*10),c(10,10,10)), NeuroSpace(c(10,10,10)))
#' mask <- as.logical(vol > .5)
#' mask.idx <- which(mask != 0)
#' grid <- index_to_coord(mask, mask.idx)
#' vox <- index_to_grid(mask, mask.idx)
#'
#' library(purrr)
#' ## partition coordinates into 50 clusters using 'kmeans'
#' kres <- kmeans(grid, centers=50, iter.max=500)
#' kvol <- ClusteredNeuroVol(mask, kres$cluster)
#' klis <- split_clusters(mask, kvol)
#' ret1 <- vol %>% split_clusters(kvol) %>% purrr::map_dbl(~ mean(values(.)))
#'
#' ## split NeuroVol with 'integer' vector of clusters.
#' indices <- numeric(prod(dim(mask)[1:3]))
#'
#' ## locations with a cluster value of 0 are ignored
#' indices[mask.idx] <- kres$cluster
#'
#' ret2 <- vol %>% split_clusters(as.integer(indices)) %>% purrr::map_dbl(~ mean(values(.)))
#' all(ret1 == ret1)
#'
setMethod(f="split_clusters", signature=signature(x="NeuroVol", clusters="ClusteredNeuroVol"),
          def = function(x,clusters) {
            f <- function(i) {
              idx <- clusters@cluster_map[[as.character(i)]]
              ROIVol(space(x), index_to_grid(x,as.numeric(idx)), x[idx])
            }

            dlis <- deferred_list(lapply(1:num_clusters(clusters), function(i) f))
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
              ROIVol(space(x), index_to_grid(x,as.numeric(idx)), x[idx])
            }

            dlis <- deferred_list(lapply(1:length(clist), function(i) f))
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
            callGeneric(x,as.vector(x@data))
          })


#' get number of clusters in a ClusteredNeuroVol
#' @rdname num_clusters-methods
#' @export
setMethod(f="num_clusters", signature=signature(x="ClusteredNeuroVol"),
          def=function(x) {
            length(x@cluster_map)
          })


#' @rdname as.dense-methods
#' @export
setMethod("as.dense", signature(x="ClusteredNeuroVol"),
          function(x) {

          })




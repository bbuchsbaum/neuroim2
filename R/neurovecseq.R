#' @include all_class.R
NULL
#' @include all_generic.R
NULL
#' @importFrom methods new
#' @importFrom assertthat assert_that
#' @importFrom purrr map map_lgl map_dbl reduce
#' @importFrom deflist deflist

#' @name NeuroVecSeq
#'
#' @title NeuroVecSeq: A Container for Sequential NeuroVec Objects
#'
#' @description
#' The NeuroVecSeq class provides a container for managing a sequence of NeuroVec objects,
#' particularly useful for handling time series or multi-session neuroimaging data where
#' each segment may have different lengths.
#'
#' @details
#' NeuroVecSeq objects store:
#' \itemize{
#'   \item A list of NeuroVec objects, each potentially with different time dimensions
#'   \item The lengths of each constituent NeuroVec
#'   \item A combined NeuroSpace object representing the total space
#' }
#'
#' The class provides methods for:
#' \itemize{
#'   \item Accessing individual time points across all vectors
#'   \item Extracting subsequences
#'   \item Computing statistics across the sequence
#'   \item Linear access to the underlying data
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{[[}{Extract a single volume at a specified time point}
#'   \item{length}{Get the total number of time points}
#'   \item{sub_vector}{Extract a subsequence of volumes}
#'   \item{linear_access}{Access data linearly across all vectors}
#' }
#'
#' @examples
#' # Create some example NeuroVec objects
#' v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
#'                space = NeuroSpace(dim = c(5, 5, 5, 2)))
#' v2 <- NeuroVec(array(1, c(5, 5, 5, 4)),
#'                space = NeuroSpace(dim = c(5, 5, 5, 4)))
#' v3 <- NeuroVec(array(2, c(5, 5, 5, 6)),
#'                space = NeuroSpace(dim = c(5, 5, 5, 6)))
#'
#' # Combine them into a sequence
#' vs <- NeuroVecSeq(v1, v2, v3)
#'
#' # Access properties
#' length(vs)  # Total time points
#' vs[[5]]     # Get the 5th volume
#'
#' # Extract a subsequence
#' sub_seq <- sub_vector(vs, 1:5)
#'
#'
#' @seealso
#' \code{\linkS4class{NeuroVec}} for the base vector class,
#' \code{\linkS4class{NeuroSpace}} for spatial information
#'
#' @export
NULL

#' Create a NeuroVecSeq Instance
#'
#' @description
#' Constructs a NeuroVecSeq object to represent a variable-length sequence of NeuroVec objects.
#' This is particularly useful for managing time series data where different segments may have
#' different lengths.
#'
#' @param ... One or more instances of type \code{\linkS4class{NeuroVec}}.
#'
#' @return A NeuroVecSeq object containing:
#' \itemize{
#'   \item The provided NeuroVec objects
#'   \item Associated space information
#'   \item Length information for each vector
#' }
#'
#' @details
#' The function performs several validations:
#' \itemize{
#'   \item Ensures all inputs are NeuroVec objects
#'   \item Verifies spatial compatibility
#'   \item Combines spatial information appropriately
#' }
#'
#' @examples
#' # Create sample vectors
#' v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
#'                space = NeuroSpace(dim = c(5, 5, 5, 2)))
#' v2 <- NeuroVec(array(0, c(5, 5, 5, 4)),
#'                space = NeuroSpace(dim = c(5, 5, 5, 4)))
#'
#' # Combine into sequence
#' vs <- NeuroVecSeq(v1, v2)
#' print(vs)
#'
#'
#' @export
NeuroVecSeq <- function(...) {
  vecs <- list(...)
  assert_that(all(map_lgl(vecs, ~ inherits(., "NeuroVec"))),
              msg = "All inputs must be NeuroVec objects")

  # Validate spatial compatibility
  sp <- space(vecs[[1]])
  base_dims <- dim(sp)[1:3]
  assert_that(
    all(map_lgl(vecs, function(x) {
      all(dim(space(x))[1:3] == base_dims)
    })),
    msg = "All NeuroVec objects must have the same spatial dimensions"
  )

  # Calculate lengths and create combined space
  lens <- map_dbl(vecs, function(x) dim(x)[4])
  sp <- add_dim(drop_dim(sp), sum(lens))

  new("NeuroVecSeq", space=sp, vecs=vecs, lens=lens)
}

#' Get length of NeuroVecSeq object
#'
#' @description Returns the total number of time points across all vectors in the sequence
#'
#' @param x A NeuroVecSeq object
#' @return Integer length (total number of time points)
#'
#' @export
#' @rdname length-methods
setMethod("length", signature=c("NeuroVecSeq"),
          def=function(x) {
            sum(map_dbl(x@vecs, ~ length(.)))
          })

#' Extract Element from NeuroVecSeq
#'
#' @description
#' Extracts a single volume from a NeuroVecSeq object at the specified time point.
#'
#' @param x A NeuroVecSeq object
#' @param i Numeric index specifying the time point to extract
#'
#' @return A NeuroVol object representing the extracted volume
#'
#' @rdname NeuroVecSeq-methods
#' @export
setMethod(f="[[", signature=signature(x="NeuroVecSeq", i="numeric"),
          def = function(x, i) {
            assert_that(length(i) == 1 && i > 0 && i <= dim(x)[4],
                        msg = "Index must be a single positive integer within bounds")

            # Calculate cumulative offsets
            offsets <- cumsum(c(1, x@lens))[1:(length(x@lens))]
            vnum <- i - offsets
            vnum[vnum < 0] <- Inf
            bucket <- which.min(vnum)
            bucket_elnum <- vnum[bucket] + 1

            x@vecs[[bucket]][[bucket_elnum]]
          })

#' Extract Subsequence from NeuroVecSeq
#'
#' @description
#' Extracts a subsequence of volumes from a NeuroVecSeq object.
#'
#' @param x A NeuroVecSeq object
#' @param i Numeric vector of indices specifying the time points to extract
#'
#' @return A NeuroVecSeq object containing the extracted subsequence
#'
#' @rdname sub_vector-methods
#' @export
setMethod(f="sub_vector", signature=signature(x="NeuroVecSeq", i="numeric"),
          def=function(x, i) {
            assertthat::assert_that(max(i) <= dim(x)[4],
                                  msg = "Indices must be within bounds")

            lens <- sapply(x@vecs, function(v) dim(v)[4])
            offset <- c(0, cumsum(lens)) + 1

            svecs <- list()
            for (j in seq_along(x@vecs)) {
              idx <- i[i >= offset[j] & i < offset[j+1]]
              if (length(idx) > 0) {
                vidx <- idx - offset[j] + 1
                svecs[[length(svecs) + 1]] <- sub_vector(x@vecs[[j]], vidx)
              }
            }

            do.call(NeuroVecSeq, svecs)
          })

#' Linear Access to NeuroVecSeq Data
#'
#' @description
#' Provides linear access to the data across all vectors in the sequence.
#'
#' @param x A NeuroVecSeq object
#' @param i Numeric vector of indices for linear access
#'
#' @return Numeric vector of accessed values
#'
#' @rdname linear_access-methods
#' @export
## @import bit64
setMethod(f = "linear_access",
          signature = signature(x = "NeuroVecSeq", i = "numeric"),
          def = function(x, i) {
            assertthat::assert_that(is.numeric(i),
                                  msg = "Index must be numeric")

            # Calculate dimensions and offsets
            nels <- prod(dim(x)[1:3])
            els <- nels * x@lens
            csum <- cumsum(nels * x@lens) + 1
            cels <- c(1, csum[-length(csum)])

            # Find which vector each index belongs to
            vnum <- findInterval(i, cels)
            offsets <- i - cels[vnum] + 1

            # Split indices by vector
            soff <- split(offsets, vnum)
            sind <- split(seq_along(vnum), vnum)

            # Access each vector
            res <- purrr::map(names(soff), function(vnum) {
              x@vecs[[as.integer(vnum)]][soff[[vnum]]]
            })

            # Combine results
            out <- numeric(length(i))
            for (j in seq_along(res)) {
              out[sind[[j]]] <- res[[j]]
            }

            out
          })


#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVecSeq", i="integer"),
          function(x, i, j, k, drop=TRUE) {

            if (missing(j) && missing(k)) {
              # i is a vector of linear voxel indices

              ts_list <- lapply(x@vecs, function(v) {
                m <- series(v, i, drop=FALSE)  # This might be [voxels x time]

                # If we want [time x voxels], transpose if needed
                # (assuming length(i) = # of voxels, and dim(v)[4] = # of timepoints)
                if (nrow(m) == length(i) && ncol(m) == dim(v)[4]) {
                  m <- t(m)  # now [time x voxels]
                }
                m
              })

              # Now row-bind them, so total rows = sum of all time
              out <- do.call(rbind, ts_list)

              # If only 1 voxel and drop=TRUE, flatten to a single vector
              if (drop && length(i) == 1) {
                out <- drop(out)  # yields length=total_time
              }
              return(out)

            } else {
              # i,j,k provided => single 3D coordinate or sets of them
              # the logic is the same: just ensure each sub-vec's series()
              # ends up as [time x coordinate], then rbind
              if (length(i) == 1 && length(j) == 1 && length(k) == 1) {

                ts_list <- lapply(x@vecs, function(v) {
                  m <- series(v, i, j, k)
                  m
                })
                #out <- do.call(rbind, ts_list)
                out <- unlist(ts_list)
                if (drop) {
                  out
                } else {
                  as.matrix(out)
                }
                return(out)
              } else {
                stop("Multiple 3D voxel coordinates not yet implemented for NeuroVecSeq.")
              }
            }
          }
)

#' @export
#' @rdname series-methods
setMethod("series", signature(x="NeuroVecSeq", i="numeric"),
          function(x,i,j,k,drop=TRUE) {
            # Just convert numeric i to integer and call the integer method
            callGeneric(x, as.integer(i), j, k, drop=drop)
          })

#' @rdname series-methods
#' @param x A NeuroVecSeq object
#' @param i A matrix of voxel coordinates (n x 3)
#' @return A matrix where each column represents a voxel's time series
#' @export
setMethod("series", signature(x="NeuroVecSeq", i="matrix"),
          def=function(x, i) {
            assertthat::assert_that(ncol(i) == 3, msg="Coordinate matrix must have 3 columns")
            # More efficient to pre-allocate and bind
            do.call(rbind, purrr::map(x@vecs, ~ series(., i)))
          })

#' @rdname series-methods
#' @param x A NeuroVecSeq object
#' @param i A matrix of ROI coordinates (n x 3)
#' @return A ROIVec object containing the time series for the specified ROI
#' @export
setMethod("series_roi", signature(x="NeuroVecSeq", i="matrix"),
          def=function(x, i) {
            assertthat::assert_that(ncol(i) == 3, msg="ROI coordinate matrix must have 3 columns")
            # Get ROIs for each vector in sequence
            rois <- purrr::map(x@vecs, ~ series_roi(., i))

            # Optimize concatenation based on number of ROIs
            if (length(rois) == 1) {
              rois[[1]]
            } else if (length(rois) == 2) {
              concat(rois[[1]], rois[[2]])
            } else {
              # For 3+ ROIs, use reduce for efficient concatenation
              purrr::reduce(rois, concat)
            }
          })


#' @export
#' @rdname vectors-methods
#' @return A deflist object where each element is a function that returns the time series
#'         for a voxel. The length of the deflist equals the total number of voxels.
setMethod("vectors", signature(x="NeuroVecSeq", subset="missing"),
          function(x) {
            # For a NeuroVec, vectors() returns a deflist of length nvox,
            # each element is a function returning the voxel time course.
            nvox <- prod(dim(x)[1:3])
            ind <- seq_len(nvox)
            f <- function(i) series(x, ind[i]) # returns a vector (the time series)
            deflist::deflist(f, nvox)
          })

#' @export
#' @rdname vectors-methods
setMethod("vectors", signature(x="NeuroVecSeq", subset="numeric"),
          function(x, subset) {
            n <- length(subset)
            f <- function(i) series(x, subset[i])
            deflist::deflist(f, n)
          })

#' @export
#' @rdname vectors-methods
setMethod("vectors", signature(x="NeuroVecSeq", subset="logical"),
          function(x, subset) {
            ind <- which(subset)
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname show-methods
setMethod("show", "NeuroVecSeq",
          def=function(object) {
            cat("\n", crayon::bold(crayon::blue("NeuroVecSeq")), " ",
                crayon::silver(paste0("(", length(object@vecs), " vectors)")), "\n", sep="")

            cat(crayon::bold("\n+= Sequence Info "), crayon::silver("---------------------------"), "\n", sep="")
            cat("| ", crayon::yellow("Length"), "        : ", length(object@vecs), "\n", sep="")
            cat("| ", crayon::yellow("Total Time"), "    : ", sum(object@lens), " points\n", sep="")

            sp <- space(object@vecs[[1]])
            cat(crayon::bold("\n+= Spatial Info "), crayon::silver("---------------------------"), "\n", sep="")
            cat("| ", crayon::yellow("Dimensions"), "    : ", paste(dim(object@vecs[[1]])[1:3], collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Spacing"), "       : ", paste(sp@spacing[1:3], collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Origin"), "        : ", paste(round(sp@origin[1:3], 2), collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Orientation"), "   : ", paste(sp@axes@i@axis, sp@axes@j@axis, sp@axes@k@axis), "\n", sep="")

            cat(crayon::bold("\n+= Vector Details "), crayon::silver("--------------------------"), "\n", sep="")
            for (i in seq_along(object@vecs)) {
              v <- object@vecs[[i]]
              vclass <- sub(".*:", "", class(v)[1])
              cat("  ", crayon::green(paste0(i, ".")), " ",
                  crayon::cyan(vclass), " ",
                  crayon::silver(paste0("(", dim(v)[4], " timepoints)")),
                  "\n", sep="")
            }
            cat("\n")
          })

#' @rdname NeuroVecSeq-methods
#' @aliases length,NeuroVecSeq-method
#'          show,NeuroVecSeq-method
#'          series,NeuroVecSeq,numeric-method
#'          vectors,NeuroVecSeq,logical-method
#'          vectors,NeuroVecSeq,missing-method
#'          vectors,NeuroVecSeq,numeric-method
#' @export

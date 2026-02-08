#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Arithmetic and Comparison Operations for Neuroimaging Objects
#'
#' @name neuro-ops
#' @description Methods for performing arithmetic and comparison operations on neuroimaging objects
NULL

#' @importFrom assertthat assert_that
#' @keywords internal
#' @noRd
checkDim <- function(e1,e2) {
  assert_that(all(dim(e1) == dim(e2)))
  assert_that(all(spacing(e1) == spacing(e2)))

}

#' Comparison Operations
#'
#' @name Compare-methods
#' @aliases Compare,DenseNeuroVol,DenseNeuroVol-method
#'          Compare,DenseNeuroVol,numeric-method
#'          Compare,numeric,DenseNeuroVol-method
#'          Compare,SparseNeuroVol,numeric-method
#'          Compare,numeric,SparseNeuroVol-method
#'          Compare,NeuroVec,NeuroVec-method
#' @description Methods for comparing neuroimaging objects.
#'   All volume comparisons return \code{\linkS4class{LogicalNeuroVol}} objects
#'   that preserve spatial metadata.
#'
#' @param e1,e2 Neuroimaging objects or numeric values.
#' @return A \code{\linkS4class{LogicalNeuroVol}} for volume comparisons.
#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="DenseNeuroVol", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1, e2)
            ret <- callGeneric(e1@.Data, e2@.Data)
            LogicalNeuroVol(ret, space(e1))
          })

#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="DenseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(e1@.Data, e2)
            LogicalNeuroVol(ret, space(e1))
          })

#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="numeric", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            ret <- callGeneric(e1, e2@.Data)
            LogicalNeuroVol(ret, space(e2))
          })

#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="SparseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(as.vector(e1@data), e2)
            LogicalNeuroVol(array(ret, dim(e1)), space(e1))
          })

#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="numeric", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            ret <- callGeneric(e1, as.vector(e2@data))
            LogicalNeuroVol(array(ret, dim(e2)), space(e2))
          })


#' Arithmetic Operations
#'
#' @name Arith-methods
#' @aliases Arith,SparseNeuroVol,SparseNeuroVol-method
#'          Arith,DenseNeuroVec,DenseNeuroVec-method
#'          Arith,SparseNeuroVol,NeuroVol-method
#'          Arith,NeuroVol,SparseNeuroVol-method
#' @description Methods for performing arithmetic operations on neuroimaging objects
#'
#' @param e1 A SparseNeuroVol object.
#' @param e2 A SparseNeuroVol object.
#' @return A SparseNeuroVol object representing the result of the arithmetic operation.
#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="SparseNeuroVol", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            ret <- callGeneric(e1@data,e2@data)
            if (is.numeric(ret)) {
              SparseNeuroVol(ret, space(e1), indices=which(ret!=0))
            } else {
              SparseNeuroVol(ret@x, space(e1), indices=ret@i)
            }

          })


#'
#' This function performs arithmetic operations on two ROIVol objects.
#'
#' @param e1 An ROIVol object.
#' @param e2 An ROIVol object.
#'
#' @return An ROIVol object resulting from the arithmetic operation.
#'
#' @export
setMethod(f="Arith", signature=signature(e1="ROIVol", e2="ROIVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)

            idx1 <- grid_to_index(e1@space, e1@coords)
            idx2 <- grid_to_index(e2@space, e2@coords)

            indices <- sort(union(idx1, idx2))
            v1 <- numeric(length(indices))
            v2 <- numeric(length(indices))
            v1[indices %in% idx1] <- e1@.Data
            v2[indices %in% idx2] <- e2@.Data
            res <- callGeneric(v1,v2)

            ROIVol(space(e1), data=res,
                   coords = index_to_grid(space(e1), indices))

          })

#' Perform arithmetic operations between two ROIVol objects
#'
#' This method performs arithmetic operations between two ROIVol objects (\code{e1} and \code{e2}) using a generic arithmetic function.
#' The dimensions of both objects are checked for compatibility before performing the operation.
#'
#' @param e1 An ROIVol object to be used in the arithmetic operation.
#' @param e2 An ROIVol object to be used in the arithmetic operation.
#' @return An ROIVol object containing the result of the arithmetic operation between \code{e1} and \code{e2}.
#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="DenseNeuroVol", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            ret <- callGeneric(e1@.Data,e2@.Data)
            bv <- DenseNeuroVol(ret, space(e1))
            bv
          })



# #' @export
# #' @name Arith
# #' @rdname Arith-methods
# #' @param e1 A SparseNeuroVec object.
# #' @param e2 A SparseNeuroVec object.
# #' @return A SparseNeuroVec object representing the result of the arithmetic operation.
# #' @description Perform an arithmetic operation between two SparseNeuroVec objects.
# #' The input SparseNeuroVec objects must have the same dimensions and NeuroSpace objects.
# #' The method computes the union of the masks and performs the arithmetic operation
# #' on the non-zero values. The result is returned as a new SparseNeuroVec object.
# setMethod(f="Arith", signature=signature(e1="SparseNeuroVec", e2="SparseNeuroVec"),
#           def=function(e1, e2) {
#             checkDim(e1, e2)
#             if (!identical(space(e1), space(e2))) {
#               stop("The NeuroSpace objects of e1 and e2 must be identical.")
#             }
#
#             mask_union <- e1@mask | e2@mask
#             indices_union <- which(mask_union)
#             data_e1 <- matrix(0, nrow(e1@data), length(indices_union))
#             data_e2 <- matrix(0, nrow(e2@data), length(indices_union))
#
#             # Fill the data matrices with their corresponding values
#             indices_e1 <- indices(e1)
#             indices_e2 <- indices(e2)
#             data_e1[, indices_union %in% indices_e1] <- e1@data[, indices_e1 %in% indices_union]
#             data_e2[, indices_union %in% indices_e2] <- e2@data[, indices_e2 %in% indices_union]
#
#             # Perform the arithmetic operation
#             result_data <- callGeneric(data_e1, data_e2)
#
#             # Create the resulting SparseNeuroVec object
#             result <- SparseNeuroVec(data=result_data, space=space(e1), mask=mask_union)
#             return(result)
#           })


#' @export
#' @rdname Arith-methods
#' @param e1 A DenseNeuroVec object.
#' @param e2 A DenseNeuroVec object.
#' @return A DenseNeuroVec object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between two DenseNeuroVec objects.
#' The input DenseNeuroVec objects must have the same dimensions and NeuroSpace objects.
#' The method computes the elementwise arithmetic operation
#' and returns a new DenseNeuroVec object.
setMethod(f="Arith", signature=signature(e1="DenseNeuroVec", e2="DenseNeuroVec"),
          def=function(e1, e2) {

            checkDim(e1,e2)
			      ret <- callGeneric(e1@.Data,e2@.Data)
            DenseNeuroVec(ret, space(e1))
          })




#' @export
#' @rdname Arith-methods
#' @param e1 A SparseNeuroVol object.
#' @param e2 A NeuroVol object.
#' @return A DenseNeuroVol object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between a SparseNeuroVol object and a NeuroVol object.
#' The input SparseNeuroVol and NeuroVol objects must have the same dimensions.
#' The method performs the arithmetic operation on the non-zero values of the SparseNeuroVol
#' and the corresponding values of the NeuroVol. The result is returned as a new DenseNeuroVol object.
setMethod(f="Arith", signature=signature(e1="SparseNeuroVol", e2="NeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            ret <- callGeneric(as.vector(e1@data), as.vector(e2@.Data))
            DenseNeuroVol(ret, space(e1))
          })




#' @export
#' @rdname Arith-methods
#' @param e1 A NeuroVol object.
#' @param e2 A SparseNeuroVol object.
#' @return A DenseNeuroVol object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between a NeuroVol object and a SparseNeuroVol object.
#' The input NeuroVol and SparseNeuroVol objects must have the same dimensions.
#' The method performs the arithmetic operation on the values of the NeuroVol and the non-zero values
#' of the SparseNeuroVol. The result is returned as a new DenseNeuroVol object.
setMethod(f="Arith", signature=signature(e1="NeuroVol", e2="SparseNeuroVol"),
            def=function(e1, e2) {
              checkDim(e1,e2)
              ret <- callGeneric(as.vector(e1@.Data), as.vector(e2@data))
              DenseNeuroVol(ret, space(e1))
      })




#' @export
#' @rdname Arith-methods
#' @param e1 A SparseNeuroVec object.
#' @param e2 A SparseNeuroVec object.
setMethod(f="Arith", signature=signature(e1="SparseNeuroVec", e2="SparseNeuroVec"),
          def=function(e1, e2) {
            checkDim(e1, e2)

            # Global voxel indices for each sparse vector
            idx1 <- indices(e1)
            idx2 <- indices(e2)
            union_idx <- sort(unique(c(idx1, idx2)))

            if (length(union_idx) == 0) {
              stop("Resulting SparseNeuroVec has no non-zero elements.")
            }

            # Allocate result-aligned matrices (rows=time, cols=union voxels)
            n_time <- nrow(e1@data)
            res1 <- matrix(0, nrow = n_time, ncol = length(union_idx))
            res2 <- matrix(0, nrow = n_time, ncol = length(union_idx))

            res1[, match(idx1, union_idx)] <- e1@data
            res2[, match(idx2, union_idx)] <- e2@data

            ret <- callGeneric(res1, res2)

            # Keep only voxels with any non-zero value across time
            keep <- which(colSums(ret != 0) > 0)
            if (length(keep) == 0) {
              stop("Resulting SparseNeuroVec has no non-zero elements.")
            }

            mask_dims <- space(e1)@dim[1:3]
            new_mask <- array(FALSE, mask_dims)
            new_mask[union_idx[keep]] <- TRUE

            SparseNeuroVec(data = ret[, keep, drop = FALSE],
                           space = space(e1),
                           mask = new_mask)
          })


#' @export
#' @rdname Arith-methods
#' @param e1 A NeuroVec object.
#' @param e2 A NeuroVec object.
#' @return A DenseNeuroVec object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between two NeuroVec objects.
#' The input NeuroVec objects must have the same dimensions.
#' The method performs the arithmetic operation on the elements of the NeuroVec objects.
#' The result is returned as a new DenseNeuroVec object.
setMethod(f="Arith", signature=signature(e1="NeuroVec", e2="NeuroVec"),
          def=function(e1, e2) {
            if (!all(dim(e1) == dim(e2))) {
              stop("cannot perform arithmetic operation on arguments with different dimensions")
            }

            D4 <- dim(e1)[4]
            vols <- list()

            for (i in 1:D4) {
              ## sub_vol(e1,i)
              vols[[i]] <- callGeneric(e1[[i]], e2[[i]])
            }

            mat <- do.call(cbind, vols)
            dspace <- add_dim(space(vols[[1]]), length(vols))
            DenseNeuroVec(mat, dspace)

          })

#' Arithmetic Operations for NeuroVec and NeuroVol
#'
#' This function performs arithmetic operations on a NeuroVec object and a
#' NeuroVol object.
#'
#' @param e1 A NeuroVec object.
#' @param e2 A NeuroVol object.
#'
#' @return A DenseNeuroVec object resulting from the arithmetic operation.
#'
#' @export
setMethod(f="Arith", signature=signature(e1="NeuroVec", e2="NeuroVol"),
		  def=function(e1, e2) {
			  if (!all(dim(e1)[1:3] == dim(e2))) {
				  stop("cannot perform arithmetic operation on arguments with different
               spatial dimensions")
			  }

			  D4 <- dim(e1)[4]
			  vols <- list()

			  for (i in 1:D4) {
			    ## sub_vol(e1,i)
				  vols[[i]] <- callGeneric(e1[[i]], e2)
			  }

			  mat <- do.call(cbind, vols)
			  dspace <- add_dim(space(vols[[1]]), length(vols))
			  DenseNeuroVec(mat, dspace)


		  })

#' Arithmetic Operations for NeuroVol and NeuroVec
#'
#' This function performs arithmetic operations on a NeuroVol object and a
#' NeuroVec object.
#'
#' @param e1 A NeuroVol object.
#' @param e2 A NeuroVec object.
#'
#' @return A DenseNeuroVec object resulting from the arithmetic operation.
#'
#' @export
setMethod(f="Arith", signature=signature(e1="NeuroVol", e2="NeuroVec"),
          def=function(e1, e2) {
            if (!all(dim(e1) == dim(e2)[1:3])) {
              stop("cannot perform arithmetic operation on arguments with
                    different spatial dimensions")
            }

            D4 <- dim(e2)[4]
            vols <- list()

            for (i in 1:D4) {
              vols[[i]] <- callGeneric(e1, e2[[i]])
            }

            mat <- do.call(cbind, vols)
            dspace <- add_dim(space(e1), D4)
            DenseNeuroVec(mat, dspace)
          })



#' Summary Methods for Neuroimaging Objects
#'
#' Methods for the \code{Summary} group generic (e.g., \code{sum}, \code{min},
#' \code{max}, \code{range}, \code{prod}, \code{any}, \code{all}) applied to
#' neuroimaging data objects.
#'
#' @name Summary-methods
#' @aliases Summary,SparseNeuroVec-method
#' @param x A neuroimaging object (SparseNeuroVec, SparseNeuroVol, or DenseNeuroVol)
#' @param ... Additional arguments passed to methods
#' @param na.rm Logical indicating whether to remove NA values before computation
#' @return The result of the summary operation
#'
#' @examples
#' # Create a simple volume
#' vol <- DenseNeuroVol(array(1:27, c(3,3,3)),
#'                      NeuroSpace(c(3L,3L,3L), c(1,1,1)))
#' sum(vol)
#' range(vol)
#'
#' @export
#' @rdname Summary-methods
setMethod(f="Summary", signature=signature(x="SparseNeuroVec"),
		def=function(x, ..., na.rm = FALSE) {
			callGeneric(x@data, ..., na.rm = na.rm)
		})

#' @export
#' @rdname Summary-methods
setMethod(f="Summary", signature=signature(x="SparseNeuroVol"),
    def=function(x, ..., na.rm = FALSE) {
      callGeneric(x@data, ..., na.rm = na.rm)
    })

#' @export
#' @rdname Summary-methods
setMethod(f="Summary", signature=signature(x="DenseNeuroVol"),
    def=function(x, ..., na.rm = FALSE) {
      callGeneric(x@.Data, ..., na.rm = na.rm)
    })

#' @export
#' @rdname Summary-methods
setMethod(f="Summary", signature=signature(x="DenseNeuroVol", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@.Data, ..., na.rm=na.rm)
    })

#' Compare two NeuroVec objects
#'
#' This method compares two NeuroVec objects (\code{e1} and \code{e2}) using a generic comparison function.
#' The dimensions of both objects are checked for compatibility before performing the comparison.
#'
#' @param e1 A NeuroVec object to be compared.
#' @param e2 A NeuroVec object to be compared.
#' @return The result of the comparison between \code{e1} and \code{e2}.
#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="NeuroVec", e2="NeuroVec"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            callGeneric(e1@.Data, e2@.Data)
          })


# ---- Scalar Arith for DenseNeuroVol ----------------------------------------

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="DenseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(e1@.Data, e2)
            DenseNeuroVol(ret, space(e1))
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="numeric", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            ret <- callGeneric(e1, e2@.Data)
            DenseNeuroVol(ret, space(e2))
          })


# ---- Scalar Arith for SparseNeuroVol ---------------------------------------

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="SparseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(as.vector(e1@data), e2)
            DenseNeuroVol(ret, space(e1))
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="numeric", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            ret <- callGeneric(e1, as.vector(e2@data))
            DenseNeuroVol(ret, space(e2))
          })


# ---- ClusteredNeuroVol Arith (warns about cluster loss) --------------------

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="ClusteredNeuroVol", e2="ClusteredNeuroVol"),
          def=function(e1, e2) {
            warning("Arithmetic on ClusteredNeuroVol: cluster structure is not preserved")
            callGeneric(as.dense(e1), as.dense(e2))
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="ClusteredNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            warning("Arithmetic on ClusteredNeuroVol: cluster structure is not preserved")
            callGeneric(as.dense(e1), e2)
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="numeric", e2="ClusteredNeuroVol"),
          def=function(e1, e2) {
            warning("Arithmetic on ClusteredNeuroVol: cluster structure is not preserved")
            callGeneric(e1, as.dense(e2))
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="ClusteredNeuroVol", e2="NeuroVol"),
          def=function(e1, e2) {
            warning("Arithmetic on ClusteredNeuroVol: cluster structure is not preserved")
            callGeneric(as.dense(e1), as.dense(e2))
          })

#' @rdname Arith-methods
#' @export
setMethod(f="Arith", signature=signature(e1="NeuroVol", e2="ClusteredNeuroVol"),
          def=function(e1, e2) {
            warning("Arithmetic on ClusteredNeuroVol: cluster structure is not preserved")
            callGeneric(as.dense(e1), as.dense(e2))
          })


# ---- Logic Operations (& and |) for NeuroVol types -------------------------

#' Logic Operations for Neuroimaging Volumes
#'
#' @name Logic-methods
#' @aliases Logic,DenseNeuroVol,DenseNeuroVol-method
#'          Logic,SparseNeuroVol,SparseNeuroVol-method
#'          Logic,SparseNeuroVol,NeuroVol-method
#'          Logic,NeuroVol,SparseNeuroVol-method
#'          Logic,NeuroVol,logical-method
#'          Logic,logical,NeuroVol-method
#' @description Methods for performing logical operations (\code{&} and
#'   \code{|}) on neuroimaging volume objects. Results are always returned as
#'   \code{\linkS4class{LogicalNeuroVol}} objects that preserve spatial metadata.
#'
#' @param e1,e2 Neuroimaging volume objects or logical values.
#' @return A \code{\linkS4class{LogicalNeuroVol}}.
#'
#' @examples
#' sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
#' v1 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
#' v2 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
#' intersection <- v1 & v2
#' union_mask  <- v1 | v2
#'
#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="DenseNeuroVol", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1, e2)
            ret <- callGeneric(e1@.Data, e2@.Data)
            LogicalNeuroVol(ret, space(e1))
          })

#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="SparseNeuroVol", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1, e2)
            ret <- callGeneric(as.vector(e1@data), as.vector(e2@data))
            LogicalNeuroVol(array(ret, dim(e1)), space(e1))
          })

#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="SparseNeuroVol", e2="NeuroVol"),
          def=function(e1, e2) {
            checkDim(e1, e2)
            ret <- callGeneric(as.vector(e1@data), as.vector(as.dense(e2)@.Data))
            LogicalNeuroVol(array(ret, dim(e1)), space(e1))
          })

#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="NeuroVol", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1, e2)
            ret <- callGeneric(as.vector(as.dense(e1)@.Data), as.vector(e2@data))
            LogicalNeuroVol(array(ret, dim(e1)), space(e1))
          })

#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="NeuroVol", e2="logical"),
          def=function(e1, e2) {
            ret <- callGeneric(as.dense(e1)@.Data, e2)
            LogicalNeuroVol(ret, space(e1))
          })

#' @rdname Logic-methods
#' @export
setMethod(f="Logic", signature=signature(e1="logical", e2="NeuroVol"),
          def=function(e1, e2) {
            ret <- callGeneric(e1, as.dense(e2)@.Data)
            LogicalNeuroVol(ret, space(e2))
          })


# ---- Logical NOT (!) for NeuroVol types ------------------------------------

#' Logical Negation for Neuroimaging Volumes
#'
#' @name not-methods
#' @aliases !,DenseNeuroVol-method !,SparseNeuroVol-method
#' @description Logical negation (\code{!}) for neuroimaging volumes. Returns a
#'   \code{\linkS4class{LogicalNeuroVol}} where non-zero voxels become
#'   \code{FALSE} and zero voxels become \code{TRUE}.
#'
#' @param x A neuroimaging volume object.
#' @return A \code{\linkS4class{LogicalNeuroVol}}.
#'
#' @examples
#' sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
#' mask <- LogicalNeuroVol(array(sample(c(TRUE, FALSE), 125, replace = TRUE),
#'                               c(5, 5, 5)), sp)
#' inv <- !mask
#'
#' @rdname not-methods
#' @export
setMethod(f="!", signature=signature(x="DenseNeuroVol"),
          def=function(x) {
            LogicalNeuroVol(!x@.Data, space(x))
          })

#' @rdname not-methods
#' @export
setMethod(f="!", signature=signature(x="SparseNeuroVol"),
          def=function(x) {
            LogicalNeuroVol(array(!as.vector(x@data), dim(x)), space(x))
          })

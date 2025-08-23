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
#' @aliases Compare,SparseNeuroVol,numeric-method
#'          Compare,numeric,SparseNeuroVol-method
#'          Compare,NeuroVec,NeuroVec-method
#' @description Methods for comparing neuroimaging objects
#'
#' @param e1 A SparseNeuroVol object containing the data to be compared.
#' @param e2 A numeric value to compare with the data of the SparseNeuroVol object.
#' @return The result of the comparison between the SparseNeuroVol object's data and the numeric value.
#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="SparseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(e1@data,e2)
          })


#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="numeric", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            callGeneric(e1@data,e2)
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
#' @return A DenseNeuroVol object representing the result of the arithmetic operation.
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
#' @param e1 A SparseNeuroVec object.
#' @param e2 A SparseNeuroVec object.
#' @return A SparseNeuroVec object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between two SparseNeuroVec objects.
#' The input SparseNeuroVec objects must have the same dimensions and NeuroSpace objects.
#' The method computes the union of the masks and performs the arithmetic operation
#' on the non-zero values. The result is returned as a new SparseNeuroVec object.
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




setMethod(f="Arith", signature=signature(e1="SparseNeuroVec", e2="SparseNeuroVec"),
          def=function(e1, e2) {
            # Ensure dimensions match
            checkDim(e1, e2)

            # Extract the fourth dimension from the space slot
            D4 <- space(e1)@dim[4]

            # Initialize lists to store results
            vols <- vector("list", D4)
            ind <- vector("list", D4)

            # Iterate over the fourth dimension
            for (i in seq_len(D4)) {
              # Perform the arithmetic operation on the i-th slice
              # Access data directly from the @data slot
              vols[[i]] <- callGeneric(e1@data[i, ], e2@data[i, ])

              # Extract non-zero indices from the result
              ind[[i]] <- which(vols[[i]] != 0)
            }

            # Combine all unique non-zero indices across all dimensions
            combined_ind <- sort(unique(unlist(ind)))

            # Handle case where there are no non-zero elements
            if (length(combined_ind) == 0) {
              stop("Resulting SparseNeuroVec has no non-zero elements.")
            }

            # Extract the non-zero data for each volume based on combined indices
            vret <- do.call(rbind, lapply(vols, function(vol) vol[combined_ind]))

            # Update the NeuroSpace object by ensuring the fourth dimension remains consistent
            dspace <- space(e1)  # Assuming space remains the same

            # Construct the new mask based on non-zero elements
            dims <- space(e1)@dim[1:3]
            new_mask <- array(FALSE, dims)
            new_mask[combined_ind] <- TRUE

            # Create the new SparseNeuroVec object
            SparseNeuroVec(data = vret,
                           space = dspace,
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



#' @export
#' @rdname Summary-methods
#' @param x A SparseNeuroVec object
#' @param ... Additional arguments passed to methods
#' @param na.rm Logical indicating whether to remove NA values before computation
#' @return the summary of the SparseNeuroVec object
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


#' Summary method for Neuroimaging objects
#' @rdname Summary-methods
#' @export
setMethod(f="Summary", signature=signature(x="DenseNeuroVol", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@.Data, ..., na.rm=na.rm)
    })

#' Arithmetic operations for SparseNeuroVec objects
#' @name Arith-methods
#' @rdname Arith-methods
#' @aliases Arith,SparseNeuroVec,SparseNeuroVec-method
setMethod(f="Arith", signature=signature(e1="SparseNeuroVec", e2="SparseNeuroVec"),
          def=function(e1, e2) {
            # Ensure dimensions match
            checkDim(e1, e2)

            # Extract the fourth dimension from the space slot
            D4 <- space(e1)@dim[4]

            # Initialize lists to store results
            vols <- vector("list", D4)
            ind <- vector("list", D4)

            # Iterate over the fourth dimension
            for (i in seq_len(D4)) {
              # Perform the arithmetic operation on the i-th slice
              # Access data directly from the @data slot
              vols[[i]] <- callGeneric(e1@data[i, ], e2@data[i, ])

              # Extract non-zero indices from the result
              ind[[i]] <- which(vols[[i]] != 0)
            }

            # Combine all unique non-zero indices across all dimensions
            combined_ind <- sort(unique(unlist(ind)))

            # Handle case where there are no non-zero elements
            if (length(combined_ind) == 0) {
              stop("Resulting SparseNeuroVec has no non-zero elements.")
            }

            # Extract the non-zero data for each volume based on combined indices
            vret <- do.call(rbind, lapply(vols, function(vol) vol[combined_ind]))

            # Update the NeuroSpace object by ensuring the fourth dimension remains consistent
            dspace <- space(e1)  # Assuming space remains the same

            # Construct the new mask based on non-zero elements
            dims <- space(e1)@dim[1:3]
            new_mask <- array(FALSE, dims)
            new_mask[combined_ind] <- TRUE

            # Create the new SparseNeuroVec object
            SparseNeuroVec(data = vret,
                           space = dspace,
                           mask = new_mask)
          })

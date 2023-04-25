#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' @importFrom assertthat assert_that
#' @keywords internal
checkDim <- function(e1,e2) {
  assert_that(all(dim(e1) == dim(e2)))
  assert_that(all(spacing(e1) == spacing(e2)))

}

#' Compare a SparseNeuroVol object with a numeric value
#'
#' This method compares the data of a SparseNeuroVol object (\code{e1}) with a numeric value (\code{e2}) using a generic comparison function.
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


#' Compare a numeric value with a SparseNeuroVol object
#'
#' This method compares a numeric value (\code{e1}) with the data of a SparseNeuroVol object (\code{e2}) using a generic comparison function.
#'
#' @param e1 A numeric value to compare with the data of the SparseNeuroVol object.
#' @param e2 A SparseNeuroVol object containing the data to be compared.
#' @return The result of the comparison between the numeric value and the SparseNeuroVol object's data.
#' @rdname Compare-methods
#' @export
setMethod(f="Compare", signature=signature(e1="numeric", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            callGeneric(e1@data,e2)
          })


#' Perform arithmetic operations between two SparseNeuroVol objects
#'
#' This method performs arithmetic operations between two SparseNeuroVol objects (\code{e1} and \code{e2}) using a generic arithmetic function.
#' The dimensions of both objects are checked for compatibility before performing the operation.
#'
#' @param e1 A SparseNeuroVol object to be used in the arithmetic operation.
#' @param e2 A SparseNeuroVol object to be used in the arithmetic operation.
#' @return A SparseNeuroVol object containing the result of the arithmetic operation between \code{e1} and \code{e2}.
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

            ROIVol(space(e1), data=res, coords = index_to_grid(space(e1), indices))

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


# setMethod(f="Compare", signature=signature(e1="NeuroVec", e2="NeuroVec"),
#           def=function(e1, e2) {
#             checkDim(e1,e2)
#
#           })


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
#' @param e1 A SparseNeuroVec object.
#' @param e2 A SparseNeuroVec object.
#' @return A SparseNeuroVec object representing the result of the arithmetic operation.
#' @description Perform an arithmetic operation between two SparseNeuroVec objects.
#' The input SparseNeuroVec objects must have the same dimensions.
#' The method performs the arithmetic operation on the non-zero values of the SparseNeuroVec objects.
#' The result is returned as a new SparseNeuroVec object.
setMethod(f="Arith", signature=signature(e1="SparseNeuroVec", e2="SparseNeuroVec"),
          def=function(e1, e2) {
				    D4 <- dim(e1)[4]
				    vols <- list()
				    ind <- list()

				    for (i in 1:D4) {
					    vols[[i]] <- callGeneric(e1[[i]], e2[[i]])
					    ind[[i]] <- vols[[i]]@data@i
				    }

				    ind <- sort(unique(unlist(ind)))
				    vret <- do.call(rbind, lapply(vols, function(vol) as.numeric(vol[ind])))

				    dspace <- add_dim(space(vols[[1]]), length(vols))
				    m <- logical(prod(dim(space(vols[[1]]))))
				    m[ind] <- TRUE
				    mask <- LogicalNeuroVol(m, space(vols[[1]]))
				    SparseNeuroVec(vret, dspace, mask=mask)

          }
)

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
#' This function performs arithmetic operations on a NeuroVec object and a NeuroVol object.
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
				  stop("cannot perform arithmetic operation on arguments with different spatial dimensions")
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


#' Summary of SparseNeuroVec
#'
#' This function computes a summary of a SparseNeuroVec object.
#'
#' @param x A SparseNeuroVec object.
#'
#' @return A summary of the input SparseNeuroVec object.
#'
#' @export
setMethod(f="Summary", signature=signature(x="SparseNeuroVec"),
		def=function(x) {
			callGeneric(x@data)
		})



#' Summary of SparseNeuroVol
#'
#' This function computes a summary of a SparseNeuroVol object.
#'
#' @param x A SparseNeuroVol object.
#' @param ... Additional arguments (currently ignored).
#' @param na.rm A logical value indicating whether NA values should be removed before computation.
#'
#' @return A summary of the input SparseNeuroVol object.
#'
#' @export
setMethod(f="Summary", signature=signature(x="SparseNeuroVol", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@data)
    })


#setMethod("sum", signature()


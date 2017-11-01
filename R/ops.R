#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' @importFrom assertthat assert_that
checkDim <- function(e1,e2) {
  assert_that(all(dim(e1) == dim(e2)))
  assert_that(all(spacing(e1) == spacing(e2)))

}

setMethod(f="Arith", signature=signature(e1="SparseNeuroVol", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            res <- callGeneric(e1@data,e2@data)
            new("SparseNeuroVol", data=res, source=e1@source, space=space(e1))

          })


setMethod(f="Arith", signature=signature(e1="ROIVol", e2="ROIVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)

            idx1 <- grid_to_index(e1@space, e1@coords)
            idx2 <- grid_to_index(e2@space, e2@coords)

            indices <- sort(union(idx1, idx2))
            v1 <- numeric(length(indices))
            v2 <- numeric(length(indices))
            v1[indices %in% idx1] <- e1@data
            v2[indices %in% idx2] <- e2@data
            res <- callGeneric(v1,v2)

            ROIVol(space(e1), data=res, coords = NeuroVol(space(e1), indices))

          })





setMethod(f="Arith", signature=signature(e1="NeuroVol", e2="NeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)

            ret <- callGeneric(e1@.Data,e2@.Data)
            bv <- DenseNeuroVol(ret, space(e1))

          })


setMethod(f="Arith", signature=signature(e1="NeuroVec", e2="NeuroVec"),
          def=function(e1, e2) {

      checkDim(e1,e2)

			if (inherits(e1, "DenseNeuroVec") && inherits(e2, "DenseNeuroVec")) {
            	ret <- callGeneric(e1@.Data,e2@.Data)
            	DenseNeuroVec(ret, space(e1))
			} else {
				D4 <- dim(e1)[4]
				vols <- list()
				for (i in 1:D4) {
					vols[[i]] <- callGeneric(takeVolume(e1,i), takeVolume(e2,i))
				}

				mat <- do.call(cbind, vols)
				dspace <- addDim(space(vols[[1]]), length(vols))
				DenseNeuroVec(mat, dspace)

			}

})


 setMethod(f="Arith", signature=signature(e1="NeuroVec", e2="NeuroVol"),
		  def=function(e1, e2) {
			  if (!all(dim(e1)[1:3] == dim(e2))) {
				  stop("cannot perform arithmetic operation on arguments with different spatial dimensions")
			  }

			  D4 <- dim(e1)[4]
			  vols <- list()

			  for (i in 1:D4) {
			    ## sub_vol(e1,i)
				  vols[[i]] <- callGeneric(takeVolume(e1,i), e2)
			  }

			  mat <- do.call(cbind, vols)
			  dspace <- add_dim(space(vols[[1]]), length(vols))
			  DenseNeuroVec(mat, dspace)


		  })


setMethod(f="Arith", signature=signature(e1="NeuroVol", e2="NeuroVec"),
		def=function(e1, e2) {
			callGeneric(e2,e1)
		})


setMethod(f="Summary", signature=signature(x="SparseNeuroVec"),
		def=function(x) {
			callGeneric(x@data)
		})


setMethod(f="Summary", signature=signature(x="SparseNeuroVol", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@data)
    })


#setMethod("sum", signature()


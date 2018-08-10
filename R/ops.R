#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' @importFrom assertthat assert_that
checkDim <- function(e1,e2) {
  assert_that(all(dim(e1) == dim(e2)))
  assert_that(all(spacing(e1) == spacing(e2)))

}



setMethod(f="Compare", signature=signature(e1="SparseNeuroVol", e2="numeric"),
          def=function(e1, e2) {
            ret <- callGeneric(e1@data,e2)
          })


setMethod(f="Compare", signature=signature(e1="numeric", e2="SparseNeuroVol"),
          def=function(e1, e2) {
            callGeneric(e1@data,e2)
          })



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


setMethod(f="Arith", signature=signature(e1="DenseNeuroVol", e2="DenseNeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            ret <- callGeneric(e1@.Data,e2@.Data)
            bv <- DenseNeuroVol(ret, space(e1))

          })


setMethod(f="Arith", signature=signature(e1="DenseNeuroVec", e2="DenseNeuroVec"),
          def=function(e1, e2) {

            checkDim(e1,e2)
			      ret <- callGeneric(e1@.Data,e2@.Data)
            DenseNeuroVec(ret, space(e1))
          })

setMethod(f="Arith", signature=signature(e1="SparseNeuroVol", e2="NeuroVol"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            ret <- callGeneric(as.vector(e1@data), as.vector(e2@.Data))
            DenseNeuroVol(ret, space(e1))
          })


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



setMethod(f="Summary", signature=signature(x="SparseNeuroVec"),
		def=function(x) {
			callGeneric(x@data)
		})


setMethod(f="Summary", signature=signature(x="SparseNeuroVol", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@data)
    })


#setMethod("sum", signature()


#' @include all_class.R
{}
#' @include all_generic.R
{}


#' neuro_slice constructor
#' @param data data vector or matrix
#' @param space an instance of class BrainSpace
#' @param indices linear indices corresponding used if \code{data} is a 1D vector.
#' @export
#' @examples
#' bspace <- NeuroSpace(c(64,64), spacing=c(1,1))
#' dat <- array(rnorm(64*64), c(64,64))
#' bslice <- NeuroSlice(dat,bspace)
#' print(bslice)
NeuroSlice <- function(data, space, indices=NULL) {
	if (ndim(space) != 2) {
		stop("incorrect dimension for neuro_slice")
	}



	if (is.null(indices)) {
		if (length(dim(data)) != 2) {
		  stopifnot(length(data) == prod(dim(space)[1:2]))
			data <- matrix(data, dim(space)[1], dim(space)[2])
		}

	  stopifnot(all(dim(data) == dim(space)))
		new("NeuroSlice", .Data=data, space=space)

	} else {
		mdat <- matrix(0, dim(space))
		mdat[indices] <- data
		new("NeuroSlice", .Data=mdat, space=space)
	}
}


#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x = "NeuroSlice", coords="matrix"),
		def=function(x, coords) {
			callGeneric(x@space, coords)
})

#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x = "NeuroSlice", coords="numeric"),
          def=function(x, coords) {
            callGeneric(x@space, coords)
          })



#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x = "NeuroSlice", idx="numeric"),
		def=function(x, idx) {
		  callGeneric(x@space, idx)
})







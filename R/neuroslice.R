#' @include all_class.R
{}
#' @include all_generic.R
{}


#' Construct a NeuroSlice object
#'
#' @param data data vector or matrix
#' @param space an instance of class NeuroSpace
#' @param indices linear indices corresponding used if \code{data} is a 1D vector.
#' @export
#' @examples
#' bspace <- NeuroSpace(c(64,64), spacing=c(1,1))
#' dat <- array(rnorm(64*64), c(64,64))
#' bslice <- NeuroSlice(dat,bspace)
#' bslice
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




#' @export
#' @importFrom graphics plot
#' @rdname plot-methods
#' @examples
#'
#' dat <- matrix(rnorm(100*100), 100, 100)
#' slice <- NeuroSlice(dat, NeuroSpace(c(100,100)))
#' #plot(slice)
setMethod("plot", signature=signature(x="NeuroSlice"),
          def=function(x,cmap=gray(seq(0,1,length.out=255)), irange=range(x, na.rm=TRUE)) {
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }



            ## map intensities to colors
            imcols <- mapToColors(x, cmap, alpha=1, irange=irange, zero_col="#000000")

            {y=value=NULL}

            cds <- index_to_coord(space(x), 1:length(x))
            df1 <- data.frame(x=cds[,1], y=cds[,2], value=as.vector(imcols))
            ggplot2::ggplot(ggplot2::aes(x=x, y=y), data=df1) + ggplot2::geom_raster(ggplot2::aes(fill=value)) +
              ggplot2::scale_fill_identity() + ggplot2::xlab("") + ggplot2::ylab("")
              ggplot2::theme_bw()

          })



#' @import assertthat
#' @keywords internal
#' @importFrom grDevices col2rgb gray heat.colors
mapToColors <- function (imslice, col = heat.colors(128, alpha = 1), zero_col = "#00000000",
                         alpha = 1, irange = range(imslice), threshold = c(0, 0)) {

  assertthat::assert_that(diff(irange) >= 0)
  drange <- diff(irange)
  mcols <- (imslice - irange[1])/diff(irange) * (length(col) -1) + 1
  mcols[mcols < 1] <- 1
  mcols[mcols > length(col)] <- length(col)
  imcols <- col[mcols]

  if (!is.vector(imslice)) {
    dim(imcols) <- dim(imslice)
  }

  imcols[imslice == 0] <- zero_col

  if (diff(threshold) > 0) {
    imcols[(imslice >= threshold[1]) & (imslice <= threshold[2])] <- "#00000000"
  }

  if (alpha < 1) {
    rgbmat <- col2rgb(imcols, alpha = TRUE)
    rgbmat <- rgbmat/255
    rgbmat[4, ] <- rgbmat[4, ] * alpha

    if (is.vector(imslice)) {
      array(t(rgbmat), c(length(imslice), 4))
    } else {
      array(t(rgbmat), c(dim(imslice), 4))
    }
  }
  else {
    imcols
  }
}


#' show a \code{NeuroSlice}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("NeuroSlice"),
          def=function(object) {
            sp <- space(object)
            cat("NeuroSlice\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(signif(sp@spacing[1:(length(sp@spacing)-1)],2), " X ", collapse=" "),
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(signif(sp@origin[1:(length(sp@origin)-1)],2), " X ", collapse=" "),
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis), "\n")

          }
)








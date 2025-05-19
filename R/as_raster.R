#' Convert neuroimaging objects to raster images
#'
#' These methods convert \code{NeuroSlice} and \code{NeuroVol} objects
#' into \code{raster} images that can be displayed with
#' \code{grid::grid.raster} or base \code{plot}.
#'
#' @inheritParams plot,NeuroSlice-method
#' @param x object to convert.
#' @param zlevel slice index for \code{NeuroVol} objects. Defaults to the middle slice.
#' @param ... additional arguments passed to methods
#' @return A \code{raster} object
#' @seealso \code{\link{plot,NeuroSlice-method}}, \code{\link{plot,NeuroVol-method}}
#' @rdname as.raster-methods
#' @export
setMethod("as.raster", "NeuroSlice",
          function(x,
                   cmap=gray(seq(0,1,length.out=255)),
                   irange=range(x, na.rm=TRUE), ...){
            cols <- mapToColors(x, col=cmap, alpha=1,
                                irange=irange, zero_col="#00000000")
            grDevices::as.raster(cols)
          })

#' @rdname as.raster-methods
#' @export
setMethod("as.raster", "NeuroVol",
          function(x, zlevel=ceiling(dim(x)[3]/2),
                   cmap=gray(seq(0,1,length.out=255)),
                   irange=range(x, na.rm=TRUE), ...){
            sl <- slice(x, zlevel=zlevel, along=3)
            as.raster(sl, cmap=cmap, irange=irange, ...)
          })

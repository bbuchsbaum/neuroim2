#' @include all_class.R
{}
#' @include all_generic.R
{}


#'
#' Create an IndexLookupVol object.
#'
#' @param space A \code{\linkS4class{NeuroSpace}} object representing the 3D space of the brain image.
#' @param indices An integer vector containing the 1D indices of the voxels in the grid.
#'
#' @return An object of class \code{\linkS4class{IndexLookupVol}} representing the index lookup volume.
#'
#' @seealso \code{\linkS4class{IndexLookupVol}}
#'
#' @export
#' @rdname IndexLookupVol-class
IndexLookupVol <- function(space, indices) {
  new("IndexLookupVol", space=space, indices=indices)
}


setMethod(f="initialize", signature=signature("IndexLookupVol"),
          def=function(.Object, space, indices) {

            #print("initializing IndexLookupVol")
            .Object@space <- space
            .Object@indices <- as.integer(indices)
            nels <- prod(dim(space)[1:3])
            map <- integer(nels)
            map[indices] <- 1:length(indices)
            .Object@map <- as.integer(map)
            .Object
          })

#' indices
#'
#' @export indices
#' @rdname indices-methods
setMethod(f="indices", signature=signature(x="IndexLookupVol"),
          def=function(x) {
            x@indices
          })

#' lookup
#'
#' @export lookup
#' @rdname lookup-methods
setMethod(f="lookup", signature=signature(x="IndexLookupVol", i="numeric"),
          def=function(x,i) {
            x@map[i]
          })


#' @export space
#' @rdname space-methods
setMethod(f="space", signature=signature(x="IndexLookupVol"),
          def=function(x) {
            x@space
          })

#' coords
#'
#' @export coords
#' @param i the index in to the lookup volume
#' @rdname coords-methods
setMethod(f="coords", signature(x="IndexLookupVol"),
          def=function(x,i) {
            idx <- lookup(x,i)
            idx <- idx[idx!=0]
            if (length(idx) == 0) {
              NA
            } else {
              index_to_grid(space(x), idx)
            }

          })




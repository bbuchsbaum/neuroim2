#' Create an IndexLookupVol Object
#'
#' @description
#' This function creates an IndexLookupVol object, which represents a lookup volume
#' for efficient indexing of voxels in a 3D brain image space.
#'
#' @param space A \code{\linkS4class{NeuroSpace}} object representing the 3D space of the brain image.
#' @param indices An integer vector containing the 1D indices of the voxels in the grid.
#'
#' @return An object of class \code{\linkS4class{IndexLookupVol}} representing the index lookup volume.
#'
#' @examples
#' \dontrun{
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' indices <- sample(1:262144, 10000)  # Random selection of 10000 voxels
#' ilv <- IndexLookupVol(space, indices)
#' }
#'
#' @seealso \code{\linkS4class{IndexLookupVol}}, \code{\linkS4class{NeuroSpace}}
#'
#' @export
#' @rdname IndexLookupVol-class
IndexLookupVol <- function(space, indices) {
  new("IndexLookupVol", space = space, indices = indices)
}

#' Initialize an IndexLookupVol Object
#'
#' @description
#' This method initializes an IndexLookupVol object with the given space and indices.
#'
#' @param .Object The IndexLookupVol object to initialize.
#' @param space A NeuroSpace object representing the 3D space.
#' @param indices An integer vector of indices.
#'
#' @return An initialized IndexLookupVol object.
#'
#' @keywords internal
setMethod(f = "initialize", signature = signature("IndexLookupVol"),
          def = function(.Object, space, indices) {
            .Object@space <- space
            .Object@indices <- as.integer(indices)
            nels <- prod(dim(space)[1:3])
            map <- integer(nels)
            map[indices] <- 1:length(indices)
            .Object@map <- as.integer(map)
            .Object
          })

#' Get Indices from an IndexLookupVol Object
#'
#' @description
#' This method retrieves the indices stored in an IndexLookupVol object.
#'
#' @param x An IndexLookupVol object.
#'
#' @return An integer vector of indices.
#'
#' @export
#' @rdname indices-methods
setMethod(f = "indices", signature = signature(x = "IndexLookupVol"),
          def = function(x) {
            x@indices
          })

#' Lookup Values in an IndexLookupVol Object
#'
#' @description
#' This method performs a lookup operation on an IndexLookupVol object.
#'
#' @param x An IndexLookupVol object.
#' @param i A numeric vector of indices to look up.
#'
#' @return A vector of lookup values corresponding to the input indices.
#'
#' @export
#' @rdname lookup-methods
setMethod(f = "lookup", signature = signature(x = "IndexLookupVol", i = "numeric"),
          def = function(x, i) {
            x@map[i]
          })

#' Get Space from an IndexLookupVol Object
#'
#' @description
#' This method retrieves the NeuroSpace object associated with an IndexLookupVol object.
#'
#' @param x An IndexLookupVol object.
#'
#' @return A NeuroSpace object.
#'
#' @export
#' @rdname space-methods
setMethod(f = "space", signature = signature(x = "IndexLookupVol"),
          def = function(x) {
            x@space
          })

#' Extract Coordinates from an IndexLookupVol Object
#'
#' @description
#' This method extracts the coordinates from an IndexLookupVol object based on a given index.
#'
#' @param x An IndexLookupVol object to extract coordinates from.
#' @param i The index into the lookup volume.
#'
#' @return The extracted coordinates corresponding to the provided index.
#'         If the index is not found, it returns NA.
#'
#' @examples
#' \dontrun{
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' indices <- sample(1:262144, 10000)
#' ilv <- IndexLookupVol(space, indices)
#' coords(ilv, 1000)
#' }
#'
#' @export
#' @rdname coords-methods
setMethod(f = "coords", signature(x = "IndexLookupVol"),
          def = function(x, i) {
            idx <- lookup(x, i)
            idx <- idx[idx != 0]
            if (length(idx) == 0) {
              NA
            } else {
              index_to_grid(space(x), idx)
            }
          })




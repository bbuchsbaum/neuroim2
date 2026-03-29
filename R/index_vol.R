#' Index Lookup Volume for 3D Neuroimaging Data
#'
#' @name IndexLookupVol-class
#' @description
#' The \code{IndexLookupVol} class provides efficient indexing and coordinate lookup
#' functionality for 3D neuroimaging data. It maintains a mapping between linear indices
#' and 3D coordinates, optimizing memory usage and access speed for sparse volumes.
#'
#' @section Implementation Details:
#' The class uses an integer mapping array for O(1) lookups between linear indices
#' and their corresponding positions in the sparse representation.
#'
#' @seealso \code{\link{IndexLookupVol}} for creating instances of this class
#'
NULL

#' Create an IndexLookupVol Object
#'
#' @title Create an Index Lookup Volume
#' @description
#' Creates an \code{\linkS4class{IndexLookupVol}} object, which provides efficient
#' bidirectional mapping between linear indices and 3D coordinates in a neuroimaging
#' volume. This is particularly useful for working with masked or sparse brain volumes.
#'
#' @param space A \code{\linkS4class{NeuroSpace}} object defining the 3D space dimensions,
#'   spacing, and orientation.
#' @param indices An integer vector containing the linear indices of the voxels to include
#'   in the lookup volume. These should be 1-based indices within the range of the space.
#'
#' @return An object of class \code{\linkS4class{IndexLookupVol}} containing:
#' \itemize{
#'   \item A mapping between linear indices and sparse positions
#'   \item The original space information
#'   \item The subset of included voxel indices
#' }
#'
#' @examples
#'
#' # Create a 64x64x64 space
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#'
#' # Create a lookup volume with random indices
#' indices <- sample(1:262144, 10000)  # Select 10000 random voxels
#' ilv <- IndexLookupVol(space, indices)
#'
#' # Look up coordinates for specific indices
#' coords <- coords(ilv, indices[1:10])
#'
#'
#' @seealso
#' \code{\link{coords}} for coordinate lookup,
#' \code{\link{lookup}} for index mapping,
#' \code{\linkS4class{NeuroSpace}} for space representation
#'
#' @importFrom methods new
#'
#' @export
#' @rdname IndexLookupVol-class
IndexLookupVol <- function(space, indices) {
  if (!inherits(space, "NeuroSpace")) {
    cli::cli_abort("{.arg space} must be a {.cls NeuroSpace} object.")
  }
  if (!is.numeric(indices) || !all(is.finite(indices))) {
    cli::cli_abort("{.arg indices} must be a vector of finite numeric values.")
  }

  # Convert to integer and validate range
  indices <- as.integer(indices)
  nels <- prod(dim(space)[1:3])
  if (!all(indices >= 1) || !all(indices <= nels)) {
    cli::cli_abort("{.arg indices} must be within the range of the space dimensions [1, {nels}].")
  }

  new("IndexLookupVol", space = space, indices = indices)
}

#' Initialize an IndexLookupVol Object
#'
#' @description
#' Internal method to initialize an IndexLookupVol object with proper validation
#' and mapping setup.
#'
#' @param .Object The IndexLookupVol object to initialize
#' @param space A \code{\linkS4class{NeuroSpace}} object
#' @param indices An integer vector of indices
#'
#' @return An initialized IndexLookupVol object
#'
#' @keywords internal
#' @noRd
setMethod(f = "initialize",
          signature = signature("IndexLookupVol"),
          def = function(.Object, space, indices) {
            # Ensure indices are unique and sorted
            indices <- sort(unique(as.integer(indices)))

            .Object@space <- space
            .Object@indices <- indices

            # Create reverse mapping
            nels <- prod(dim(space)[1:3])
            map <- integer(nels)
            map[indices] <- seq_along(indices)
            .Object@map <- map

            .Object
          })

#' Get Indices from an IndexLookupVol Object
#'
#' @description
#' Retrieves the vector of indices that are included in the lookup volume.
#'
#' @param x An \code{\linkS4class{IndexLookupVol}} object
#'
#'
#' @examples
#'
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' ilv <- IndexLookupVol(space, c(1:100))
#' idx <- indices(ilv)  # Get included indices
#'
#'
#' @export
#' @return the indices of the lookup volume
#' @rdname indices-methods
setMethod(f = "indices",
          signature = signature(x = "IndexLookupVol"),
          def = function(x) {
            x@indices
          })

#' Lookup Values in an IndexLookupVol Object
#'
#' @description
#' Performs a lookup operation on an IndexLookupVol object.
#'
#' @param x An \code{\linkS4class{IndexLookupVol}} object
#' @param i A numeric vector of indices to look up
#'
#' @return the values of the lookup volume
#'
#'
#' @examples
#'
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' ilv <- IndexLookupVol(space, c(1:100))
#' lookup(ilv, c(1, 2, 3))  # Look up values for indices 1, 2, and 3
#'
#'
#' @export
#' @rdname lookup-methods
setMethod(f = "lookup",
          signature = signature(x = "IndexLookupVol", i = "numeric"),
          def = function(x, i) {
            if (!is.numeric(i) || !all(is.finite(i))) {
              cli::cli_abort("{.arg i} must be a vector of finite numeric values.")
            }

            # Convert to integer and validate range
            i <- as.integer(i)
            nels <- prod(dim(x@space)[1:3])
            if (!all(i >= 1) || !all(i <= nels)) {
              cli::cli_abort("{.arg i} must be within the range of the space dimensions [1, {nels}].")
            }

            x@map[i]
          })

#' Get Space from an IndexLookupVol Object
#'
#' @description
#' Retrieves the NeuroSpace object associated with an IndexLookupVol object.
#'
#' @param x An \code{\linkS4class{IndexLookupVol}} object
#'
#'
#' @examples
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' ilv <- IndexLookupVol(space, c(1:100))
#' space(ilv)  # Get the associated NeuroSpace object
#'
#'
#' @export
#' @rdname space-methods
setMethod(f = "space",
          signature = signature(x = "IndexLookupVol"),
          def = function(x) {
            x@space
          })

#' Extract Coordinates from an IndexLookupVol Object
#'
#' @description
#' Extracts the coordinates from an IndexLookupVol object based on a given index.
#'
#' @param x An \code{\linkS4class{IndexLookupVol}} object to extract coordinates from
#' @param i The index into the lookup volume
#'
#'
#' @examples
#'
#' space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
#' ilv <- IndexLookupVol(space, c(1:100))
#' coords(ilv, 1)  # Extract coordinates for index 1
#'
#'
#' @export
#' @rdname coords-methods
setMethod(f = "coords",
          signature = signature(x = "IndexLookupVol"),
          def = function(x, i) {
            # Input validation
            if (!is.numeric(i) || !all(is.finite(i))) {
              stop("'i' must be a vector of finite numeric values")
            }

            # Convert to integer and validate range
            i <- as.integer(i)
            nels <- prod(dim(x@space)[1:3])
            if (any(i < 1) || any(i > nels)) {
              stop("'i' must be within the range of the space dimensions")
            }

            idx <- lookup(x, i)
            idx <- idx[idx != 0]
            if (length(idx) == 0) {
              NA
            } else {
              index_to_grid(x@space, idx)
            }
          })


#' @export
#' @rdname show-methods
setMethod("show", "IndexLookupVol", function(object) {
  n <- length(object@indices)
  nz <- sum(object@indices != 0)
  show_header("IndexLookupVol", paste(nz, "/", n, "mapped"))
  show_field("Total Indices", format(n, big.mark = ","))
  show_field("Non-zero", paste0(format(nz, big.mark = ","),
             " (", round(100 * nz / max(n, 1), 1), "%)"))
  show_field("Map Entries", length(object@map))
})

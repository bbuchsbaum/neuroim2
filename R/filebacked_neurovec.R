#' Create a FileBackedNeuroVec Object
#'
#' @description
#' Constructs a \code{FileBackedNeuroVec} instance, which represents a file-backed
#' neuroimaging vector object.
#'
#' @param file_name A character string specifying the name of the image file.
#' @return A new instance of class \code{FileBackedNeuroVec}.
#'
#' @details
#' This function reads the header of the specified file, creates a NeuroSpace object,
#' and initializes a FileBackedNeuroVec with the appropriate metadata.
#'
#' @examples
#' \dontrun{
#' fbvec <- FileBackedNeuroVec("path/to/image.nii")
#' }
#'
#' @seealso \code{\link{NeuroSpace}}, \code{\link{read_header}}
#'
#' @export
FileBackedNeuroVec <- function(file_name) {
  meta <- read_header(file_name)
  assert_that(length(meta@dims) == 4)
  sp <- NeuroSpace(meta@dims, meta@spacing, meta@origin,
                   meta@spatial_axes, trans(meta))

  new("FileBackedNeuroVec",
      space = sp,
      meta = meta)
}

#' Extract Sub-vector from FileBackedNeuroVec
#'
#' @param x A FileBackedNeuroVec object.
#' @param i A numeric vector specifying the indices to extract.
#'
#' @return A DenseNeuroVec object containing the extracted sub-vector.
#'
#' @details
#' This method reads the specified volumes from the file-backed data and
#' returns them as a DenseNeuroVec object.
#'
#' @export
#' @rdname sub_vector-methods
setMethod(f = "sub_vector", signature = signature(x = "FileBackedNeuroVec", i = "numeric"),
          def = function(x, i) {
            mat <- read_mapped_vols(x@meta, i)
            sp <- add_dim(drop_dim(space(x)), length(i))
            DenseNeuroVec(mat, sp)
          })

#' Convert FileBackedNeuroVec to List
#'
#' @description
#' Converts a FileBackedNeuroVec object to a list of DenseNeuroVol objects.
#'
#' @param x A FileBackedNeuroVec object.
#'
#' @return A deferred list of DenseNeuroVol objects.
#'
#' @details
#' This method creates a deferred list, where each element is a DenseNeuroVol
#' object representing a single volume from the FileBackedNeuroVec.
#'
#' @export
#' @rdname as.list-methods
setMethod(f = "as.list", signature = signature(x = "FileBackedNeuroVec"),
          def = function(x) {
            D4 <- dim(x)[4]
            f <- function(i) {
              drop(sub_vector(x, i))
            }
            deflist::deflist(f, D4)
          })


#' Convert FileBackedNeuroVec to Matrix
#'
#' @description
#' This method converts a FileBackedNeuroVec object to a matrix.
#'
#' @param from A FileBackedNeuroVec object to be converted.
#'
#' @return A matrix representation of the FileBackedNeuroVec object.
#'
#' @details
#' The resulting matrix will have rows representing time points (or the 4th dimension)
#' and columns representing voxels. The voxels are arranged in a linear order.
#'
#' @noRd
setAs(from = "FileBackedNeuroVec", to = "matrix",
      def = function(from) {
        len <- prod(from@meta@dims[1:3])
        t(series(from, seq(1, len)))
      })

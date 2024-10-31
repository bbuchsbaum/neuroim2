#' Construct a MappedNeuroVecSource Object
#'
#' @description
#' This function constructs a MappedNeuroVecSource object, which represents a source
#' for memory-mapped neuroimaging data.
#'
#' @param file_name A character string specifying the name of the image file to be memory mapped.
#'
#' @return A new instance of class \code{\linkS4class{MappedNeuroVecSource}}.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads the header of the specified file using \code{\link{read_header}}.
#' 2. Checks that the dimensions of the meta information are at least 3D.
#' 3. Creates and returns a new MappedNeuroVecSource object with the meta information.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a valid neuroimaging file named "brain.nii"
#' source <- MappedNeuroVecSource("brain.nii")
#' }
#'
#' @seealso 
#' \code{\link{MappedNeuroVec}}, \code{\link{read_header}}
#'
#' @export
#' @rdname MappedNeuroVecSource-class
MappedNeuroVecSource <- function(file_name) {
  meta <- read_header(file_name)
  stopifnot(length(dim(meta)) >= 3)
  new("MappedNeuroVecSource", meta_info = meta)
}

#' Construct a MappedNeuroVec Object
#'
#' @description
#' This function constructs a MappedNeuroVec object, which represents memory-mapped
#' neuroimaging data as a vector.
#'
#' @param file_name A character string specifying the name of the 4D image file 
#'   containing the memory-mapped data source.
#'
#' @return A new instance of class \code{\linkS4class{MappedNeuroVec}}.
#'
#' @details
#' The function performs the following steps:
#' 1. Creates a MappedNeuroVecSource object from the input file.
#' 2. Loads the data from the source using the \code{load_data} method.
#'
#' This approach allows for efficient handling of large neuroimaging datasets
#' by keeping the data on disk and only loading it into memory when needed.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a valid 4D neuroimaging file named "fmri.nii"
#' mapped_vec <- MappedNeuroVec("fmri.nii")
#' }
#'
#' @seealso 
#' \code{\link{MappedNeuroVecSource}}, \code{\link{load_data}}
#'
#' @export
#' @rdname MappedNeuroVec-class
MappedNeuroVec <- function(file_name) {
  src <- MappedNeuroVecSource(file_name)
  load_data(src)
}

#' Load Data from MappedNeuroVecSource
#'
#' @description
#' This method loads data from a MappedNeuroVecSource object and creates a MappedNeuroVec object.
#'
#' @param x A MappedNeuroVecSource object.
#'
#' @return A new instance of class \code{\linkS4class{MappedNeuroVec}}.
#'
#' @details
#' This method performs the following steps:
#' 1. Extracts meta information from the source.
#' 2. Memory maps the data file using the mmap package.
#' 3. Calculates the data offset.
#' 4. Creates a NeuroSpace object from the meta information.
#' 5. Constructs and returns a new MappedNeuroVec object.
#'
#' @keywords internal
#' @noRd
setMethod(f = "load_data", signature = c(x = "MappedNeuroVecSource"),
          def = function(x) {
            # Method body...
          })

#' Linear Access to MappedNeuroVec Data
#'
#' @description
#' This method provides linear access to the data in a MappedNeuroVec object.
#'
#' @param x A MappedNeuroVec object.
#' @param i A numeric vector of indices to access.
#'
#' @return A numeric vector of values corresponding to the requested indices.
#'
#' @details
#' This method adjusts the input indices by the offset stored in the MappedNeuroVec object
#' and then retrieves the corresponding values from the memory-mapped file.
#'
#' @keywords internal
#' @noRd
setMethod(f = "linear_access", signature = signature(x = "MappedNeuroVec", i = "numeric"),
          def = function(x, i) {
            # Method body...
          })



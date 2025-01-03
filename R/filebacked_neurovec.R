#' File-Backed Neuroimaging Vector Class
#'
#' @name FileBackedNeuroVec-class
#' @description
#' The \code{FileBackedNeuroVec} class represents a memory-efficient vector of neuroimaging data
#' that is stored on disk rather than in memory. This is particularly useful for large datasets
#' where memory constraints are a concern.
#'
#' @section Memory Management:
#' Data is read from disk on-demand, reducing memory usage compared to in-memory storage.
#' The trade-off is slightly slower access times due to disk I/O operations.
#'
NULL

#' Create a FileBackedNeuroVec Object
#'
#' @title Create a File-Backed Neuroimaging Vector
#' @description
#' Constructs a \code{\linkS4class{FileBackedNeuroVec}} instance, which represents a file-backed
#' neuroimaging vector object. This constructor provides memory-efficient access to large
#' neuroimaging datasets by keeping the data on disk until needed.
#'
#' @param file_name A character string specifying the path to the neuroimaging file.
#'   Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).
#' @param label Optional character string providing a label for the vector
#'
#' @return A new instance of class \code{\linkS4class{FileBackedNeuroVec}}.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Reads the header information from the specified file
#'   \item Validates the dimensionality (must be 4D data)
#'   \item Creates a \code{\linkS4class{NeuroSpace}} object with appropriate metadata
#'   \item Initializes the file-backed vector with minimal memory footprint
#' }
#'
#' @examples
#' \dontrun{
#' # Create a file-backed vector from a NIFTI file
#' fbvec <- FileBackedNeuroVec("path/to/fmri_data.nii")
#'
#' # Access specific volumes without loading entire dataset
#' first_vol <- sub_vector(fbvec, 1)
#' }
#'
#' @seealso 
#' \code{\link{NeuroSpace}} for spatial metadata management
#' \code{\link{read_header}} for header information extraction
#' \code{\link{sub_vector}} for data access methods
#'
#' @export
FileBackedNeuroVec <- function(file_name, label = basename(file_name)) {
  if (!is.character(file_name) || length(file_name) != 1) {
    stop("'file_name' must be a single character string")
  }
  if (!file.exists(file_name)) {
    stop("File '", file_name, "' does not exist")
  }
  
  meta <- read_header(file_name)
  if (length(meta@dims) != 4) {
    stop("Input data must be 4-dimensional (3D + time)")
  }
  
  sp <- NeuroSpace(meta@dims, meta@spacing, meta@origin,
                   meta@spatial_axes, trans(meta))

  new("FileBackedNeuroVec",
      space = sp,
      meta = meta,
      label = label)
}

#' Extract Sub-vector from FileBackedNeuroVec
#'
#' @description
#' Extracts a subset of volumes from a file-backed neuroimaging vector and returns
#' them as a dense (in-memory) vector.
#'
#' @param x A \code{\linkS4class{FileBackedNeuroVec}} object.
#' @param i A numeric vector specifying the indices of volumes to extract.
#'
#' @return A \code{\linkS4class{DenseNeuroVec}} object containing the extracted volumes.
#'
#' @details
#' This method efficiently reads only the requested volumes from disk, converting them
#' to an in-memory representation. The spatial metadata is preserved but adjusted to
#' reflect the new number of volumes.
#'
#' Memory usage is proportional to the number of volumes requested, not the size of
#' the full dataset.
#'
#' @examples
#' \dontrun{
#' fbvec <- FileBackedNeuroVec("fmri_data.nii")
#' 
#' # Extract first 10 volumes
#' subset <- sub_vector(fbvec, 1:10)
#' 
#' # Extract specific timepoints
#' volumes <- sub_vector(fbvec, c(1, 5, 10))
#' }
#'
#' @export
#' @rdname sub_vector-methods
setMethod(f = "sub_vector", 
          signature = signature(x = "FileBackedNeuroVec", i = "numeric"),
          def = function(x, i) {
            if (!is.numeric(i) || any(is.na(i))) {
              stop("Index 'i' must be a numeric vector without NA values")
            }
            if (any(i < 1) || any(i > dim(x)[4])) {
              stop("Index out of bounds")
            }
            
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
#' @examples
#' \dontrun{
#' fbvec <- FileBackedNeuroVec("fmri_data.nii")
#' 
#' # Convert to list of volumes
#' vol_list <- as.list(fbvec)
#' }
#'
#' @export
#' @rdname as.list-methods
setMethod(f = "as.list", 
          signature = signature(x = "FileBackedNeuroVec"),
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
#' @examples
#' \dontrun{
#' fbvec <- FileBackedNeuroVec("fmri_data.nii")
#' 
#' # Convert to matrix
#' mat <- as.matrix(fbvec)
#' }
#'
#' @noRd
setAs(from = "FileBackedNeuroVec", to = "matrix",
      def = function(from) {
        len <- prod(from@meta@dims[1:3])
        t(series(from, seq(1, len)))
      })

#' Linear Access Method for FileBackedNeuroVec
#'
#' @param x A FileBackedNeuroVec object
#' @param i A numeric vector of indices
#' @return A numeric vector of values
#' @export
setMethod(
  f = "linear_access",
  signature = signature(x = "FileBackedNeuroVec", i = "numeric"),
  def = function(x, i) {
    if (!is.numeric(i) || any(is.na(i))) {
      stop("Index 'i' must be a numeric vector without NA values")
    }
    if (any(i < 1) || any(i > prod(dim(x)))) {
      stop("Index out of bounds")
    }
    
    dims <- dim(x)
    spatial_nels <- prod(dims[1:3])
    
    # Calculate timepoints and spatial offsets
    timepoints <- ((i - 1) %/% spatial_nels) + 1
    spatial_offsets <- ((i - 1) %% spatial_nels) + 1
    
    # Get unique timepoints to minimize file reads
    unique_timepoints <- sort(unique(timepoints))
    
    # Read the required volumes
    mat <- read_mapped_vols(x@meta, unique_timepoints)  # Returns [time, voxels]
    
    # Create lookup table for timepoint indices
    time_lookup <- match(timepoints, unique_timepoints)
    
    # Extract values using the computed indices
    values <- numeric(length(i))
    for (idx in seq_along(i)) {
      values[idx] <- mat[time_lookup[idx], spatial_offsets[idx]]
    }
    
    values
  }
)

#' Memory-Mapped Neuroimaging Vector Class
#'
#' @name MappedNeuroVec-class
#' @description
#' The \code{MappedNeuroVec} class provides memory-efficient access to large neuroimaging
#' datasets through memory mapping. This allows processing of datasets larger than available
#' RAM by keeping data on disk and only loading requested portions into memory.
#'
#' @section Implementation Details:
#' The class uses the \code{mmap} package to establish a memory mapping between the file
#' and memory space. Key features include:
#' \itemize{
#'   \item Zero-copy access to file data
#'   \item Automatic memory management
#'   \item Support for large datasets
#'   \item Efficient random access
#' }
#'
NULL

#' Create a Memory-Mapped Source for Neuroimaging Data
#'
#' @title Create a Memory-Mapped Data Source
#' @description
#' Creates a \code{\linkS4class{MappedNeuroVecSource}} object that manages the memory
#' mapping between a neuroimaging file and memory space. This is typically used internally
#' by \code{\link{MappedNeuroVec}} but can be created directly for custom access patterns.
#'
#' @param file_name Character string specifying the path to the neuroimaging file.
#'   Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).
#'
#' @return A new \code{\linkS4class{MappedNeuroVecSource}} object containing:
#' \itemize{
#'   \item Meta information about the dataset
#'   \item File format details
#'   \item Dimensional information
#' }
#'
#' @details
#' The function performs several important checks:
#' \itemize{
#'   \item Validates file existence and permissions
#'   \item Reads and validates header information
#'   \item Ensures proper dimensionality (>= 3D)
#'   \item Verifies file format compatibility
#' }
#'
#' @examples
#' \dontrun{
#' # Create source from NIFTI file
#' source <- MappedNeuroVecSource("path/to/fmri.nii")
#' 
#' # Check dimensions
#' dim(source@meta_info)
#' 
#' # View header information
#' str(source@meta_info)
#' }
#'
#' @seealso 
#' \code{\link{MappedNeuroVec}} for the main user interface
#' \code{\link{read_header}} for header reading details
#'
#' @export
#' @rdname MappedNeuroVecSource-class
MappedNeuroVecSource <- function(file_name) {
  # Input validation
  if (!is.character(file_name) || length(file_name) != 1) {
    stop("'file_name' must be a single character string")
  }
  if (!file.exists(file_name)) {
    stop("File '", file_name, "' does not exist")
  }
  if (!file.access(file_name, mode = 4) == 0) {
    stop("File '", file_name, "' is not readable")
  }
  
  # Read and validate header
  meta <- read_header(file_name)
  if (length(dim(meta)) < 3) {
    stop("Input data must be at least 3-dimensional")
  }
  
  new("MappedNeuroVecSource", meta_info = meta)
}

#' Create a Memory-Mapped Neuroimaging Vector
#'
#' @title Create a Memory-Mapped Vector for Neuroimaging Data
#' @description
#' Creates a \code{\linkS4class{MappedNeuroVec}} object that provides efficient,
#' memory-mapped access to large neuroimaging datasets. This allows processing of
#' data larger than available RAM by keeping it on disk and only loading requested
#' portions into memory.
#'
#' @param file_name Character string specifying the path to the neuroimaging file.
#'   Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).
#' @param label Optional character string providing a label for the vector
#'
#' @return A new \code{\linkS4class{MappedNeuroVec}} object providing:
#' \itemize{
#'   \item Memory-mapped access to the data
#'   \item Spatial and temporal indexing
#'   \item Efficient data extraction
#'   \item Automatic memory management
#' }
#'
#' @details
#' The function implements several key features:
#' \itemize{
#'   \item Zero-copy access to file data
#'   \item Automatic memory management
#'   \item Support for large datasets
#'   \item Efficient random access
#'   \item Proper cleanup on object deletion
#' }
#'
#' Memory mapping is particularly useful when:
#' \itemize{
#'   \item Working with large datasets
#'   \item Only portions of data are needed at once
#'   \item Random access is required
#'   \item Multiple processes need to share data
#' }
#'
#' @examples
#' \dontrun{
#' # Create mapped vector from NIFTI file
#' mvec <- MappedNeuroVec("path/to/fmri.nii")
#' 
#' # Extract first volume
#' vol1 <- mvec[[1]]
#' 
#' # Get dimensions
#' dim(mvec)
#' 
#' # Access specific timepoint
#' timepoint <- mvec[, , , 10]
#' }
#'
#' @seealso 
#' \code{\link{MappedNeuroVecSource}} for low-level access
#' \code{\link{mmap}} for memory mapping details
#'
#' @export
#' @rdname MappedNeuroVec-class
MappedNeuroVec <- function(file_name, label = basename(file_name)) {
  # Input validation
  if (!is.character(file_name) || length(file_name) != 1) {
    stop("'file_name' must be a single character string")
  }
  if (!file.exists(file_name)) {
    stop("File '", file_name, "' does not exist")
  }
  
  # Create source and load data
  src <- MappedNeuroVecSource(file_name)
  load_data(src)
}

#' Load Data from Memory-Mapped Source
#'
#' @description
#' Internal method to load data from a MappedNeuroVecSource and create a MappedNeuroVec.
#'
#' @param x A MappedNeuroVecSource object
#'
#' @return A new MappedNeuroVec object
#'
#' @keywords internal
#' @noRd
setMethod(f = "load_data", 
          signature = "MappedNeuroVecSource",
          def = function(x) {
            # Extract meta information
            meta <- x@meta_info
            
            # Create memory mapping
            tryCatch({
              fmap <- mmap::mmap(meta@data_file, 
                                mode = .getMMapMode(meta@data_type),
                                prot = mmap::mmapFlags("PROT_READ"))
            }, error = function(e) {
              stop("Failed to create memory mapping: ", e$message)
            })
            
            # Calculate data offset
            offset <- meta@data_offset / .getDataSize(meta@data_type)
            
            # Create space object
            bspace <- NeuroSpace(dim(meta), 
                               meta@spacing, 
                               meta@origin,
                               meta@spatial_axes, 
                               trans = trans(meta))
            
            # Create mapped vector with basename as label
            new("MappedNeuroVec", 
                space = bspace, 
                filemap = fmap,
                offset = as.integer(offset),
                label = basename(meta@data_file))
          })

#' Linear Access to Memory-Mapped Data
#'
#' @description
#' Internal method providing linear access to memory-mapped data.
#'
#' @param x A MappedNeuroVec object
#' @param i Numeric vector of indices
#'
#' @return Numeric vector of values
#'
#' @keywords internal
#' @noRd
setMethod(f = "linear_access", 
          signature = signature(x = "MappedNeuroVec", i = "numeric"),
          def = function(x, i) {
            if (!is.numeric(i) || any(is.na(i))) {
              stop("Index 'i' must be a numeric vector without NA values")
            }
            
            # Calculate adjusted indices
            idx <- i + x@offset
            
            # Bounds checking
            if (any(idx < 1) || any(idx > length(x@filemap))) {
              stop("Index out of bounds")
            }
            
            # Access mapped data
            x@filemap[idx]
          })

#' Show Method for MappedNeuroVec Objects
#'
#' @description
#' Displays a formatted summary of a MappedNeuroVec object.
#'
#' @param object A MappedNeuroVec object
#'
#' @importFrom crayon bold red green blue yellow silver
#' @importFrom utils object.size
#'
#' @keywords internal
setMethod("show", "MappedNeuroVec",
          function(object) {
            # Header
            cat("\n", crayon::bold(crayon::blue("MappedNeuroVec Object")), "\n")
            cat(crayon::silver("══════════════════════\n"))
            
            # Dimensions
            dims <- dim(object)
            spatial_dims <- paste(dims[1:3], collapse=" × ")
            temporal_dim <- if(length(dims) > 3) dims[4] else 1
            
            cat("\n", crayon::yellow("Dimensions:"), "\n")
            cat(" ", crayon::silver("•"), " Spatial: ", crayon::green(spatial_dims), "\n")
            cat(" ", crayon::silver("•"), " Temporal: ", crayon::green(temporal_dim), "\n")
            
            # Memory Mapping Information
            cat("\n", crayon::yellow("Memory Mapping:"), "\n")
            cat(" ", crayon::silver("•"), " Offset: ", crayon::green(object@offset), " elements\n")
            
            # Space Information
            sp <- space(object)
            spacing_str <- paste(round(spacing(sp), 2), collapse=" × ")
            origin_str <- paste(round(origin(sp), 2), collapse=" × ")
            
            cat("\n", crayon::yellow("Space Information:"), "\n")
            cat(" ", crayon::silver("•"), " Spacing: ", crayon::green(spacing_str), "\n")
            cat(" ", crayon::silver("•"), " Origin: ", crayon::green(origin_str), "\n")
            
            # Label
            if (!is.null(object@label) && object@label != "") {
              cat("\n", crayon::yellow("Label:"), "\n")
              cat(" ", crayon::silver("•"), " ", crayon::green(object@label), "\n")
            }
            
            cat("\n")
          })

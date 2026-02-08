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
#' @seealso \code{\link{MappedNeuroVec}} for creating instances of this class
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
#' \donttest{
#' # Create source from NIFTI file
#' source <- MappedNeuroVecSource(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Check dimensions
#' dim(source@meta_info)
#'
#' # View header information
#' str(source@meta_info)
#' }
#'
#' @seealso
#' \code{\link{MappedNeuroVec}} for the main user interface,
#' \code{\link{read_header}} for header reading details
#'
#' @importFrom methods new
#' @importFrom assertthat assert_that
#'
#' @export
#' @rdname MappedNeuroVecSource-class
MappedNeuroVecSource <- function(file_name) {
  assert_that(is.character(file_name) && length(file_name) == 1,
              msg = "'file_name' must be a single character string")
  assert_that(file.exists(file_name),
              msg = paste("File", file_name, "does not exist"))
  assert_that(file.access(file_name, mode = 4) == 0,
              msg = paste("File", file_name, "is not readable"))

  # Read and validate header
  meta <- read_header(file_name)
  assert_that(length(dim(meta)) >= 3,
              msg = "Input data must be at least 3-dimensional")

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
#' # Create mapped vector from NIFTI file
#' mvec <- MappedNeuroVec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Extract first volume
#' vol1 <- mvec[[1]]
#'
#' # Get dimensions
#' dim(mvec)
#'
#' # Access specific timepoint
#' timepoint <- mvec[, , , 2]
#'
#'
#' @seealso
#
#' \code{\link[mmap]{mmap}} for memory mapping details
#'
#' @importFrom mmap mmap mmapFlags
#' @importFrom methods new
#' @importFrom assertthat assert_that
#'
#' @export
#' @rdname MappedNeuroVec-class
MappedNeuroVec <- function(file_name, label = basename(file_name)) {
  assert_that(is.character(file_name) && length(file_name) == 1,
              msg = "'file_name' must be a single character string")
  assert_that(file.exists(file_name),
              msg = paste("File", file_name, "does not exist"))

  # Create source and load data
  src <- MappedNeuroVecSource(file_name)
  load_data(src)
}


#' @rdname load_data-methods
#' @export
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
#'
#' @rdname linear_access-methods
setMethod(f = "linear_access",
          signature = signature(x = "MappedNeuroVec", i = "numeric"),
          def = function(x, i) {
            assert_that(is.numeric(i) && !any(is.na(i)),
                        msg = "Index 'i' must be a numeric vector without NA values")

            # Calculate adjusted indices
            idx <- i + x@offset

            # Bounds checking
            assert_that(all(idx >= 1) && all(idx <= length(x@filemap)),
                        msg = "Index out of bounds")

            # Access mapped data
            x@filemap[idx]
          })


#' @export
#' @rdname show-methods
setMethod("show", "MappedNeuroVec",
          function(object) {
            # Header
            cat("\n", crayon::bold(crayon::blue("MappedNeuroVec Object")), "\n")
            cat(crayon::silver("====================\n"))

            # Dimensions
            dims <- dim(object)
            spatial_dims <- paste(dims[1:3], collapse=" x ")
            temporal_dim <- if(length(dims) > 3) dims[4] else 1

            cat("\n", crayon::yellow("Dimensions:"), "\n")
            cat(" ", crayon::silver("."), " Spatial: ", crayon::green(spatial_dims), "\n")
            cat(" ", crayon::silver("."), " Temporal: ", crayon::green(temporal_dim), "\n")

            # Memory Mapping Information
            cat("\n", crayon::yellow("Memory Mapping:"), "\n")
            cat(" ", crayon::silver("."), " Offset: ", crayon::green(object@offset), " elements\n")

            # Space Information
            sp <- space(object)
            spacing_str <- paste(round(spacing(sp), 2), collapse=" x ")
            origin_str <- paste(round(origin(sp), 2), collapse=" x ")

            cat("\n", crayon::yellow("Space Information:"), "\n")
            cat(" ", crayon::silver("."), " Spacing: ", crayon::green(spacing_str), "\n")
            cat(" ", crayon::silver("."), " Origin: ", crayon::green(origin_str), "\n")

            # Label
            if (!is.null(object@label) && object@label != "") {
              cat("\n", crayon::yellow("Label:"), "\n")
              cat(" ", crayon::silver("."), " ", crayon::green(object@label), "\n")
            }

            cat("\n")
          })

#' @export
#' @rdname as.matrix-methods
setMethod("as.matrix", "MappedNeuroVec",
  function(x) {
    dm <- dim(x)
    d123 <- prod(dm[1:3])
    d4 <- dm[4]
    tot <- d123 * d4
    vals <- linear_access(x, seq_len(tot))
    matrix(vals, d123, d4)
  }
)

#' @rdname mask-methods
#' @export
setMethod("mask", "MappedNeuroVec",
          function(x) {
            x@mask
          })


#' Convert to memory-mapped NeuroVec
#'
#' @description
#' Methods for the \code{\link{as_mmap}} generic, which convert various
#' neuroimaging vector types to a \code{\linkS4class{MappedNeuroVec}} backed
#' by an on-disk NIfTI file.
#'
#' @param x A neuroimaging vector (\code{NeuroVec}, \code{MappedNeuroVec},
#'   or \code{FileBackedNeuroVec}).
#' @param file Optional output file name. If \code{NULL}, a temporary file
#'   with extension \code{.nii} is created.
#' @param ... Additional arguments passed to methods (e.g. \code{data_type},
#'   \code{overwrite}).
#'
#' @return A \code{\linkS4class{MappedNeuroVec}} (or \code{x} itself if it is
#'   already memory-mapped).
#'
#' @name as_mmap-methods
#' @aliases as_mmap,MappedNeuroVec-method
#'          as_mmap,FileBackedNeuroVec-method
#'          as_mmap,NeuroVec-method
#'          as_mmap,SparseNeuroVec-method
#' @export
#' @rdname as_mmap
setMethod("as_mmap", signature(x = "MappedNeuroVec"),
          function(x, file = NULL, ...) {
            x
          })

#' @export
#' @rdname as_mmap
setMethod("as_mmap", signature(x = "FileBackedNeuroVec"),
          function(x, file = NULL, ...) {
            fn <- x@meta@data_file
            if (endsWith(fn, ".gz")) {
              stop("as_mmap: cannot memory-map gzipped file; re-save as uncompressed NIfTI with write_vec().")
            }
            MappedNeuroVec(fn)
          })

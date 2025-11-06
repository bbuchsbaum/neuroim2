#' Summarize Image Metadata
#'
#' @title Lightweight metadata for neuroimaging files
#' @description
#' `meta_info()` provides a simple, CRAN-friendly way to retrieve essential
#' image metadata without teaching S4 details up front. It accepts a file
#' path or a `FileMetaInfo` object and returns a normalized list containing
#' common fields like dimensions, spacing, origin, and transform.
#'
#' The function does not read image data; it only parses header information.
#'
#' @param x A character file path (e.g., `"image.nii.gz"`) or an object of
#'   class \code{\linkS4class{FileMetaInfo}}.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item `dim` Integer vector of image dimensions.
#'   \item `spacing` Numeric voxel spacing (mm).
#'   \item `origin` Numeric coordinate origin.
#'   \item `trans` 4x4 transformation matrix mapping grid to world (mm).
#'   \item `path` Data file path.
#'   \item `filename` Basename of `path`.
#'   \item `format` File format label (e.g., "NIFTI", "AFNI").
#'   \item `dtype` Storage data type label.
#'   \item `bytes_per_element` Bytes per element.
#'   \item `nvox` Number of voxels in the spatial volume (prod of first 3 dims).
#'   \item `nvol` Number of volumes (4th dim if present, else 1).
#'   \item `size_bytes` Approximate uncompressed size in bytes (`nvox * nvol * bytes_per_element`).
#'   \item `time_step` Time step (TR in seconds) if available for NIfTI, else `NA_real_`.
#' }
#'
#' @seealso \code{\link{read_header}}, \code{\link{trans}},
#'   \code{\linkS4class{FileMetaInfo}}, \code{\linkS4class{NIFTIMetaInfo}}
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
#' mi <- meta_info(f)
#' mi$dim
#' mi$spacing
#' mi$origin
#' mi$filename
#' # 4x4 transform
#' mi$trans
#' }
#'
#' @export
setGeneric("meta_info", function(x) standardGeneric("meta_info"))


# Internal helper to construct the normalized list from a FileMetaInfo
.meta_info_list <- function(mi) {
  # Basic fields
  dims <- as.integer(dim(mi))
  spacing <- as.numeric(mi@spacing)
  origin <- as.numeric(mi@origin)
  tx <- trans(mi)
  path <- mi@data_file
  filename <- basename(path)

  # Helpful extras
  fmt <- mi@descriptor@file_format
  dtype <- mi@data_type
  bpe <- as.integer(mi@bytes_per_element)

  # Handle dimensions robustly (ensure at least 3 for nvox)
  d3 <- dims[seq_len(min(3L, length(dims)))]
  # Use 1 for any missing spatial dimensions
  if (length(d3) < 3L) d3 <- c(d3, rep(1L, 3L - length(d3)))
  nvox <- as.integer(prod(d3))
  nvol <- if (length(dims) >= 4L) as.integer(dims[4L]) else 1L
  size_bytes <- as.numeric(nvox) * as.numeric(nvol) * as.numeric(bpe)

  # Optional: TR for NIfTI if available
  time_step <- NA_real_
  if (methods::is(mi, "NIFTIMetaInfo")) {
    hdr <- mi@header
    if (!is.null(hdr$pixdim) && is.numeric(hdr$pixdim) && length(hdr$pixdim) >= 5L) {
      time_step <- as.numeric(hdr$pixdim[5L])
    }
  }

  list(
    dim = dims,
    spacing = spacing,
    origin = origin,
    trans = tx,
    path = path,
    filename = filename,
    format = fmt,
    dtype = dtype,
    bytes_per_element = bpe,
    nvox = nvox,
    nvol = nvol,
    size_bytes = size_bytes,
    time_step = time_step
  )
}


#' @rdname meta_info
#' @export
setMethod("meta_info", signature(x = "FileMetaInfo"), function(x) {
  .meta_info_list(x)
})


#' @rdname meta_info
#' @export
setMethod("meta_info", signature(x = "character"), function(x) {
  if (length(x) == 0L) {
    stop("'x' must be a non-empty character vector of file paths")
  }
  if (length(x) == 1L) {
    mi <- read_header(x)
    return(.meta_info_list(mi))
  }
  # Vectorized: return a list of summaries
  lapply(x, function(p) {
    mi <- read_header(p)
    .meta_info_list(mi)
  })
})

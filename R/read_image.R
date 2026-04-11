#' read_image
#'
#' @description
#' Convenience wrapper that inspects the file header(s) and dispatches to the
#' appropriate specialized reader: \code{\link{read_vol}} for a single 3D file,
#' \code{\link{read_vec}} for 4D data or for any multi-file input, or
#' \code{\link{read_hyper_vec}} for 5D data.
#'
#' @details
#' Auto-dispatch (\code{type = "auto"}) uses the following rules:
#' \itemize{
#'   \item \code{length(file_name) > 1}: always routed to \code{\link{read_vec}},
#'         so the result is a \code{\linkS4class{NeuroVecSeq}} regardless of
#'         whether the individual files are 3D, 4D, or a mix. See
#'         \code{\link{read_vec}} for the return-type details.
#'   \item Single file with a 5th dimension \code{> 1}: routed to
#'         \code{\link{read_hyper_vec}}, returning a \code{\linkS4class{NeuroHyperVec}}.
#'   \item Single file with a 4th dimension \code{> 1}: routed to
#'         \code{\link{read_vec}}, returning a \code{\linkS4class{NeuroVec}}.
#'         If \code{index} is supplied (and \code{indices} is not), it is forwarded
#'         as \code{indices} so you can pull out a subset while still getting a
#'         \code{NeuroVec} back.
#'   \item Single effectively-3D file (either truly 3D or 4D with \code{dim[4] == 1}):
#'         routed to \code{\link{read_vol}}, returning a \code{\linkS4class{DenseNeuroVol}}.
#' }
#'
#' Explicit \code{type} values bypass header inspection:
#' \itemize{
#'   \item \code{type = "vol"}: requires a single file; always returns a
#'         \code{\linkS4class{DenseNeuroVol}}. The \code{index} argument picks the
#'         sub-volume when the file is 4D.
#'   \item \code{type = "vec"}: forwards to \code{\link{read_vec}}. Multi-file input
#'         yields a \code{\linkS4class{NeuroVecSeq}}. A single 3D file is promoted
#'         to a \code{NeuroVec} with \code{dim[4] == 1}.
#'   \item \code{type = "hyper"}: requires a single file; returns a
#'         \code{\linkS4class{NeuroHyperVec}}.
#' }
#'
#' @param file_name Character vector of one or more file paths.
#' @param type One of \code{"auto"}, \code{"vol"}, \code{"vec"}, or \code{"hyper"}
#'   to override header-based dispatch.
#' @param index Integer volume index. Used as the single-volume selector for
#'   \code{\link{read_vol}} (when dispatching to it) and, when \code{indices} is
#'   \code{NULL}, forwarded to \code{\link{read_vec}} as \code{indices} so you can
#'   pull out a subset of time points while still receiving a \code{NeuroVec}.
#' @param indices Optional integer vector of sub-volume indices forwarded to
#'   \code{\link{read_vec}}. Takes precedence over \code{index}.
#' @param mask Optional spatial mask passed through to the vector or hyper-vector
#'   readers. Ignored by \code{read_vol}.
#' @param mode IO mode forwarded to \code{\link{read_vec}}; see that function for
#'   details. Ignored for 3D and 5D dispatch.
#'
#' @return The return type depends on dispatch:
#' \itemize{
#'   \item \strong{3D dispatch} (single effectively-3D file, or \code{type = "vol"}):
#'         a \code{\linkS4class{DenseNeuroVol}}.
#'   \item \strong{4D dispatch} (single 4D file, or \code{type = "vec"} with one file):
#'         a \code{\linkS4class{NeuroVec}} (concrete subclass depends on \code{mode}).
#'   \item \strong{Multi-file dispatch} (\code{length(file_name) > 1}, or
#'         \code{type = "vec"} with multiple files): a
#'         \code{\linkS4class{NeuroVecSeq}}, which extends \code{NeuroVec} but stores
#'         each file as a distinct segment rather than concatenating into a single
#'         contiguous 4D array. Mixed 3D/4D inputs are allowed; each 3D file
#'         contributes one time point to the sequence.
#'   \item \strong{5D dispatch} (single 5D file, or \code{type = "hyper"}):
#'         a \code{\linkS4class{NeuroHyperVec}}.
#' }
#'
#' @seealso \code{\link{read_vol}}, \code{\link{read_vec}}, \code{\link{read_hyper_vec}}
#'
#' @examples
#' # 3D file -> DenseNeuroVol
#' vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
#'
#' # 4D file -> NeuroVec
#' vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Multiple files (any mix of 3D and 4D) -> NeuroVecSeq
#' fn <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
#' seq_vec <- read_image(c(fn, fn))
#' class(seq_vec)  # "NeuroVecSeq"
#'
#' @export
read_image <- function(file_name,
                       type = c("auto", "vol", "vec", "hyper"),
                       index = 1,
                       indices = NULL,
                       mask = NULL,
                       mode = c("normal", "mmap", "bigvec", "filebacked")) {
  type <- match.arg(type)

  if (type == "vol") {
    if (length(file_name) != 1L) {
      stop("read_image: type='vol' expects a single file_name")
    }
    return(read_vol(file_name, index = index))
  }

  if (type == "vec") {
    if (length(file_name) == 1L && is.null(indices) && !missing(index)) {
      indices <- index
    }
    return(read_vec(file_name, indices = indices, mask = mask, mode = mode))
  }

  if (type == "hyper") {
    if (length(file_name) != 1L) {
      stop("read_image: type='hyper' expects a single file_name")
    }
    return(read_hyper_vec(file_name, mask = mask))
  }

  if (length(file_name) > 1L) {
    return(read_vec(file_name, indices = indices, mask = mask, mode = mode))
  }

  meta <- read_header(file_name)
  dims <- dim(meta)
  extra_hyper_dims <- if (length(dims) > 4L) dims[-(1:4)] else integer(0)
  extra_dims <- if (length(dims) > 3L) dims[-(1:3)] else integer(0)
  has_hyper_dim <- length(extra_hyper_dims) > 0L && any(extra_hyper_dims > 1L)
  has_vec_dim <- length(extra_dims) > 0L && any(extra_dims > 1L)

  if (has_hyper_dim) {
    read_hyper_vec(file_name, mask = mask)
  } else if (has_vec_dim) {
    if (is.null(indices) && !missing(index)) {
      indices <- index
    }
    read_vec(file_name, indices = indices, mask = mask, mode = mode)
  } else {
    read_vol(file_name, index = index)
  }
}


#' Read a 5D image as a NeuroHyperVec
#'
#' @param file_name Path to a single NIfTI file.
#' @param mask Optional spatial mask (logical array/vector, \code{NeuroVol}, or
#'   \code{LogicalNeuroVol}). When provided, only masked voxels are stored.
#'
#' @return A \code{\linkS4class{NeuroHyperVec}}.
#' @export
read_hyper_vec <- function(file_name, mask = NULL) {
  if (!is.character(file_name) || length(file_name) != 1L) {
    stop("'file_name' must be a single character string")
  }
  if (!file.exists(file_name)) {
    stop("file does not exist: ", file_name)
  }

  meta <- read_header(file_name)
  dims <- dim(meta)
  if (length(dims) != 5L) {
    stop("read_hyper_vec currently requires a 5D image, got dimensions: ",
         paste(dims, collapse = " x "))
  }

  spatial_dims <- dims[1:3]
  ntrials <- dims[4]
  nfeatures <- dims[5]

  spatial_space <- NeuroSpace(
    spatial_dims,
    spacing = meta@spacing,
    origin = meta@origin,
    axes = meta@spatial_axes,
    trans = trans(meta)
  )

  mask_arr <- NULL
  if (is.null(mask)) {
    mask_arr <- array(TRUE, dim = spatial_dims)
  } else if (inherits(mask, "LogicalNeuroVol")) {
    mask_arr <- as(mask, "array")
  } else if (inherits(mask, "NeuroVol")) {
    mask_arr <- as(mask, "array") != 0
  } else if (is.array(mask)) {
    mask_arr <- as.logical(mask)
  } else if (is.vector(mask)) {
    if (length(mask) != prod(spatial_dims)) {
      stop("mask vector length does not match spatial dimensions")
    }
    mask_arr <- array(as.logical(mask), dim = spatial_dims)
  } else {
    stop("'mask' must be NULL, an array/vector, NeuroVol, or LogicalNeuroVol")
  }

  if (!all(dim(mask_arr) == spatial_dims)) {
    stop("mask dimensions must match spatial dimensions: ",
         paste(spatial_dims, collapse = " x "))
  }
  if (any(is.na(mask_arr))) {
    stop("mask cannot contain NA values")
  }

  mask_vol <- LogicalNeuroVol(mask_arr, spatial_space)
  mask_idx <- which(mask_vol@.Data)
  nvox <- length(mask_idx)

  data <- array(0, dim = c(nfeatures, ntrials, nvox))
  nels <- prod(spatial_dims)
  reader <- data_reader(meta, offset = 0)
  on.exit(close(reader), add = TRUE)

  vol_index <- 1L
  for (f in seq_len(nfeatures)) {
    for (t in seq_len(ntrials)) {
      vals <- read_elements(reader, nels)
      vals <- .apply_data_scaling(vals, meta, index = vol_index)
      data[f, t, ] <- vals[mask_idx]
      vol_index <- vol_index + 1L
    }
  }

  full_space <- NeuroSpace(
    dim = dims,
    spacing = meta@spacing,
    origin = meta@origin,
    axes = meta@spatial_axes,
    trans = trans(meta)
  )

  NeuroHyperVec(data, full_space, mask_vol)
}

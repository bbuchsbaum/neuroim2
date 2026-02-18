#' read_image
#'
#' Convenience wrapper that inspects the file metadata and dispatches to
#' \code{\link{read_vol}} for 3D data, \code{\link{read_vec}} for 4D data,
#' or \code{\link{read_hyper_vec}} for 5D data.
#'
#' @param file_name Character vector of file paths.
#' @param type One of \code{"auto"}, \code{"vol"}, \code{"vec"}, or \code{"hyper"}
#'   to override dispatch.
#' @param index Volume index to use when returning a \code{\linkS4class{NeuroVol}} or
#'   when you want to load a subset of volumes while still returning a \code{NeuroVec}.
#' @param indices Optional vector of indices passed through to \code{\link{read_vec}}.
#' @param mask Optional spatial mask passed through to vector/hyper-vector readers.
#' @param mode IO mode forwarded to \code{\link{read_vec}}.
#'
#' @return A \code{\linkS4class{NeuroVol}} when the input is effectively 3D (or when
#'   \code{type = "vol"}), a \code{\linkS4class{NeuroVec}}/\code{\linkS4class{NeuroVecSeq}}
#'   for 4D input, or a \code{\linkS4class{NeuroHyperVec}} for 5D input.
#'
#' @examples
#' vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
#' vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
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

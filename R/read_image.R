#' read_image
#'
#' Convenience wrapper that inspects the file metadata and dispatches to
#' \code{\link{read_vol}} for 3D data or \code{\link{read_vec}} for 4D data.
#'
#' @param file_name Character vector of file paths.
#' @param type One of \code{"auto"}, \code{"vol"}, or \code{"vec"} to override dispatch.
#' @param index Volume index to use when returning a \code{\linkS4class{NeuroVol}} or
#'   when you want to load a subset of volumes while still returning a \code{NeuroVec}.
#' @param indices Optional vector of indices passed through to \code{\link{read_vec}}.
#' @param mask Optional mask passed to \code{\link{read_vec}}.
#' @param mode IO mode forwarded to \code{\link{read_vec}}.
#'
#' @return A \code{\linkS4class{NeuroVol}} when the input is effectively 3D (or when
#'   \code{type = "vol"}), otherwise a \code{\linkS4class{NeuroVec}}/\code{\linkS4class{NeuroVecSeq}}.
#'
#' @examples
#' vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
#' vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' @export
read_image <- function(file_name,
                       type = c("auto", "vol", "vec"),
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

  if (length(file_name) > 1L) {
    return(read_vec(file_name, indices = indices, mask = mask, mode = mode))
  }

  meta <- read_header(file_name)
  dims <- dim(meta)
  extra_dims <- if (length(dims) > 3L) dims[-(1:3)] else integer(0)
  has_vec_dim <- length(extra_dims) > 0L && any(extra_dims > 1L)

  if (has_vec_dim) {
    if (is.null(indices) && !missing(index)) {
      indices <- index
    }
    read_vec(file_name, indices = indices, mask = mask, mode = mode)
  } else {
    read_vol(file_name, index = index)
  }
}

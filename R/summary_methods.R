#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Summary of Neuroimaging Objects
#'
#' @description
#' Provides a concise summary of neuroimaging volume and vector objects,
#' including spatial metadata and data statistics.
#'
#' @param object A neuroimaging object.
#' @param ... Additional arguments (currently ignored).
#' @return An object of class \code{"summary.NeuroVol"} or
#'   \code{"summary.NeuroVec"} (invisibly), with a print method.
#'
#' @examples
#' vol <- DenseNeuroVol(array(rnorm(27), c(3,3,3)),
#'                      NeuroSpace(c(3L,3L,3L), c(1,1,1)))
#' summary(vol)
#'
#' @name summary-neuro-methods
#' @rdname summary-neuro-methods
#' @aliases summary,NeuroVol-method summary,DenseNeuroVec-method summary,SparseNeuroVec-method
#' @export
setMethod("summary", signature(object = "NeuroVol"), function(object, ...) {
  d <- dim(object)
  vals <- as.numeric(as.dense(object)@.Data)
  sp <- space(object)

  out <- list(
    class = class(object)[1],
    dim = d,
    spacing = spacing(object),
    origin = origin(object),
    orientation = safe_axcodes(sp),
    min = min(vals, na.rm = TRUE),
    max = max(vals, na.rm = TRUE),
    mean = mean(vals, na.rm = TRUE),
    sd = stats::sd(vals, na.rm = TRUE),
    NAs = sum(is.na(vals)),
    zeros = sum(vals == 0, na.rm = TRUE),
    nonzeros = sum(vals != 0, na.rm = TRUE)
  )
  class(out) <- "summary.NeuroVol"
  out
})

#' @rdname summary-neuro-methods
#' @export
setMethod("summary", signature(object = "DenseNeuroVec"), function(object, ...) {
  d <- dim(object)
  M <- matrix(object@.Data, nrow = prod(d[1:3]), ncol = d[4])
  sp <- space(object)

  rmeans <- rowMeans(M)
  rsds <- if (d[4] > 1) sqrt(rowSums((M - rmeans)^2) / (d[4] - 1L)) else rep(0, nrow(M))

  out <- list(
    class = class(object)[1],
    dim = d,
    spacing = spacing(object),
    origin = origin(object),
    orientation = safe_axcodes(sp),
    time_points = d[4],
    global_min = min(M, na.rm = TRUE),
    global_max = max(M, na.rm = TRUE),
    global_mean = mean(M, na.rm = TRUE),
    temporal_mean_range = range(rmeans, na.rm = TRUE),
    temporal_sd_range = range(rsds, na.rm = TRUE),
    NAs = sum(is.na(M)),
    nonzero_voxels = sum(rmeans != 0, na.rm = TRUE),
    total_voxels = prod(d[1:3])
  )
  class(out) <- "summary.NeuroVec"
  out
})

#' @rdname summary-neuro-methods
#' @export
setMethod("summary", signature(object = "SparseNeuroVec"), function(object, ...) {
  d <- dim(object)
  sp <- space(object)

  cmeans <- colMeans(object@data)
  csds <- if (d[4] > 1) {
    sqrt(colSums((object@data - rep(cmeans, each = nrow(object@data)))^2) / (d[4] - 1L))
  } else {
    rep(0, ncol(object@data))
  }

  out <- list(
    class = class(object)[1],
    dim = d,
    spacing = spacing(object),
    origin = origin(object),
    orientation = safe_axcodes(sp),
    time_points = d[4],
    global_min = min(object@data, na.rm = TRUE),
    global_max = max(object@data, na.rm = TRUE),
    global_mean = mean(object@data, na.rm = TRUE),
    temporal_mean_range = range(cmeans, na.rm = TRUE),
    temporal_sd_range = range(csds, na.rm = TRUE),
    NAs = sum(is.na(object@data)),
    nonzero_voxels = length(indices(object)),
    total_voxels = prod(d[1:3])
  )
  class(out) <- "summary.NeuroVec"
  out
})

#' @method print summary.NeuroVol
#' @export
print.summary.NeuroVol <- function(x, ...) {
  cat(sprintf("<%s> [%s]\n", x$class, paste(x$dim, collapse = " x ")))
  cat(sprintf("  Spacing    : %s mm\n", paste(round(x$spacing, 3), collapse = " x ")))
  cat(sprintf("  Origin     : %s\n", paste(round(x$origin, 2), collapse = ", ")))
  cat(sprintf("  Orientation: %s\n", x$orientation))
  cat(sprintf("  Range      : [%.4g, %.4g]\n", x$min, x$max))
  cat(sprintf("  Mean (SD)  : %.4g (%.4g)\n", x$mean, x$sd))
  if (x$NAs > 0) cat(sprintf("  NAs        : %d\n", x$NAs))
  cat(sprintf("  Non-zero   : %d / %d voxels\n", x$nonzeros, x$nonzeros + x$zeros))
  invisible(x)
}

#' @method print summary.NeuroVec
#' @export
print.summary.NeuroVec <- function(x, ...) {
  cat(sprintf("<%s> [%s]\n", x$class, paste(x$dim, collapse = " x ")))
  cat(sprintf("  Spacing     : %s mm\n", paste(round(x$spacing, 3), collapse = " x ")))
  cat(sprintf("  Origin      : %s\n", paste(round(x$origin, 2), collapse = ", ")))
  cat(sprintf("  Orientation : %s\n", x$orientation))
  cat(sprintf("  Time points : %d\n", x$time_points))
  cat(sprintf("  Global range: [%.4g, %.4g]\n", x$global_min, x$global_max))
  cat(sprintf("  Global mean : %.4g\n", x$global_mean))
  cat(sprintf("  Temporal mean range: [%.4g, %.4g]\n", x$temporal_mean_range[1], x$temporal_mean_range[2]))
  cat(sprintf("  Temporal SD range  : [%.4g, %.4g]\n", x$temporal_sd_range[1], x$temporal_sd_range[2]))
  if (x$NAs > 0) cat(sprintf("  NAs         : %d\n", x$NAs))
  cat(sprintf("  Non-zero    : %d / %d voxels\n", x$nonzero_voxels, x$total_voxels))
  invisible(x)
}

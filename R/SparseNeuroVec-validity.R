#' @include all_class.R
NULL

# Validity checks for SparseNeuroVec
setValidity("SparseNeuroVec", function(object) {
  d  <- dim(object)
  if (length(d) != 4L) {
    return("SparseNeuroVec requires a 4D space (x,y,z,time).")
  }
  nt   <- d[4L]
  nvox <- prod(d[1:3])

  # data should be a matrix-like object: rows = time, cols = nonzero voxels
  dat <- object@data
  if (!(is.matrix(dat) || inherits(dat, "Matrix"))) {
    return("SparseNeuroVec data must be a matrix (or Matrix::Matrix).")
  }
  if (nrow(dat) != nt) {
    return(sprintf("Data/time mismatch: nrow(data) = %d, but space timepoints = %d.", nrow(dat), nt))
  }

  # mask must be logical 3D volume with length equal to nvox
  mv <- mask(object)  # LogicalNeuroVol
  mlog <- as.array(mv)
  if (length(mlog) != nvox) {
    return(sprintf("Mask/space mismatch: length(mask) = %d, but prod(space[1:3]) = %d.", length(mlog), nvox))
  }
  n_nonzero <- sum(mlog)
  if (ncol(dat) != n_nonzero) {
    return(sprintf("Data/mask mismatch: ncol(data) = %d, but sum(mask) = %d.", ncol(dat), n_nonzero))
  }
  TRUE
})


#' Normalize a mask to a logical array
#'
#' Coerces various mask representations into a logical 3D array of the
#' requested dimensions. Accepts: logical arrays, integer index vectors,
#' \code{LogicalNeuroVol}, numeric \code{NeuroVol} (thresholded > 0),
#' logical vectors, or \code{NULL} (all \code{TRUE}).
#'
#' @param mask The mask input (see description).
#' @param target_dim Integer vector of length 3 giving the target dimensions.
#' @return A logical array of dimension \code{target_dim}.
#' @keywords internal
normalize_mask <- function(mask, target_dim) {
  D <- as.integer(target_dim)

  if (is.null(mask)) {
    return(array(TRUE, D))
  }

  # Already a LogicalNeuroVol — extract the array

  if (inherits(mask, "LogicalNeuroVol")) {
    if (!identical(dim(mask)[1:3], D)) {
      stop(sprintf("mask dimensions [%s] do not match target [%s]",
                   paste(dim(mask), collapse = "x"), paste(D, collapse = "x")))
    }
    return(as.array(mask))
  }

  # Numeric NeuroVol — threshold > 0
  if (inherits(mask, "NeuroVol")) {
    if (!identical(dim(mask)[1:3], D)) {
      stop(sprintf("mask dimensions [%s] do not match target [%s]",
                   paste(dim(mask), collapse = "x"), paste(D, collapse = "x")))
    }
    return(as.array(mask) > 0)
  }

  # Short numeric vector — treat as 1D indices into the volume
  if (is.vector(mask) && is.numeric(mask) && !is.logical(mask) && length(mask) < prod(D)) {
    m <- array(FALSE, D)
    m[mask] <- TRUE
    return(m)
  }

  # Array with matching dims
  if (is.array(mask) && identical(dim(mask), D)) {
    return(array(as.logical(mask), D))
  }

  # Flat vector of length prod(D)
  if (is.vector(mask) && length(mask) == prod(D)) {
    return(array(as.logical(mask), D))
  }

  stop(sprintf("Cannot coerce mask of class '%s' (length %d) to logical array [%s]",
               class(mask)[1], length(mask), paste(D, collapse = "x")))
}

#' @keywords internal
#' @noRd
.same_spatial_space <- function(x, y) {
  identical(as.integer(dim(x)), as.integer(dim(y))) &&
    isTRUE(all.equal(as.numeric(spacing(x)),
                     as.numeric(spacing(y)),
                     check.attributes = FALSE)) &&
    isTRUE(all.equal(as.numeric(origin(x)),
                     as.numeric(origin(y)),
                     check.attributes = FALSE)) &&
    isTRUE(all.equal(unname(as.matrix(trans(x))),
                     unname(as.matrix(trans(y))),
                     check.attributes = FALSE))
}

#' @keywords internal
#' @noRd
.coerce_spatial_mask <- function(mask, target_space) {
  if (!inherits(target_space, "NeuroSpace") || ndim(target_space) != 3L) {
    cli::cli_abort("{.arg target_space} must be a 3D {.cls NeuroSpace}.")
  }

  if (inherits(mask, "NeuroVol")) {
    if (!.same_spatial_space(space(mask), target_space)) {
      cli::cli_abort("Spatial geometry of {.arg mask} must match the target image space.")
    }
  }

  LogicalNeuroVol(normalize_mask(mask, dim(target_space)), target_space)
}

#' @keywords internal
#' @noRd
.representative_volume_from_matrix <- function(mat, dims3, representative = c("median", "mean_abs")) {
  representative <- match.arg(representative)
  vals <- representative_volume_cpp(unname(as.matrix(mat)), representative)

  array(as.numeric(vals), dim = as.integer(dims3))
}

#' @keywords internal
#' @noRd
.afni_clip_level_numeric <- function(x, mfrac = 0.5, nhist = 10000L) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals) & vals > 0]

  if (length(vals) <= 222L) {
    return(0)
  }

  if (!is.finite(mfrac) || mfrac <= 0 || mfrac >= 0.99) {
    mfrac <- 0.5
  }

  vmax <- max(vals)
  if (!is.finite(vmax) || vmax < 1e-100) {
    return(0)
  }

  is_integerish <- all(abs(vals - round(vals)) < 1e-8)
  use_integer_hist <- is_integerish && vmax <= 32767

  if (use_integer_hist) {
    nhist_eff <- if (vmax <= 255) 255L else 32767L
    sfac <- 1
    bins <- as.integer(round(vals))
    bins <- bins[bins >= 0L & bins <= nhist_eff]
  } else {
    nhist_eff <- as.integer(nhist)
    sfac <- nhist_eff / vmax
    bins <- as.integer(floor(sfac * vals + 0.499))
    bins <- bins[bins >= 0L & bins <= nhist_eff]
  }

  if (length(bins) <= 222L) {
    return(0)
  }

  hist <- tabulate(bins + 1L, nbins = nhist_eff + 1L)
  npos <- sum(hist)
  dsum <- sum(as.double(bins) * as.double(bins))

  qq <- as.integer(0.65 * npos)
  ib <- as.integer(round(0.5 * sqrt(dsum / npos)))
  acc <- 0L
  ii <- nhist_eff - 1L

  while (ii >= ib && acc < qq) {
    acc <- acc + hist[ii + 1L]
    ii <- ii - 1L
  }

  ncut <- max(ii, 0L)
  iter <- 0L

  repeat {
    start <- max(ncut, 0L) + 1L
    npos_cut <- if (start <= nhist_eff) sum(hist[start:nhist_eff]) else 0L
    nhalf <- npos_cut %/% 2L

    acc <- 0L
    ii <- max(ncut, 0L)
    while (ii < nhist_eff && acc < nhalf) {
      acc <- acc + hist[ii + 1L]
      ii <- ii + 1L
    }

    nold <- ncut
    ncut <- as.integer(mfrac * ii)
    iter <- iter + 1L

    if (iter >= 66L || ncut == nold) {
      break
    }
  }

  out <- as.numeric(ncut) / sfac
  min(out, 1e38)
}

#' @keywords internal
#' @noRd
.cmass_3d_zero_based <- function(arr) {
  dims <- dim(arr)
  w <- as.numeric(arr)
  w[!is.finite(w) | w < 0] <- 0

  if (!any(w > 0)) {
    return((as.numeric(dims) - 1) / 2)
  }

  idx <- arrayInd(seq_along(w), dims)
  c(
    stats::weighted.mean(idx[, 1] - 1, w),
    stats::weighted.mean(idx[, 2] - 1, w),
    stats::weighted.mean(idx[, 3] - 1, w)
  )
}

#' @keywords internal
#' @noRd
.afni_gradual_clip_array <- function(arr, mfrac = 0.5) {
  dims <- as.integer(dim(arr))
  stopifnot(length(dims) == 3L)

  nx <- dims[1]
  ny <- dims[2]
  nz <- dims[3]

  cm <- .cmass_3d_zero_based(arr)
  ic <- min(max(as.integer(round(cm[1])), 0L), nx - 1L)
  jc <- min(max(as.integer(round(cm[2])), 0L), ny - 1L)
  kc <- min(max(as.integer(round(cm[3])), 0L), nz - 1L)

  it <- nx - 1L
  jt <- ny - 1L
  kt <- nz - 1L

  val_floor <- 0.333 * .afni_clip_level_numeric(arr, mfrac)

  ox <- max(1L, as.integer(round(0.01 * nx)))
  oy <- max(1L, as.integer(round(0.01 * ny)))
  oz <- max(1L, as.integer(round(0.01 * nz)))

  icm <- max(ic - ox, 0L)
  icp <- min(ic + ox, it)
  jcm <- max(jc - oy, 0L)
  jcp <- min(jc + oy, jt)
  kcm <- max(kc - oz, 0L)
  kcp <- min(kc + oz, kt)

  octclip <- function(xa, xb, ya, yb, za, zb) {
    max(
      .afni_clip_level_numeric(
        arr[(xa + 1L):(xb + 1L), (ya + 1L):(yb + 1L), (za + 1L):(zb + 1L)],
        mfrac = mfrac
      ),
      val_floor
    )
  }

  c000 <- octclip(0L,  icp, 0L,  jcp, 0L,  kcp)
  c100 <- octclip(icm, it,  0L,  jcp, 0L,  kcp)
  c010 <- octclip(0L,  icp, jcm, jt,  0L,  kcp)
  c110 <- octclip(icm, it,  jcm, jt,  0L,  kcp)
  c001 <- octclip(0L,  icp, 0L,  jcp, kcm, kt)
  c101 <- octclip(icm, it,  0L,  jcp, kcm, kt)
  c011 <- octclip(0L,  icp, jcm, jt,  kcm, kt)
  c111 <- octclip(icm, it,  jcm, jt,  kcm, kt)

  x0 <- 0.5 * ic
  x1 <- 0.5 * (ic + it)
  y0 <- 0.5 * jc
  y1 <- 0.5 * (jc + jt)
  z0 <- 0.5 * kc
  z1 <- 0.5 * (kc + kt)

  dxi <- if (x1 > x0) 1 / (x1 - x0) else 0
  dyi <- if (y1 > y0) 1 / (y1 - y0) else 0
  dzi <- if (z1 > z0) 1 / (z1 - z0) else 0

  xw1 <- pmin(1, pmax(0, ((0:(nx - 1L)) - x0) * dxi))
  yw1 <- pmin(1, pmax(0, ((0:(ny - 1L)) - y0) * dyi))
  zw1 <- pmin(1, pmax(0, ((0:(nz - 1L)) - z0) * dzi))
  xw0 <- 1 - xw1
  yw0 <- 1 - yw1
  zw0 <- 1 - zw1

  outer3 <- function(x, y, z) {
    xy <- outer(x, y)
    array(rep(as.vector(xy), times = length(z)), dim = c(length(x), length(y), length(z))) *
      array(rep(z, each = length(x) * length(y)), dim = c(length(x), length(y), length(z)))
  }

  c000 * outer3(xw0, yw0, zw0) +
    c100 * outer3(xw1, yw0, zw0) +
    c010 * outer3(xw0, yw1, zw0) +
    c110 * outer3(xw1, yw1, zw0) +
    c001 * outer3(xw0, yw0, zw1) +
    c101 * outer3(xw1, yw0, zw1) +
    c011 * outer3(xw0, yw1, zw1) +
    c111 * outer3(xw1, yw1, zw1)
}

#' @keywords internal
#' @noRd
.shift_array_clamp <- function(arr, dx = 0L, dy = 0L, dz = 0L) {
  dims <- dim(arr)
  ix <- pmin(pmax(seq_len(dims[1]) + dx, 1L), dims[1])
  iy <- pmin(pmax(seq_len(dims[2]) + dy, 1L), dims[2])
  iz <- pmin(pmax(seq_len(dims[3]) + dz, 1L), dims[3])
  arr[ix, iy, iz, drop = FALSE]
}

#' @keywords internal
#' @noRd
.neighbor_count_18 <- function(mask) {
  offsets <- as.matrix(expand.grid(dx = -1:1, dy = -1:1, dz = -1:1))
  shell <- rowSums(abs(offsets))
  offsets <- offsets[shell > 0 & shell <= 2, , drop = FALSE]

  counts <- array(0L, dim(mask))
  for (ii in seq_len(nrow(offsets))) {
    counts <- counts + .shift_array_clamp(mask, offsets[ii, 1], offsets[ii, 2], offsets[ii, 3])
  }
  counts
}

#' @keywords internal
#' @noRd
.largest_component_mask <- function(mask, connect = c("26-connect", "18-connect", "6-connect")) {
  connect <- match.arg(connect)
  if (!any(mask)) {
    return(array(FALSE, dim(mask)))
  }

  cc <- conn_comp_3D(mask, connect = connect)
  labs <- cc$index
  ids <- labs[labs > 0]

  if (length(ids) == 0L) {
    return(array(FALSE, dim(mask)))
  }

  keep <- which.max(tabulate(ids))
  array(labs == keep, dim = dim(mask))
}

#' @keywords internal
#' @noRd
.fill_holes_mask <- function(mask) {
  if (!any(mask)) {
    return(mask)
  }

  outside <- !mask
  if (!any(outside)) {
    return(mask)
  }

  cc <- conn_comp_3D(outside, connect = "6-connect")
  labs <- cc$index
  dims <- dim(mask)
  edge_ids <- unique(c(
    labs[1, , ], labs[dims[1], , ],
    labs[, 1, ], labs[, dims[2], ],
    labs[, , 1], labs[, , dims[3]]
  ))
  edge_ids <- edge_ids[edge_ids > 0]

  outside_keep <- array(labs %in% edge_ids, dim = dim(mask))
  array(!outside_keep, dim = dim(mask))
}

#' @keywords internal
#' @noRd
.peel_restore_mask <- function(mask, peels = 1L, peel_threshold = 17L) {
  peels <- as.integer(peels)
  peel_threshold <- as.integer(peel_threshold)

  if (peels < 1L || !any(mask)) {
    return(mask)
  }

  real_peel_threshold <- min(18L, max(1L, peel_threshold))
  marks <- array(0L, dim(mask))
  out <- mask

  for (pp in seq_len(peels)) {
    counts <- .neighbor_count_18(out)
    to_erode <- out & (counts < real_peel_threshold)
    marks[to_erode] <- pp
    out[to_erode] <- FALSE
  }

  for (pp in seq.int(peels, 1L)) {
    counts <- .neighbor_count_18(out)
    bth <- if (pp == peels) 0L else 1L
    to_restore <- (marks >= pp) & !out & (counts > bth)
    out[to_restore] <- TRUE
  }

  out
}

#' @keywords internal
#' @noRd
.automask_array <- function(arr,
                            mfrac = 0.5,
                            gradual = TRUE,
                            peels = 1L,
                            peel_threshold = 17L,
                            connect = c("26-connect", "18-connect", "6-connect")) {
  connect <- match.arg(connect)
  arr <- as.array(arr)
  dims <- dim(arr)

  if (length(dims) != 3L) {
    cli::cli_abort("{.arg arr} must be a 3D array.")
  }

  clip <- .afni_clip_level_numeric(arr, mfrac = mfrac)
  if (!is.finite(clip) || clip <= 0) {
    return(array(FALSE, dims))
  }

  thr <- if (isTRUE(gradual)) .afni_gradual_clip_array(arr, mfrac = mfrac) else array(clip, dims)
  mask <- arr >= thr

  if (!any(mask)) {
    return(mask)
  }

  mask <- .largest_component_mask(mask, connect = connect)

  if (peels > 0L && any(mask)) {
    mask <- .peel_restore_mask(mask, peels = peels, peel_threshold = peel_threshold)
    if (any(mask)) {
      mask <- .largest_component_mask(mask, connect = connect)
    }
  }

  if (any(mask)) {
    mask <- .fill_holes_mask(mask)
    mask <- .largest_component_mask(mask, connect = connect)
  }

  mask
}


#' @export
#' @rdname read_columns-methods
setMethod(f="read_columns", signature=c(x="ColumnReader", column_indices="integer"),
          def=function(x, column_indices) {
            x@reader(column_indices)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="integer", FUN="function"),
          def=function(x, fac, FUN) {
            callGeneric(x,as.factor(fac), FUN)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="factor", FUN="missing"),
          def=function(x, fac) {
            split_reduce(x, fac, mean)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            if (length(fac) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }

            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(names(ind), function(lev) {
              future.apply::future_apply(x[ind[[lev]],,drop=FALSE], 2, FUN)
            }))

            row.names(out) <- levels(fac)
            out
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "NeuroVec", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            split_by_voxel <- if (length(fac) == prod(dim(x)[1:3])) {
              TRUE
            } else if (length(fac) == dim(x)[4]) {
              FALSE
            } else {
              stop("length of 'fac' must be equal to number of voxels or to number of volumes")
            }

            if (split_by_voxel) {
              ind <- split(seq_along(fac), fac)
              out <- do.call(rbind, lapply(names(ind), function(lev) {
                mat <- series(x, ind[[lev]])
                apply(mat, 1, FUN)
              }))

              row.names(out) <- levels(fac)
              out
            } else {
              m <- as.matrix(x)
              callGeneric(m, fac, FUN)
            }
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "NeuroVec", fac="factor", FUN="missing"),
          def=function(x, fac, FUN) {
            split_by_voxel <- if (length(fac) == prod(dim(x)[1:3])) {
              TRUE
            } else if (length(fac) == dim(x)[4]) {
              FALSE
            } else {
              stop("length of 'fac' must be equal to number of voxels or to number of volumes")
            }

            if (split_by_voxel) {
              ind <- split(seq_along(fac), fac)
              out <- do.call(rbind, lapply(names(ind), function(lev) {
                rowMeans(series(x, ind[[lev]]))
              }))

              row.names(out) <- levels(fac)
              out
            } else {
              m <- as.matrix(x)
              callGeneric(m, fac)
            }
          })




#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "matrix", f="factor", center="logical", scale="logical"),
          def=function(x, f, center=TRUE, scale=TRUE) {
            if (length(f) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }

            out <- matrix(0, nrow(x), ncol(x))
            ind <- split(seq_along(f), f)

            for (lev in names(ind)) {
              keep <- ind[[lev]]
              xs <- base::scale(x[keep,,drop=FALSE], center=center, scale=scale)
              out[keep,] <- xs
            }

            out
          })



#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "matrix", f="factor", center="missing", scale="missing"),
          def=function(x, f) {
            callGeneric(x,f, TRUE, TRUE)
          })


#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="logical", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, center, TRUE)
          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="missing", scale="missing"),
          def=function(x, f) {
            callGeneric(x, f, TRUE, TRUE)

          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="logical", scale="missing"),
          def=function(x, f, center) {
            callGeneric(x, f, center, TRUE)

          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="logical", scale="logical"),
          def=function(x, f, center, scale) {
            m <- callGeneric(t(as.matrix(x)), f, center, scale)
            NeuroVec(m, space(x))
          })



#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="DenseNeuroVec", center="logical", scale="logical"),
          def=function(x, center, scale) {
            d <- dim(x)
            nv <- prod(d[1:3])
            nt <- d[4]
            # Reshape directly to voxels x time — no transpose needed
            M <- matrix(x@.Data, nrow = nv, ncol = nt)
            if (center) {
              M <- M - rowMeans(M)
            }
            if (scale) {
              rsd <- sqrt(rowSums(M * M) / (nt - 1L))
              rsd[rsd == 0] <- 1
              M <- M / rsd
            }
            dim(M) <- d
            new("DenseNeuroVec", .Data = M, space = space(x),
                label = x@label, volume_labels = volume_labels(x))
          })

#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="SparseNeuroVec", center="logical", scale="logical"),
          def=function(x, center, scale) {
            # x@data is T x K (time x masked_voxels)
            M <- x@data
            if (center) {
              means <- colMeans(M)
              M <- sweep(M, 2, means, "-")
            }
            if (scale) {
              if (nrow(M) <= 1L) {
                sds <- rep(1, ncol(M))
              } else if (center) {
                sds <- sqrt(colSums(M * M) / (nrow(M) - 1L))
              } else {
                means <- colMeans(M)
                centered <- M - rep(means, each = nrow(M))
                sds <- sqrt(colSums(centered * centered) / (nrow(M) - 1L))
              }
              sds[!is.finite(sds) | sds == 0] <- 1
              M <- sweep(M, 2, sds, "/")
            }
            M <- unname(as.matrix(M))
            SparseNeuroVec(M, space(x), mask(x))
          })

#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="logical", scale="logical"),
          def=function(x, center, scale) {
            M <- as.matrix(x)
            Ms <- base::scale(t(M), center, scale)
            NeuroVec(Ms, space(x))

          })


#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="missing", scale="logical"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, scale)
          })



#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="missing", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, TRUE)
          })


#' .isExtension
#' @keywords internal
#' @noRd
.isExtension <- function(fname, extension) {
  last <- substr(fname, nchar(fname)+1 - nchar(extension), nchar(fname))
  return(last==extension)
}

#' .concat4D
#' @keywords internal
#' @noRd
.concat4D <- function(x, y, ...) {
  rest <- list(...)
  objs <- c(list(x, y), rest)

  D <- dim(x)[1:3]

  rvols <- lapply(rest, function(z) {
    stopifnot(length(dim(z)) >= 3)
    stopifnot(identical(D, dim(z)[1:3]))
    as(z, "matrix")
  })

  clist <- if (length(rvols) > 0) {
    c(list(as(x, "matrix"), as(y, "matrix")), rvols)
  } else {
    list(as(x, "matrix"), as(y, "matrix"))
  }


  out <- do.call(cbind, clist)

  nvols <- ncol(out)
  new.dim <- c(D, nvols)

  nspace <-
    NeuroSpace(
      new.dim,
      origin = origin(x@space),
      spacing = spacing(x@space),
      axes = axes(x@space),
      trans = trans(x@space)
    )

  volume_labels <- {
    labs <- unlist(lapply(objs, function(obj) {
      if (inherits(obj, "NeuroVec")) {
        cur <- volume_labels(obj)
        if (length(cur) == 0L) rep("", dim(obj)[4]) else cur
      } else if (inherits(obj, "NeuroVol")) {
        ""
      } else {
        character()
      }
    }), use.names = FALSE)

    if (any(nzchar(labs))) labs else character()
  }

  DenseNeuroVec(out, nspace, volume_labels = volume_labels)

}




#' .gridToIndex3D
#' @keywords internal
#' @noRd
.gridToIndex3D <- function(dimensions, voxmat) {
  if (length(dimensions) != 3) {
    cli::cli_abort("{.arg dimensions} must have length 3, not {length(dimensions)}.")
  }
  if (is.vector(voxmat)) {
    if (length(voxmat) != 3) {
      cli::cli_abort("{.arg voxmat} vector must have length 3, not {length(voxmat)}.")
    }
    voxmat <- matrix(voxmat, 1,3)
  }

  if (ncol(voxmat) != 3) {
    cli::cli_abort("{.arg voxmat} matrix must have 3 columns, not {ncol(voxmat)}.")
  }
  gridToIndex3DCpp(dimensions, voxmat)

}

#' .gridToIndex
#' @keywords internal
#' @noRd
.gridToIndex <- function(dimensions, vmat) {
  vmat <- as.matrix(vmat)
  if (length(dimensions) != ncol(vmat)) {
    cli::cli_abort("length(dimensions) not equal to ncol(vmat): {length(dimensions)} != {ncol(vmat)}.")
  }

  gridToIndexCpp(as.integer(dimensions), vmat)
}

#' .indexToGrid
#' @keywords internal
#' @noRd
.indexToGrid <- function(idx, array.dim) {
  if (!all(idx > 0 & idx <= prod(array.dim))) {
    cli::cli_abort("{.arg idx} contains out-of-bounds values (must be in [1, {prod(array.dim)}]).")
  }
  if (length(array.dim) > 5) {
    cli::cli_abort("{.arg array.dim} must have length <= 5, not {length(array.dim)}.")
  }
  indexToGridCpp(idx, array.dim)

}



#' .getRStorage
#' @keywords internal
#' @noRd
.getRStorage <- function(data_type) {
  dtype_upper <- toupper(data_type)
  if (any(dtype_upper == c("BINARY", "BYTE", "UBYTE", "SHORT", "INTEGER", "INT", "LONG"))) {
    "integer"
  } else if (any(dtype_upper == c("FLOAT", "DOUBLE"))) {
    "double"
  } else {
	  stop(paste("unrecognized data type", data_type))
  }
}


#' @noRd
.isSigned <- function(data_type) {
  if (data_type == "UBYTE") {
    FALSE
  } else {
    TRUE
  }
}



#' @keywords internal
#' @importFrom mmap int8 uint8 int16 int32 real32 real64
#' @importFrom mmap mmap char mmapFlags munmap
#' @noRd
.getMMapMode <- function(code) {
	if (code == "UNKNOWN") {
		stop(paste(".getMMapMode: no memory map mode for UNKNOWN data type: ", code))
	} else if (code == "BINARY") {
		mmap::int8()
	} else if (code == "UBYTE") {
	  mmap::uint8()
	} else if(code == "SHORT") {
	  mmap::int16()
	} else if(code == "INT") {
	  mmap::int32()
	} else if (code == "FLOAT") {
	  mmap::real32()
	} else if (code == "DOUBLE") {
	  mmap::real64()
	} else {
		stop(paste(".getMMapMode: unsupported data type: ", code))
	}
}


# ---- Data-type lookup tables ------------------------------------------------

#' Named vectors mapping between NIfTI data-type codes, names, and byte sizes.
#' @keywords internal
#' @noRd
.DATA_CODE_TO_STORAGE <- c(
  "0"  = "UNKNOWN", "1"  = "BINARY", "2"  = "UBYTE",
  "4"  = "SHORT",   "8"  = "INT",    "16" = "FLOAT",
  "64" = "DOUBLE"
)

.DATA_STORAGE_TO_CODE <- c(
  UNKNOWN = 0L, BINARY = 1L, UBYTE = 2L, SHORT = 4L,
  INT = 8L,     FLOAT = 16L, DOUBLE = 64L
)

.DATA_TYPE_SIZE <- c(
  BINARY = 1L, BYTE = 1L, UBYTE = 1L, SHORT = 2L,
  INTEGER = 4L, INT = 4L, FLOAT = 4L, DOUBLE = 8L, LONG = 8L
)

#' .getDataStorage
#' @keywords internal
#' @noRd
.getDataStorage <- function(code) {
  key <- as.character(code)
  res <- .DATA_CODE_TO_STORAGE[key]
  if (is.na(res)) {
    cli::cli_abort("Unsupported NIfTI data-type code: {.val {code}}.")
  }
  unname(res)
}

#' .getDataCode
#' @keywords internal
#' @noRd
.getDataCode <- function(data_type) {
  res <- .DATA_STORAGE_TO_CODE[data_type]
  if (is.na(res)) {
    cli::cli_abort("Unsupported data type: {.val {data_type}}.")
  }
  unname(res)
}

#' .getDataSize
#' @keywords internal
#' @noRd
.getDataSize <- function(data_type) {
  res <- .DATA_TYPE_SIZE[data_type]
  if (is.na(res)) {
    cli::cli_abort("Unrecognized data type: {.val {data_type}}.")
  }
  unname(res)
}

#' .getEndian
#' @keywords internal
#' @noRd
.getEndian <- function(conn) {
  #try little endian
  endian <- "little"

  hsize <- readBin(conn, integer(), 1, endian=endian)
  if (hsize != 348) {
    # might be bigendian
    endian <- "big"
    seek(conn, 0)
    hsize <- readBin(conn, integer(), 1, endian=endian)
    if (hsize != 348) {
      stop("nifti(getEndian): header size is not 348, invalid header.")
    }
  }

  return(endian)
}



#' @keywords internal
#' @noRd
.niftiExt <- function(filetype) {

  extensions <- list()

  if (filetype == "nifti-single") {
    extensions[["header"]]  <- "nii"
    extensions[["data"]] <- "nii"
  }
  else if (filetype == "nifti-pair") {
    extensions[["header"]]  <- "hdr"
    extensions[["data"]] <- "img"
  }
  else if (filetype == "nifti-gz") {
    extensions[["header"]]  <- "nii.gz"
    extensions[["data"]] <- "nii.gz"
  } else {
    stop(paste("unsupported filetype: ", filetype))
  }

  return(extensions)
}

#' Convert a Transformation Matrix to a Quaternion Representation
#'
#' @description
#' Extracts the rotation and scaling components from a 3x3 (or 4x4) transformation
#' matrix, normalizes them, and computes the corresponding quaternion parameters
#' and a sign factor (`qfac`) indicating whether the determinant is negative.
#'
#' @details
#' This function first checks and corrects for zero-length axes in the upper-left
#' corner of the matrix, then normalizes each column to extract the pure rotation.
#' If the determinant of the rotation submatrix is negative, the \code{qfac} is set
#' to \code{-1}, and the third column is negated. Finally, the quaternion parameters
#' (\eqn{a, b, c, d}) are computed following standard NIfTI-1 conventions for
#' representing the rotation in 3D.
#'
#' @param mat A numeric matrix with at least the top-left 3x3 portion containing
#'   rotation/scaling. Often a 4x4 affine transform, but only the 3x3 top-left
#'   submatrix is used in practice.
#'
#' @return A named \code{list} with two elements:
#'   \describe{
#'     \item{\code{quaternion}}{A numeric vector of length 3, \eqn{(b, c, d)},
#'       which—together with \eqn{a} derived internally—represents the rotation.}
#'     \item{\code{qfac}}{Either \code{+1} or \code{-1}, indicating whether the
#'       determinant of the rotation submatrix is positive or negative, respectively.}
#'   }
#'
#' @seealso
#' \code{\link{quaternToMatrix}} for the inverse operation, converting
#' quaternion parameters back to a transform matrix.
#'
#' @references
#' - Cox RW. *Analysis of Functional NeuroImages* (AFNI) and NIfTI-1 quaternion
#'   conventions. \url{https://afni.nimh.nih.gov}
#'
#' @export
matrixToQuatern <- function(mat) {
  xd <- sqrt(drop(crossprod(mat[1:3,1])))
  yd <- sqrt(drop(crossprod(mat[1:3,2])))
  zd <- sqrt(drop(crossprod(mat[1:3,3])))

  if (xd == 0) { mat[1,1] = 1; mat[2:3,1] = 0; xd = 1; }
  if (yd == 0) { mat[2,2] = 1; mat[c(1,3),2] = 0; yd = 1; }
  if (zd == 0) { mat[3,3] = 1; mat[1:2,3] = 0; zd = 1; }

  rmat = mat[1:3, 1:3]
  rmat[,1] = rmat[,1]/xd
  rmat[,2] = rmat[,2]/yd
  rmat[,3] = rmat[,3]/zd

  ####### skipping orthogonalization of columns

  #################################################

  zd = det(rmat)
  qfac = 1

  if (zd > 0) {
    qfac = 1
  } else {
    qfac = -1
    rmat[1:3,3] = -rmat[1:3,3]
  }

  # compute quaternion parameters

  a = rmat[1,1] + rmat[2,2] + rmat[3,3] + 1

  if (a > .5) {
    a = .5 * sqrt(a)
    b = 0.25 * (rmat[3,2]-rmat[2,3]) / a
    c = 0.25 * (rmat[1,3]-rmat[3,1]) / a
    d = 0.25 * (rmat[2,1]-rmat[1,2]) / a
   } else {
     xd = 1.0 + rmat[1,1] - (rmat[2,2]+rmat[3,3])
     yd = 1.0 + rmat[2,2] - (rmat[1,1]+rmat[3,3])
     zd = 1.0 + rmat[3,3] - (rmat[1,1]+rmat[2,2])
     if( xd > 1.0 ){
       b = 0.5 * sqrt(xd)
       c = 0.25* (rmat[1,2]+rmat[2,1]) / b
       d = 0.25* (rmat[1,3]+rmat[3,1]) / b
       a = 0.25* (rmat[3,2]-rmat[2,3]) / b
     } else if( yd > 1.0 ){
       c = 0.5 * sqrt(yd) ;
       b = 0.25* (rmat[1,2]+rmat[2,1]) / c
       d = 0.25* (rmat[2,3]+rmat[3,2]) / c
       a = 0.25* (rmat[1,3]-rmat[3,1]) / c
     } else {
       d = 0.5 * sqrt(zd) ;
       b = 0.25* (rmat[1,3]+rmat[3,1]) / d
       c = 0.25* (rmat[2,3]+rmat[3,2]) / d
       a = 0.25* (rmat[2,1]-rmat[1,2]) / d
     }
     if( a < 0.0 ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

  return(list(quaternion=c(b,c,d), qfac=qfac))

}


#' Convert Quaternion Parameters to a Transformation Matrix
#'
#' @description
#' Given a quaternion \code{(b, c, d)}, a scalar offset (origin), voxel step sizes,
#' and the \code{qfac} sign, reconstructs a 4x4 affine matrix representing rotation,
#' scaling, and translation as used in NIfTI-1.
#'
#' @details
#' This function uses the quaternion formalism common in neuroimaging, adding the
#' offset (translation) into the 4th column, and applying the voxel sizes along
#' each axis. If \code{qfac} is \code{-1}, the \eqn{z} scale is negated. The
#' resulting 4x4 matrix is typically used as an affine transform for voxel-to-world
#' coordinate mapping.
#'
#' @param quat A numeric vector of length 3 containing the quaternion parameters
#'   \eqn{(b, c, d)}. The scalar part \eqn{a} is computed internally.
#' @param origin A numeric vector of length 3 specifying the translation components
#'   (often the real-space origin or offset).
#' @param stepSize A numeric vector of length 3 giving the voxel dimensions along
#'   each axis (e.g., \code{(dx, dy, dz)}).
#' @param qfac Either \code{+1} or \code{-1}, indicating the sign from the
#'   determinant check in \code{\link{matrixToQuatern}}.
#'
#' @return A 4x4 numeric affine transformation matrix. The top-left 3x3 submatrix
#'   encodes rotation and scaling, and the 4th column encodes translation.
#'
#' @seealso
#' \code{\link{matrixToQuatern}} for converting a matrix back to quaternion form.
#'
#' @export
quaternToMatrix <- function(quat, origin, stepSize, qfac) {
  mat <- matrix(0, 4,4)
  mat[4,] <- c(0,0,0,1)

  a <- 1 - sum(quat^2)
  if (a < 1e-07) {
    a <- 1 /(sqrt(sum(quat^2)))
    quat <- quat*a
    a <- 0
  } else {
    a <- sqrt(a)
  }

  stepSize <- ifelse(stepSize > 0, stepSize, 1)
  xd <- stepSize[1]
  yd <- stepSize[2]
  zd <- stepSize[3]

  if (qfac < 0) {
    zd <- -zd
  }

  b <- quat[1]
  c <- quat[2]
  d <- quat[3]


  mat[1,1] <- (a*a+b*b-c*c-d*d) * xd
  mat[1,2] <- 2 * (b*c-a*d)     * yd
  mat[1,3] <- 2 * (b*d+a*c)     * zd
  mat[2,1] <- 2 * (b*c+a*d)     * xd
  mat[2,2] <- (a*a+c*c-b*b-d*d) * yd
  mat[2,3] <- 2 * (c*d-a*b)     * zd
  mat[3,1] <- 2 * (b*d-a*c)     * xd
  mat[3,2] <- 2 * (c*d+a*b)     * yd
  mat[3,3] <- (a*a+d*d-c*c-b*b) * zd

  mat[1:3,4] <- origin

  return(mat)
}

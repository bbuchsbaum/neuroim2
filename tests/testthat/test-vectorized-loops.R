library(testthat)
library(neuroim2)

# Phase 4 replaced five scalar R loops with vectorised equivalents. Each test
# below pins the rewrite against a naive reference that mirrors the loop it
# replaced, so the reference -- not the implementation -- defines correctness.

# ---------------------------------------------------------------- ClusteredNeuroVec

make_clustered <- function(d, ncl, holey = FALSE, nt = 11L, seed = 1L) {
  set.seed(seed)
  nvox <- prod(d)
  cl <- sample.int(ncl, nvox, TRUE)
  if (holey) {
    cl[sample.int(nvox, nvox %/% 3L)] <- 0L
    mk <- array(cl > 0, d)
    cvol <- suppressWarnings(
      ClusteredNeuroVol(LogicalNeuroVol(mk, NeuroSpace(d)), cl[cl > 0]))
  } else {
    cvol <- suppressWarnings(
      ClusteredNeuroVol(LogicalNeuroVol(array(TRUE, d), NeuroSpace(d)), cl))
  }
  vec <- NeuroVec(array(rnorm(nvox * nt), c(d, nt)), NeuroSpace(c(d, nt)))
  ClusteredNeuroVec(vec, cvol)
}

# The loop this replaced, verbatim in behaviour.
clustered_ref <- function(x, i, j, k, m, drop) {
  dims <- dim(space(x@cvol))
  idx <- ((k - 1) * dims[1] * dims[2]) + ((j - 1) * dims[1]) + i
  cids <- x@cl_map[idx]
  result <- matrix(NA_real_, length(m), length(idx))
  for (t_idx in seq_along(m)) {
    for (v_idx in seq_along(idx)) {
      if (cids[v_idx] > 0) result[t_idx, v_idx] <- x@ts[m[t_idx], cids[v_idx]]
    }
  }
  if (drop) base::drop(result) else result
}

test_that("ClusteredNeuroVec [ matches the scalar reference", {
  for (holey in c(FALSE, TRUE)) {
    for (ncl in c(3L, 25L)) {
      x <- make_clustered(c(7L, 6L, 5L), ncl, holey = holey)
      ntp <- nrow(x@ts)
      set.seed(2)
      for (nq in c(1L, 3L, 25L)) {
        i <- sample.int(7L, nq, TRUE); j <- sample.int(6L, nq, TRUE)
        k <- sample.int(5L, nq, TRUE)
        # single, several, all, reversed, and repeated timepoints
        for (m in list(1L, c(2L, 4L), seq_len(ntp), rev(seq_len(ntp)), c(3L, 3L, 1L))) {
          for (dr in c(TRUE, FALSE)) {
            lbl <- sprintf("holey=%s ncl=%d nq=%d nm=%d drop=%s",
                           holey, ncl, nq, length(m), dr)
            expect_equal(x[i, j, k, m, drop = dr],
                         clustered_ref(x, i, j, k, m, dr), info = lbl)
          }
        }
      }
    }
  }
})

test_that("ClusteredNeuroVec [ returns NA for voxels outside any cluster", {
  x <- make_clustered(c(6L, 6L, 6L), 5L, holey = TRUE)
  outside <- which(x@cl_map == 0)
  # The fixture masks out a third of the volume, so this must exist; skipping
  # here would let the assertion below quietly stop running.
  expect_gt(length(outside), 0)
  g <- arrayInd(outside[1], c(6L, 6L, 6L))
  got <- x[g[1], g[2], g[3], 1:3, drop = FALSE]
  expect_true(all(is.na(got)))
  expect_equal(dim(got), c(3L, 1L))
})

test_that("ClusteredNeuroVec [ rejects coordinates outside the volume", {
  x <- make_clustered(c(6L, 6L, 6L), 5L)
  # (7, 1, 1) on a 6-wide axis used to reach linear index 43 -- in range -- and
  # silently return voxel (1, 3, 2).
  expect_error(x[7, 1, 1, 1], "outside the volume")
  expect_error(x[1, 0, 1, 1], "outside the volume")
  expect_error(x[1, 1, -1, 1], "outside the volume")
  expect_error(x[1, 1, 1, 1], NA)   # in range: no error

  # Fractional coordinates truncate as R's own array indexing does, and only
  # the truncated value is bounds-checked -- so 2.7 and 6.5 behave alike, rather
  # than one being accepted and the other rejected.
  expect_equal(x[2.7, 2, 2, 1:2, drop = FALSE], x[2, 2, 2, 1:2, drop = FALSE])
  expect_equal(x[6.5, 2, 2, 1:2, drop = FALSE], x[6, 2, 2, 1:2, drop = FALSE])
  expect_error(x[0.5, 2, 2, 1], "outside the volume")
})

test_that("ClusteredNeuroVec [ rejects timepoints it cannot index", {
  x <- make_clustered(c(6L, 6L, 6L), 5L)
  nt <- nrow(x@ts)
  # A negative m would be read as R's drop-these-rows indexing by the matrix
  # subset; the scalar loop this replaced could not express that.
  expect_error(x[1, 1, 1, c(-1, -2)], "must index timepoints")
  expect_error(x[1, 1, 1, 0], "must index timepoints")
  expect_error(x[1, 1, 1, nt + 1L], "must index timepoints")
  expect_error(x[1, 1, 1, NA_integer_], "must index timepoints")
  expect_equal(dim(x[1, 1, 1, integer(0), drop = FALSE]), c(0L, 1L))
})

# ---------------------------------------------------------------------------- mapf

# The convolution the loop expressed, with out-of-volume taps skipped.
mapf_ref <- function(v, kern, mask = NULL) {
  d <- dim(v)
  arr <- as.array(v)
  hw <- vapply(kern@width, function(w) ceiling(w / 2 - 1), numeric(1)) + 1
  grid <- if (is.null(mask)) {
    as.matrix(expand.grid(i = hw[1]:(d[1] - hw[1]), j = hw[2]:(d[2] - hw[2]),
                          k = hw[3]:(d[3] - hw[3])))
  } else {
    index_to_grid(mask, which(mask != 0))
  }
  out <- array(0, d)
  for (r in seq_len(nrow(grid))) {
    ii <- grid[r, 1] + kern@voxels[, 1]
    jj <- grid[r, 2] + kern@voxels[, 2]
    kk <- grid[r, 3] + kern@voxels[, 3]
    ok <- ii >= 1 & ii <= d[1] & jj >= 1 & jj <= d[2] & kk >= 1 & kk <= d[3]
    out[grid[r, 1], grid[r, 2], grid[r, 3]] <-
      sum(arr[cbind(ii[ok], jj[ok], kk[ok])] * kern@weights[ok])
  }
  out
}

test_that("mapf matches the scalar reference, unmasked and masked", {
  set.seed(3)
  for (d in list(c(9L, 8L, 10L), c(12L, 12L, 12L))) {
    v <- NeuroVol(array(rnorm(prod(d)), d), NeuroSpace(d))
    for (kw in list(c(3, 3, 3), c(5, 5, 5), c(3, 5, 3))) {
      kern <- Kernel(kw, c(1, 1, 1))
      lbl <- sprintf("dim=%s kernel=%s", paste(d, collapse = "x"),
                     paste(kw, collapse = "x"))
      expect_equal(as.numeric(mapf(v, kern)),
                   as.numeric(mapf_ref(v, kern)), info = lbl)

      # A mask that deliberately includes both faces of the volume: the old
      # per-centre gather could not handle either. A coordinate past the far
      # face raised "subscript out of bounds"; one before the near face made
      # the gather return a shorter vector that recycled against the weights.
      mk <- array(FALSE, d)
      mk[c(1L, d[1]), , ] <- TRUE
      mk[, c(1L, d[2]), ] <- TRUE
      mk[, , c(1L, d[3])] <- TRUE
      mk[d[1] %/% 2L, d[2] %/% 2L, d[3] %/% 2L] <- TRUE
      mv <- LogicalNeuroVol(mk, NeuroSpace(d))
      expect_equal(as.numeric(mapf(v, kern, mv)),
                   as.numeric(mapf_ref(v, kern, mv)),
                   info = paste("boundary mask,", lbl))
      expect_false(anyNA(as.numeric(mapf(v, kern, mv))))
    }
  }
})

test_that("mapf clips taps on a volume thinner than the kernel", {
  # `hwidth[d]:(dims[d] - hwidth[d])` counts *backwards* when the volume is
  # thinner than twice the kernel half-width, so the unmasked grid contains
  # centres inside the margin and taps really do leave the volume -- the one
  # unmasked case where the boundary behaviour changed. On a 10 x 3 x 10 slab
  # with a 3x3x3 kernel the j range is 2:1.
  set.seed(31)
  d <- c(10L, 3L, 10L)
  v <- NeuroVol(array(rnorm(prod(d)), d), NeuroSpace(d))
  kern <- Kernel(c(3, 3, 3), c(1, 1, 1))
  hw <- vapply(kern@width, function(w) ceiling(w / 2 - 1), numeric(1)) + 1
  expect_gt(hw[2], d[2] - hw[2])              # the range really does reverse
  got <- mapf(v, kern)
  expect_equal(as.numeric(got), as.numeric(mapf_ref(v, kern)))
  expect_false(anyNA(as.numeric(got)))
  # The affected centres get a genuinely clipped convolution: at j = 1 only
  # 18 of the 27 taps are inside the volume.
  arr <- as.array(v)
  ii <- 5L + kern@voxels[, 1]; jj <- 1L + kern@voxels[, 2]
  kk <- 5L + kern@voxels[, 3]
  ok <- ii >= 1 & ii <= d[1] & jj >= 1 & jj <= d[2] & kk >= 1 & kk <= d[3]
  expect_equal(sum(ok), 18L)
  expect_equal(as.array(got)[5, 1, 5],
               sum(arr[cbind(ii[ok], jj[ok], kk[ok])] * kern@weights[ok]))
})

test_that("mapf handles the NeuroVol subclasses the old gather did", {
  set.seed(32)
  d <- c(9L, 9L, 9L)
  kern <- Kernel(c(3, 3, 3), c(1, 1, 1))
  arr <- array(rnorm(prod(d)), d)
  dense <- NeuroVol(arr, NeuroSpace(d))
  logical_vol <- LogicalNeuroVol(array(arr > 0, d), NeuroSpace(d))
  idx <- which(arr > 0.5)
  sparse <- SparseNeuroVol(arr[idx], NeuroSpace(d), indices = idx)
  for (v in list(dense, logical_vol, sparse)) {
    expect_equal(as.numeric(mapf(v, kern)), as.numeric(mapf_ref(v, kern)),
                 info = class(v)[1])
  }
})

test_that("mapf returns an all-zero volume for an empty mask", {
  d <- c(6L, 6L, 6L)
  v <- NeuroVol(array(rnorm(prod(d)), d), NeuroSpace(d))
  empty <- LogicalNeuroVol(array(FALSE, d), NeuroSpace(d))
  # min()/max() over the empty centre grid must not warn.
  expect_silent(out <- mapf(v, Kernel(c(3, 3, 3), c(1, 1, 1)), empty))
  expect_equal(dim(out), d)
  expect_true(all(as.numeric(out) == 0))
})

test_that("mapf leaves the interior untouched by the boundary change", {
  # Unmasked centres are held clear of the border, so no tap is ever clipped
  # and the vectorised result must equal a reference that does no clipping.
  set.seed(4)
  d <- c(11L, 10L, 9L)
  v <- NeuroVol(array(rnorm(prod(d)), d), NeuroSpace(d))
  kern <- Kernel(c(3, 3, 3), c(1, 1, 1))
  arr <- as.array(v)
  hw <- vapply(kern@width, function(w) ceiling(w / 2 - 1), numeric(1)) + 1
  grid <- as.matrix(expand.grid(i = hw[1]:(d[1] - hw[1]), j = hw[2]:(d[2] - hw[2]),
                                k = hw[3]:(d[3] - hw[3])))
  ref <- array(0, d)
  for (r in seq_len(nrow(grid))) {
    ii <- grid[r, 1] + kern@voxels[, 1]; jj <- grid[r, 2] + kern@voxels[, 2]
    kk <- grid[r, 3] + kern@voxels[, 3]
    # no clipping at all -- would error if any tap left the volume
    ref[grid[r, 1], grid[r, 2], grid[r, 3]] <-
      sum(arr[cbind(ii, jj, kk)] * kern@weights)
  }
  expect_equal(as.numeric(mapf(v, kern)), as.numeric(ref))
})

test_that("mapf reports a mask of the wrong shape", {
  d <- c(8L, 8L, 8L)
  v <- NeuroVol(array(0, d), NeuroSpace(d))
  bad <- LogicalNeuroVol(array(TRUE, c(4L, 4L, 4L)), NeuroSpace(c(4L, 4L, 4L)))
  # The guard was `!all.equal(...)`, which throws "invalid argument type" on a
  # mismatch because all.equal returns a character description, not FALSE.
  expect_error(mapf(v, Kernel(c(3, 3, 3), c(1, 1, 1)), bad),
               "same dimensions")
})

# ------------------------------------------------------------------ NeuroHyperVec

make_hyper <- function(d, frac, ntr, nf, seed = 1L) {
  set.seed(seed)
  mk <- array(runif(prod(d)) < frac, d)
  mk[1] <- TRUE
  nvm <- sum(mk)
  dat <- array(rnorm(nf * ntr * nvm), c(nf, ntr, nvm))
  NeuroHyperVec(dat, NeuroSpace(c(d, ntr, nf)),
                LogicalNeuroVol(mk, NeuroSpace(d)))
}

dense_ref <- function(x) {
  dims <- dim(x); sd3 <- dims[1:3]; ntr <- dims[4]; nf <- dims[5]
  dense <- array(0, dim = dims)
  mask_idx <- which(x@mask@.Data)
  for (f in seq_len(nf)) for (t in seq_len(ntr)) {
    vol <- numeric(prod(sd3))
    vol[mask_idx] <- x@data[f, t, ]
    dense[, , , t, f] <- array(vol, dim = sd3)
  }
  dense
}

test_that("dense_array_5d matches the per-(feature, trial) reference", {
  for (d in list(c(5L, 6L, 4L), c(8L, 8L, 8L))) {
    for (frac in c(0.25, 1.0)) for (ntr in c(1L, 3L)) for (nf in c(1L, 4L)) {
      hv <- make_hyper(d, frac, ntr, nf)
      expect_identical(neuroim2:::dense_array_5d(hv), dense_ref(hv),
                       info = sprintf("dim=%s frac=%g trials=%d features=%d",
                                      paste(d, collapse = "x"), frac, ntr, nf))
    }
  }
})

hyper_idx_ref <- function(x, i, j, k, l, m) {
  sd3 <- dim(x)[1:3]
  g <- expand.grid(i = i, j = j, k = k)
  lin <- (g$k - 1) * (sd3[1] * sd3[2]) + (g$j - 1) * sd3[1] + g$i
  lu <- x@lookup_map[lin]
  out <- array(0, c(length(i), length(j), length(k), length(l), length(m)))
  nz <- lu > 0
  if (any(nz)) {
    pos <- which(nz)
    ai <- arrayInd(pos, .dim = c(length(i), length(j), length(k)))
    vi <- lu[nz]
    for (idx in seq_along(vi)) {
      vd <- aperm(x@data[m, l, vi[idx], drop = FALSE], c(2, 1, 3))[, , 1, drop = FALSE]
      out[ai[idx, 1], ai[idx, 2], ai[idx, 3], , ] <- vd
    }
  }
  out
}

test_that("NeuroHyperVec [ matches the per-voxel reference", {
  for (d in list(c(6L, 5L, 7L), c(9L, 9L, 9L))) {
    for (frac in c(0.3, 1.0)) for (ntr in c(1L, 4L)) for (nf in c(1L, 3L)) {
      hv <- make_hyper(d, frac, ntr, nf)
      for (nq in c(1L, 2L, 5L)) {
        ii <- seq_len(min(nq, d[1])); jj <- seq_len(min(nq, d[2]))
        kk <- seq_len(min(nq, d[3]))
        # forward, reversed and repeated trial/feature selectors
        for (sel in list(list(seq_len(ntr), seq_len(nf)),
                         list(rev(seq_len(ntr)), rev(seq_len(nf))),
                         list(rep(1L, ntr), rep(1L, nf)))) {
          lbl <- sprintf("dim=%s frac=%g t=%d f=%d nq=%d",
                         paste(d, collapse = "x"), frac, ntr, nf, nq)
          expect_equal(hv[ii, jj, kk, sel[[1]], sel[[2]], drop = FALSE],
                       hyper_idx_ref(hv, ii, jj, kk, sel[[1]], sel[[2]]),
                       info = lbl)
        }
      }
    }
  }
})

test_that("NeuroHyperVec [ leaves out-of-mask positions at zero", {
  d <- c(6L, 6L, 6L)
  set.seed(9)
  mk <- array(FALSE, d); mk[2, 2, 2] <- TRUE
  dat <- array(rnorm(2 * 3 * 1), c(2, 3, 1))
  hv <- NeuroHyperVec(dat, NeuroSpace(c(d, 3L, 2L)),
                      LogicalNeuroVol(mk, NeuroSpace(d)))
  got <- hv[1:3, 1:3, 1:3, 1:3, 1:2, drop = FALSE]
  expect_equal(dim(got), c(3L, 3L, 3L, 3L, 2L))
  expect_equal(got[2, 2, 2, , ], t(dat[, , 1]))
  got[2, 2, 2, , ] <- 0
  expect_true(all(got == 0))
})

# ------------------------------------------------------------- random_searchlight

test_that("random_searchlight preserves its structural invariants", {
  for (d in list(c(9L, 9L, 9L), c(12L, 10L, 8L))) {
    for (frac in c(0.4, 1.0)) {
      set.seed(11)
      mk <- array(runif(prod(d)) < frac, d)
      if (sum(mk) < 10) next
      mask <- LogicalNeuroVol(mk, NeuroSpace(d, c(2, 2, 2)))
      midx <- sort(which(mk))
      for (r in c(2.5, 5)) for (nz in c(TRUE, FALSE)) {
        set.seed(23)
        sl <- random_searchlight(mask, radius = r, nonzero = nz)
        lbl <- sprintf("dim=%s frac=%g r=%g nonzero=%s",
                       paste(d, collapse = "x"), frac, r, nz)
        expect_gt(length(sl), 0)
        centres <- vapply(sl, function(s) attr(s, "mask_index"), integer(1))
        # A voxel is used as a centre at most once ...
        expect_false(anyDuplicated(centres) > 0, label = paste("duplicate centre:", lbl))
        expect_true(all(centres >= 1L & centres <= length(midx)), label = lbl)
        lin <- unlist(lapply(sl, function(s) {
          cc <- coords(s)
          cc[, 1] + (cc[, 2] - 1) * d[1] + (cc[, 3] - 1) * d[1] * d[2]
        }))
        # ... and a masked voxel is claimed by at most one searchlight.
        inmask <- lin[lin %in% midx]
        expect_false(anyDuplicated(inmask) > 0, label = paste("double claim:", lbl))
        if (nz) {
          expect_true(all(lin %in% midx), label = paste("stayed in mask:", lbl))
        }
        expect_true(all(vapply(sl, function(s) is(s, "ROIVolWindow"), logical(1))),
                    label = lbl)
      }
    }
  }
})

test_that("random_searchlight claims every voxel a full-mask sweep can reach", {
  # With nonzero = TRUE and a fully populated mask, no voxel can be dropped
  # through the "no available neighbours" path before being claimed, so the
  # searchlights partition the mask exactly. This is what the free-list rewrite
  # had to preserve.
  d <- c(10L, 10L, 10L)
  mask <- LogicalNeuroVol(array(TRUE, d), NeuroSpace(d, c(2, 2, 2)))
  for (seed in 1:5) {
    set.seed(seed)
    sl <- random_searchlight(mask, radius = 4, nonzero = TRUE)
    lin <- unlist(lapply(sl, function(s) {
      cc <- coords(s)
      cc[, 1] + (cc[, 2] - 1) * d[1] + (cc[, 3] - 1) * d[1] * d[2]
    }))
    expect_equal(sort(lin), seq_len(prod(d)), label = sprintf("seed=%d", seed))
  }
})

test_that("random_searchlight depends on the seed and repeats under it", {
  d <- c(8L, 8L, 8L)
  mask <- LogicalNeuroVol(array(TRUE, d), NeuroSpace(d, c(2, 2, 2)))
  set.seed(4); a <- random_searchlight(mask, radius = 3)
  set.seed(4); b <- random_searchlight(mask, radius = 3)
  expect_identical(lapply(a, coords), lapply(b, coords))
  # ... and that it is genuinely random: a different seed must give a different
  # partition. Without this the test above holds for any implementation that
  # merely draws from R's RNG, including a broken one.
  set.seed(5); c5 <- random_searchlight(mask, radius = 3)
  expect_false(identical(lapply(a, coords), lapply(c5, coords)))
  # The first centre is drawn before any removal, so it cannot depend on the
  # free list's ordering -- which is what keeps existing seeded tests valid.
  set.seed(6); p <- random_searchlight(mask, radius = 3)
  set.seed(6); q <- random_searchlight(mask, radius = 3, nonzero = FALSE)
  expect_identical(attr(p[[1]], "mask_index"), attr(q[[1]], "mask_index"))
})

# ---------------------------------------------------------------------- dead code

test_that("unreferenced compiled entry points are gone", {
  ns <- asNamespace("neuroim2")
  for (nm in c("radius_search_3d_nonisotropic", "radius_search_3d_direct",
               "radius_search_3d_precomputed", "local_spheres",
               "kernel_filt_3d_cpp")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE),
                 label = sprintf("%s should have been removed", nm))
  }
  # local_sphere stays: it is the reference the offset-template equivalence
  # test in test-roi-series-fastpaths.R compares against.
  expect_true(exists("local_sphere", envir = ns, inherits = FALSE))

  # The R wrappers going away is not the same as the native symbols going away,
  # so check the registration table too.
  routines <- names(getDLLRegisteredRoutines("neuroim2")$.Call)
  for (nm in c("_neuroim2_radius_search_3d_nonisotropic",
               "_neuroim2_radius_search_3d_direct",
               "_neuroim2_radius_search_3d_precomputed",
               "_neuroim2_local_spheres", "_neuroim2_kernel_filt_3d_cpp")) {
    expect_false(nm %in% routines, label = sprintf("%s still registered", nm))
  }
  expect_true("_neuroim2_local_sphere" %in% routines)
})

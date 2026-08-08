library(testthat)
library(neuroim2)

context("fast paths preserve behaviour")

# These tests pin the behaviour of three optimisations that must be
# indistinguishable from the implementations they replaced:
#   * .new_roi_vol_window() vs new("ROIVolWindow", ...)
#   * the cached sphere-offset template vs local_sphere()
#   * the compiled series() gathers vs the previous R loops
# They are behavioural equivalence tests, not performance tests.

# ---------------------------------------------------------------- constructor

test_that(".new_roi_vol_window is identical to new('ROIVolWindow', ...)", {
  sp <- NeuroSpace(c(20L, 20L, 20L), c(2, 2, 2))
  set.seed(42)

  for (n in c(0L, 1L, 5L, 257L)) {
    coords <- matrix(as.numeric(sample.int(20, max(n, 1L) * 3, TRUE)),
                     ncol = 3)[seq_len(n), , drop = FALSE]
    vals <- as.numeric(seq_len(n))

    a <- new("ROIVolWindow", vals, space = sp, coords = coords,
             center_index = 1L, parent_index = 5L)
    b <- neuroim2:::.new_roi_vol_window(vals, sp, coords, 1L, 5L)

    expect_identical(a, b, info = paste("n =", n))
    expect_true(validObject(b))
  }
})

test_that("the constructor coerces .Data the way new() does", {
  # new() coerces the data part to the basic type the class extends and drops
  # attributes; slot assignment does neither. Without .as_roi_data() a named or
  # logical `fill` produced an object that was NOT identical() to new()'s.
  sp <- NeuroSpace(c(10L, 10L, 10L))
  cds <- matrix(as.numeric(1:9), ncol = 3)

  for (d in list(c(a = 1, b = 2, c = 3), c(TRUE, FALSE, TRUE), 1:3,
                 c(1.5, 2.5, 3.5), list(1, 2, 3))) {
    a <- suppressWarnings(new("ROIVolWindow", d, space = sp, coords = cds,
                              center_index = 2L, parent_index = 5L))
    b <- suppressWarnings(neuroim2:::.new_roi_vol_window(d, sp, cds, 2L, 5L))
    expect_identical(a, b, info = paste(class(d), collapse = "/"))
  }

  # names must not leak through the public API either
  bvol <- NeuroVol(array(1, c(10, 10, 10)), NeuroSpace(c(10L, 10L, 10L)))
  r <- spherical_roi(bvol, c(5, 5, 5), radius = 2, fill = c(a = 1))
  expect_null(names(r@.Data))

  # non-double fills still work, as they did before the fast path
  expect_silent(rl <- spherical_roi(bvol, c(5, 5, 5), 2, fill = TRUE))
  expect_true(all(rl@.Data == 1))
  expect_identical(typeof(rl@.Data), "double")

  # factor fill: the old new()-based path was inconsistent about this
  # (sometimes an integer .Data carrying a stray "levels" attribute). Pin the
  # deterministic behaviour: integer codes, as a clean double, no extra attrs.
  rf <- suppressWarnings(spherical_roi(bvol, c(5, 5, 5), 2, fill = factor("x")))
  expect_identical(typeof(rf@.Data), "double")
  expect_true(all(rf@.Data == 1))
  expect_false("levels" %in% names(attributes(rf)))
})

test_that(".new_roi_vol_window enforces the class invariants", {
  sp <- NeuroSpace(c(10L, 10L, 10L))
  expect_error(neuroim2:::.new_roi_vol_window(c(1, 2), sp, matrix(1, 3, 3), 1L, 1L),
               "length of data vector")
  expect_error(neuroim2:::.new_roi_vol_window(c(1, 2), sp, matrix(1, 2, 2), 1L, 1L),
               "3 columns")
})

test_that("the cached ROI prototype is never mutated by construction", {
  sp <- NeuroSpace(c(10L, 10L, 10L))
  before <- neuroim2:::.roi_prototype("ROIVolWindow")

  for (i in 1:20) {
    obj <- neuroim2:::.new_roi_vol_window(rep(1, 4), sp, matrix(1, 4, 3), 1L, 2L)
    obj@.Data <- rep(99, 4)          # mutate the result
    attr(obj, "scribble") <- i
  }

  after <- neuroim2:::.roi_prototype("ROIVolWindow")
  expect_identical(before, after)
  expect_length(after@.Data, 0L)
  expect_equal(nrow(after@coords), 0L)
})

# ------------------------------------------------------------ sphere geometry

test_that("the offset template reproduces local_sphere() exactly", {
  set.seed(11)
  spacings <- list(c(1, 1, 1), c(2, 2, 2), c(1, 1, 4), c(0.8, 0.8, 3), c(3, 1, 2))
  dims <- list(c(12L, 12L, 12L), c(9L, 14L, 7L), c(6L, 20L, 11L))

  for (sp in spacings) {
    for (d in dims) {
      centroids <- rbind(c(1L, 1L, 1L), d, pmax(1L, d %/% 2L),
                         cbind(sample.int(d[1], 4, TRUE),
                               sample.int(d[2], 4, TRUE),
                               sample.int(d[3], 4, TRUE)))
      storage.mode(centroids) <- "integer"

      for (radius in c(min(sp), max(sp), 2.5, 4, 7)) {
        if (radius < min(sp)) next
        off <- neuroim2:::sphere_offsets_cpp(radius, sp)

        for (ci in seq_len(nrow(centroids))) {
          cent <- centroids[ci, ]
          ref <- neuroim2:::local_sphere(cent[1] - 1L, cent[2] - 1L, cent[3] - 1L,
                                         radius, sp, d)
          new <- neuroim2:::sphere_at_cpp(off, cent, d, base0 = TRUE)

          expect_equal(new, ref,
                       info = sprintf("spacing %s dim %s r %g centroid %s",
                                      paste(sp, collapse = "/"),
                                      paste(d, collapse = "x"), radius,
                                      paste(cent, collapse = ",")))
        }
      }
    }
  }
})

test_that("the compiled sphere path emits lexicographically sorted rows", {
  set.seed(12)
  for (sp in list(c(1, 1, 1), c(2, 2, 2), c(1, 1, 4), c(0.7, 1.3, 2.9))) {
    d <- c(15L, 13L, 11L)
    off <- neuroim2:::sphere_offsets_cpp(6, sp)
    for (i in 1:15) {
      cent <- c(sample.int(d[1], 1), sample.int(d[2], 1), sample.int(d[3], 1))
      storage.mode(cent) <- "integer"
      g <- neuroim2:::sphere_at_cpp(off, cent, d, base0 = FALSE)
      if (nrow(g) == 0) next
      expect_identical(order(g[, 1], g[, 2], g[, 3]), seq_len(nrow(g)))
    }
  }
})

test_that("spherical_roi keeps sorting the dbscan fallback path", {
  # The compiled path emits sorted rows natively, so spherical_roi() only
  # re-sorts when use_cpp = FALSE. That sort must still happen: dbscan::frNN
  # returns neighbours in its own order.
  sp1 <- NeuroSpace(c(12L, 12L, 12L), c(2, 2, 2))
  b <- spherical_roi(sp1, c(6L, 6L, 6L), 4, use_cpp = FALSE)
  g <- b@coords
  expect_identical(order(g[, 1], g[, 2], g[, 3]), seq_len(nrow(g)))

  # NOTE: the two paths select the same NUMBER of voxels but not the same
  # voxels -- a long-standing discrepancy in the legacy dbscan branch that is
  # unrelated to the offset template (both branches reproduce their own
  # pre-optimisation output exactly). Pinned so it stays visible; the compiled
  # path is the default and the one to trust.
  a <- spherical_roi(sp1, c(6L, 6L, 6L), 4, use_cpp = TRUE)
  expect_equal(nrow(a@coords), nrow(b@coords))
  expect_true(any(a@coords[, 1] == 6 & a@coords[, 2] == 6 & a@coords[, 3] == 6))
  ka <- apply(unname(a@coords), 1, paste, collapse = ",")
  kb <- apply(unname(b@coords), 1, paste, collapse = ",")
  expect_true(length(setdiff(ka, kb)) > 0)
})

test_that("the offset cache does not confuse distinct radii or spacings", {
  sp1 <- NeuroSpace(c(20L, 20L, 20L), c(2, 2, 2))
  r3 <- spherical_roi(sp1, c(10L, 10L, 10L), 3)
  r6 <- spherical_roi(sp1, c(10L, 10L, 10L), 6)
  r3b <- spherical_roi(sp1, c(10L, 10L, 10L), 3)

  expect_equal(r3@coords, r3b@coords)
  expect_true(nrow(r6@coords) > nrow(r3@coords))

  # near-identical radii must not collide in the cache
  a <- neuroim2:::sphere_offsets_cpp(4, c(2, 2, 2))
  b <- neuroim2:::sphere_offsets_cpp(4 + 1e-12, c(2, 2, 2))
  expect_equal(nrow(a), nrow(b))

  # same radius, different spacing
  s1 <- neuroim2:::sphere_offsets_cpp(6, c(1, 1, 1))
  s2 <- neuroim2:::sphere_offsets_cpp(6, c(1, 1, 4))
  expect_false(identical(s1, s2))
})

# ---------------------------------------------------------------- series()

test_that("series(DenseNeuroVec, matrix) matches the per-timepoint R loop", {
  ref_dense <- function(x, i) {
    d <- dim(x)
    i <- matrix(as.integer(i), ncol = 3)
    lin <- (i[, 3] - 1L) * d[1] * d[2] + (i[, 2] - 1L) * d[1] + i[, 1]
    nt <- d[4]
    nels <- prod(d[1:3])
    out <- matrix(0, nt, nrow(i))
    for (t in seq_len(nt)) out[t, ] <- x@.Data[lin + (t - 1L) * nels]
    out
  }

  set.seed(21)
  for (d4 in list(c(8L, 7L, 6L, 1L), c(8L, 7L, 6L, 11L), c(5L, 5L, 5L, 4L))) {
    arr <- array(rnorm(prod(d4)), d4)
    arr[1] <- NA_real_; arr[2] <- NaN; arr[3] <- Inf   # non-finite data must pass through
    dv <- NeuroVec(arr, NeuroSpace(d4, c(1, 1, 1)))

    cds <- rbind(c(1L, 1L, 1L), d4[1:3], c(2L, 3L, 4L),
                 cbind(sample.int(d4[1], 6, TRUE),
                       sample.int(d4[2], 6, TRUE),
                       sample.int(d4[3], 6, TRUE)))
    storage.mode(cds) <- "integer"

    expect_equal(series(dv, cds), ref_dense(dv, cds))
    expect_equal(series(dv, cds[1, , drop = FALSE]), ref_dense(dv, cds[1, , drop = FALSE]))
  }
})

test_that("series(SparseNeuroVec, ...) matches the zero-fill-and-scatter reference", {
  ref_sparse <- function(x, idx) {
    mapped <- lookup(x, idx)
    out <- matrix(0, nrow = dim(x)[4], ncol = length(idx))
    nz <- which(mapped > 0)
    if (length(nz) > 0) out[, nz] <- matricized_access(x, mapped[nz])
    out
  }

  set.seed(22)
  d <- c(9L, 8L, 7L); nt <- 6L
  m <- array(runif(prod(d)) > 0.4, d)
  mask <- LogicalNeuroVol(m, NeuroSpace(d, c(1, 1, 1)))
  dat <- matrix(rnorm(sum(m) * nt), nt, sum(m))
  dat[1, 1] <- NA_real_; dat[2, 1] <- NaN
  sv <- SparseNeuroVec(dat, NeuroSpace(c(d, nt), c(1, 1, 1)), mask)

  in_mask <- which(m)
  out_mask <- which(!m)

  cases <- list(
    all_in    = as.integer(head(in_mask, 12)),
    all_out   = as.integer(head(out_mask, 8)),
    mixed     = as.integer(c(head(in_mask, 5), head(out_mask, 5))),
    single_in = as.integer(in_mask[1]),
    single_out= as.integer(out_mask[1])
  )

  for (nm in names(cases)) {
    idx <- cases[[nm]]
    expect_equal(series(sv, idx, drop = FALSE), ref_sparse(sv, idx), info = nm)
  }

  # drop semantics for a single voxel
  expect_true(is.null(dim(series(sv, as.integer(in_mask[1])))))
  expect_equal(dim(series(sv, as.integer(in_mask[1]), drop = FALSE)), c(nt, 1L))

  # i/j/k form
  g <- index_to_grid(mask, as.integer(head(in_mask, 4)))
  expect_equal(series(sv, g), ref_sparse(sv, as.integer(head(in_mask, 4))))
})

test_that("lazy sparse subclasses still route through matricized_access", {
  # The compiled gather reads @data directly, so it must not engage for a
  # subclass whose real values live behind an overridden accessor.
  setClass("FastPathProbeVec",
           slots = c(data = "matrix", backing = "matrix"),
           contains = c("AbstractSparseNeuroVec"))
  on.exit(try(removeClass("FastPathProbeVec"), silent = TRUE), add = TRUE)

  setMethod("matricized_access", signature(x = "FastPathProbeVec", i = "integer"),
            function(x, i, ...) x@backing[, i, drop = FALSE])
  setMethod("matricized_access", signature(x = "FastPathProbeVec", i = "numeric"),
            function(x, i, ...) x@backing[, i, drop = FALSE])

  dims <- c(2L, 2L, 1L, 3L)
  mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = dims[1:3])
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(dims[1:3]))
  backing <- matrix(c(10, 20, 30, 40, 50, 60), nrow = dims[4])

  x <- new("FastPathProbeVec",
           space = NeuroSpace(dims), label = "probe", mask = mask,
           map = IndexLookupVol(space(mask), as.integer(which(mask_arr))),
           data = matrix(-999, nrow = dims[4], ncol = sum(mask_arr)),
           backing = backing)

  # -999 is the placeholder; seeing it means the fast path wrongly engaged.
  expect_equal(series(x, 1L), c(10, 20, 30))
  expect_false(any(series(x, 1L) == -999))
  expect_null(neuroim2:::.sparse_series_fast(x, 1L))
})

test_that("the dense gather rejects a dim that overflows, rather than reading OOB", {
  # dim() on a NeuroVec reports the NeuroSpace slot, which R does not keep in
  # sync with the length of the underlying array -- so dim can claim more
  # elements than the data holds. Computing prod(dim) to check that can
  # overflow R_xlen_t and wrap negative, letting a bogus dim past the guard and
  # segfaulting on the read. These are the exact overflow points.
  msg <- "shorter than prod\\(dim\\)"

  # d0*d1*d2 == 2^63 exactly
  D <- 2097152L
  dv <- DenseNeuroVec(array(as.numeric(1:8), c(2, 2, 2, 1)),
                      NeuroSpace(c(D, D, D, 1L)))
  expect_error(series(dv, cbind(1L, 1L, 2L)), msg)

  # overflow in the nt stride rather than the spatial one
  dv2 <- DenseNeuroVec(array(as.numeric(1:8), c(2, 2, 2, 1)),
                       NeuroSpace(c(208064L, 208064L, 208064L, 1024L)))
  expect_error(series(dv2, cbind(1L, 1L, 1L)), msg)

  # reached by a single slot assignment on an otherwise valid object
  dv3 <- DenseNeuroVec(array(rnorm(128), c(4, 4, 4, 2)), NeuroSpace(c(4, 4, 4, 2)))
  dv3@space <- NeuroSpace(c(D, D, D, 1L))
  expect_error(series(dv3, cbind(1L, 1L, 2L)), msg)

  # one below the overflow point must still be caught by the plain guard
  D2 <- 2097151L
  dv4 <- DenseNeuroVec(array(as.numeric(1:8), c(2, 2, 2, 1)),
                       NeuroSpace(c(D2, D2, D2, 1L)))
  expect_error(series(dv4, cbind(1L, 1L, 2L)), msg)

  # and a well-formed object is unaffected
  ok <- DenseNeuroVec(array(rnorm(128), c(4, 4, 4, 2)), NeuroSpace(c(4, 4, 4, 2)))
  expect_equal(dim(series(ok, cbind(1L, 1L, 1L))), c(2L, 1L))
})

test_that("the sparse fast path is gated on the exact class, package included", {
  d <- c(2L, 2L, 1L)
  m <- array(c(TRUE, FALSE, TRUE, FALSE), dim = d)
  mask <- LogicalNeuroVol(m, NeuroSpace(d))
  sv <- SparseNeuroVec(matrix(c(1, 2, 3, 4, 5, 6), 3, 2), NeuroSpace(c(d, 3L)), mask)

  # a genuine SparseNeuroVec takes the fast path
  expect_false(is.null(neuroim2:::.sparse_series_fast(sv, 1L)))

  # an object whose class attribute merely *says* SparseNeuroVec (different
  # package) must not: class(x)[1L] would drop the package qualifier.
  bad <- sv
  attr(bad, "class") <- structure("SparseNeuroVec", package = ".GlobalEnv")
  expect_null(neuroim2:::.sparse_series_fast(bad, 1L))

  # a store whose row count disagrees with dim(x)[4] must not either: the
  # implementation this replaced returned a dim(x)[4]-row result.
  odd <- sv
  odd@data <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_null(neuroim2:::.sparse_series_fast(odd, 1L))
})

test_that("the dense gather does not scale with volume size", {
  # Guards against re-introducing a whole-volume duplicate on every call:
  # passing x@.Data (rather than x) into .Call copies the entire image.
  skip_on_cran()

  roi <- matrix(as.integer(c(5, 5, 5)), 1, 3)
  timing <- function(d4) {
    dv <- NeuroVec(array(0, d4), NeuroSpace(d4, c(1, 1, 1)))
    series(dv, roi)
    median(replicate(3, system.time(for (k in 1:20) series(dv, roi))[["elapsed"]])) / 20
  }

  small <- timing(c(20L, 20L, 20L, 30L))
  large <- timing(c(80L, 80L, 80L, 30L))   # 64x more voxels

  # Allow generous slack for noise, but a full copy would be ~64x slower.
  expect_lt(large, max(small * 10, 0.01))
})

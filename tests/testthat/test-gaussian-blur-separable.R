library(testthat)
library(neuroim2)

# gaussian_blur() dispatches between two kernels that must agree bit-for-bit up
# to floating-point accumulation order: the dense (2w+1)^3 kernel evaluated at
# mask voxels, and three 1-D passes over the whole volume. These tests pin the
# equivalence, the dispatch rule, and the separable kernel's argument checking.

dense_blur <- function(...) neuroim2:::gaussian_blur_cpp(...)
sep_blur   <- function(...) neuroim2:::gaussian_blur_sep_cpp(...)
prefers    <- neuroim2:::.gaussian_blur_prefers_separable

test_that("separable and dense kernels agree across windows, spacings and sigmas", {
  set.seed(41)
  D <- c(11L, 9L, 13L)
  arr <- array(rnorm(prod(D)), D)
  idx <- sort(sample.int(prod(D), prod(D) %/% 2L))

  for (sp in list(c(1, 1, 1), c(2, 2, 2), c(3, 1.5, 0.75), c(0.5, 0.5, 4))) {
    for (sg in c(0.5, 1, 2, 5)) {
      for (w in 1:4) {
        for (nz in c(TRUE, FALSE)) {
          a <- dense_blur(arr, idx, w, sg, sp, nz)
          b <- sep_blur(arr, idx, w, sg, sp, nz)
          lbl <- sprintf("sp=%s sigma=%g window=%d normalize=%s",
                         paste(sp, collapse = ","), sg, w, nz)
          expect_identical(dim(b), dim(a), info = lbl)
          # Absolute agreement: the values here are O(1), and the only
          # difference permitted is summation order.
          expect_lt(max(abs(a - b)), 1e-12, label = lbl)
        }
      }
    }
  }
})

test_that("separable kernel handles degenerate and full-volume masks", {
  set.seed(42)
  for (D in list(c(1L, 10L, 10L), c(4L, 4L, 1L), c(1L, 1L, 1L), c(7L, 7L, 7L))) {
    arr <- array(rnorm(prod(D)), D)
    for (idx in list(seq_len(prod(D)), prod(D) %/% 2L + 1L)) {
      for (nz in c(TRUE, FALSE)) {
        a <- dense_blur(arr, idx, 2L, 1.5, c(1, 1, 1), nz)
        b <- sep_blur(arr, idx, 2L, 1.5, c(1, 1, 1), nz)
        expect_lt(max(abs(a - b)), 1e-12,
                  label = sprintf("dim=%s n_mask=%d normalize=%s",
                                  paste(D, collapse = "x"), length(idx), nz))
      }
    }
  }
})

test_that("separable kernel does not leak out-of-mask NaN into the mask", {
  set.seed(44)
  D <- c(14L, 14L, 14L)
  mk <- array(FALSE, D); mk[5:10, , ] <- TRUE
  idx <- which(mk)

  # A *structured* in-mask signal. A constant would be reproduced by any
  # normalized kernel regardless of its shape, padding or axis order, so it
  # cannot distinguish a correct kernel from a broken one.
  arr <- array(NaN, D)                            # NaN everywhere outside
  arr[idx] <- rnorm(length(idx))

  for (nonfinite in list(NaN, NA_real_, Inf, -Inf)) {
    a <- arr; a[-idx] <- nonfinite
    out <- sep_blur(a, idx, 3L, 1.27, c(2, 2, 2), TRUE)
    lbl <- sprintf("outside = %s", format(nonfinite))
    # Nothing outside the mask may reach an in-mask output, whatever it holds.
    expect_false(any(is.nan(out[idx])), label = lbl)
    expect_true(all(is.finite(out[idx])), label = lbl)
    expect_true(all(out[-idx] == 0), label = lbl)
    # The answer must not depend on what is outside the mask at all.
    expect_lt(max(abs(out[idx] - sep_blur(arr, idx, 3L, 1.27, c(2, 2, 2), TRUE)[idx])),
              1e-12, label = lbl)
    # ... and must match the dense kernel, where the guarantee is defined.
    expect_lt(max(abs(out - dense_blur(a, idx, 3L, 1.27, c(2, 2, 2), TRUE))),
              1e-12, label = lbl)
  }
})

test_that("a constant in-mask signal survives renormalization at mask edges", {
  D <- c(14L, 14L, 14L)
  mk <- array(FALSE, D); mk[5:10, , ] <- TRUE
  arr <- array(NaN, D); arr[mk] <- 1
  idx <- which(mk)
  # Renormalizing by the in-mask weight means edge voxels are not biased toward
  # the exterior, so a constant comes back unchanged right up to the boundary.
  out <- sep_blur(arr, idx, 3L, 1.27, c(2, 2, 2), TRUE)
  expect_equal(out[idx], rep(1, length(idx)), tolerance = 1e-12)
})

test_that("normalize=TRUE preserves a single-voxel mask exactly", {
  D <- c(9L, 9L, 9L)
  v <- 365L                                        # interior voxel (5,5,5)
  # Every neighbour is out of mask and hostile: the renormalized kernel must
  # reduce to the centre weight alone, so the value passes through untouched.
  arr <- array(c(NaN, NA_real_, Inf, -Inf), D)
  arr[v] <- 0.5
  out <- sep_blur(arr, v, 2L, 2, c(1, 1, 1), TRUE)
  expect_equal(out[v], 0.5, tolerance = 1e-14)
  expect_equal(sum(out != 0), 1L)
  expect_equal(out[v], dense_blur(arr, v, 2L, 2, c(1, 1, 1), TRUE)[v], tolerance = 1e-14)
})

test_that("gaussian_blur() end-to-end matches the dense kernel on both branches", {
  set.seed(43)
  D <- c(18L, 20L, 16L)
  sp <- NeuroSpace(D, spacing = c(2, 2, 2))
  arr <- array(rnorm(prod(D)), D)
  vol <- NeuroVol(arr, sp)

  # A ~30% mask at window 1 takes the separable branch; a ~2% mask at window 1
  # stays on the dense one. Both must produce the dense kernel's answer.
  big <- array(FALSE, D); big[which(runif(prod(D)) < 0.3)] <- TRUE
  small <- array(FALSE, D); small[sample.int(prod(D), round(0.02 * prod(D)))] <- TRUE

  for (mk in list(big, small)) {
    idx <- which(mk)
    for (w in 1:3) {
      got <- gaussian_blur(vol, LogicalNeuroVol(mk, sp), sigma = 2, window = w)
      ref <- dense_blur(arr, idx, w, 2, c(2, 2, 2), TRUE)
      expect_lt(max(abs(as.numeric(got) - as.numeric(ref))), 1e-12,
                label = sprintf("n_mask=%d window=%d", length(idx), w))
    }
  }

  # No mask at all: the whole volume is the mask, so this is the separable
  # branch for every window.
  for (w in 1:3) {
    got <- gaussian_blur(vol, sigma = 2, window = w)
    ref <- dense_blur(arr, seq_len(prod(D)), w, 2, c(2, 2, 2), TRUE)
    expect_lt(max(abs(as.numeric(got) - as.numeric(ref))), 1e-12,
              label = sprintf("unmasked window=%d", w))
  }
})

test_that("the dispatch rule matches its calibration", {
  n <- 1e6

  # The thresholds themselves, which is the part of this change that was fitted
  # to measurement rather than derived. If the constant or the exponent in
  # .gaussian_blur_prefers_separable() moves, this fails and the calibration
  # needs re-measuring -- see the comment block in R/spat_filter.R.
  thresh <- function(w, passes) 0.15 * passes / (2 * w + 1)^1.5
  expect_equal(vapply(1:4, thresh, numeric(1), passes = 6),
               c(0.17320508, 0.08049845, 0.04859611, 0.03333333), tolerance = 1e-6)
  expect_equal(vapply(1:4, thresh, numeric(1), passes = 3),
               c(0.08660254, 0.04024922, 0.02429805, 0.01666667), tolerance = 1e-6)
  for (w in 1:4) for (nz in c(TRUE, FALSE)) {
    t0 <- thresh(w, if (nz) 6 else 3)
    expect_false(prefers(w, (t0 - 1e-6) * n, n, nz),
                 label = sprintf("just below threshold w=%d normalize=%s", w, nz))
    expect_true(prefers(w, (t0 + 1e-6) * n, n, nz),
                label = sprintf("just above threshold w=%d normalize=%s", w, nz))
  }

  # Degenerate arguments choose the dense path rather than erroring: this
  # function only picks an engine, and validation belongs to the kernels.
  for (bad in list(0L, -1L, NA_integer_, NA_real_, NULL, numeric(0), c(1, 2), Inf, NaN)) {
    expect_false(prefers(bad, n, n, TRUE), label = sprintf("window=%s", format(bad)[1]))
  }
  expect_false(prefers(1L, 0, 0, TRUE))       # empty volume
  expect_false(prefers(1L, 0, n, TRUE))       # empty mask
  expect_false(prefers(1L, n, NA, TRUE))

  # A whole-volume mask always prefers separable; a 1-voxel mask never does.
  for (w in 1:4) for (nz in c(TRUE, FALSE)) {
    expect_true(prefers(w, n, n, nz))
    expect_false(prefers(w, 1, n, nz))
  }

  # Monotone in the mask fraction and in the window, and normalize=TRUE (twice
  # the passes) never prefers separable where normalize=FALSE does not.
  for (w in 1:4) for (nz in c(TRUE, FALSE)) {
    p <- vapply(seq(0, 1, by = 0.05), function(f) prefers(w, f * n, n, nz), logical(1))
    expect_false(any(diff(p) < 0),
                 label = sprintf("monotone in fraction w=%d normalize=%s", w, nz))
  }
  for (f in c(0.05, 0.2, 0.5)) {
    p <- vapply(1:8, function(w) prefers(w, f * n, n, TRUE), logical(1))
    expect_false(any(diff(p) < 0), label = sprintf("monotone in window frac=%g", f))
    expect_true(all(prefers(2L, f * n, n, FALSE) >= prefers(2L, f * n, n, TRUE)))
  }

  # A brain-sized mask at the default window: the case the calibration exists
  # for, and the one the tap-count rule got wrong.
  expect_true(prefers(1L, 0.30 * n, n, TRUE))
  expect_true(prefers(1L, 0.30 * n, n, FALSE))
})

test_that("the engine agrees with the dense kernel on awkward arguments", {
  set.seed(45)
  D <- c(8L, 7L, 9L)
  nvox <- prod(D)
  arr <- array(rnorm(nvox), D)
  eng <- neuroim2:::.gaussian_blur_engine

  # Windows past the volume, where the kernel is entirely clipped on some axes.
  # Nothing else in this file exercises window > 4.
  for (w in c(5L, 9L, 20L)) for (nz in c(TRUE, FALSE)) {
    expect_lt(max(abs(eng(arr, seq_len(nvox), w, 2, c(1, 1, 1), nz) -
                      dense_blur(arr, seq_len(nvox), w, 2, c(1, 1, 1), nz))), 1e-12,
              label = sprintf("window=%d normalize=%s", w, nz))
  }

  # Duplicated, unsorted, and double-typed indices, and an integer `arr`.
  idx <- sort(sample.int(nvox, nvox %/% 2L))
  for (mi in list(rev(idx), c(idx, idx), as.numeric(idx))) {
    for (nz in c(TRUE, FALSE)) {
      expect_lt(max(abs(eng(arr, mi, 2L, 2, c(1, 1, 1), nz) -
                        dense_blur(arr, as.integer(mi), 2L, 2, c(1, 1, 1), nz))), 1e-12,
                label = sprintf("n_idx=%d normalize=%s", length(mi), nz))
    }
  }
  iarr <- array(as.integer(round(arr * 10)), D)
  for (nz in c(TRUE, FALSE)) {
    expect_lt(max(abs(eng(iarr, seq_len(nvox), 2L, 2, c(1, 1, 1), nz) -
                      dense_blur(iarr, seq_len(nvox), 2L, 2, c(1, 1, 1), nz))), 1e-12,
              label = sprintf("integer arr normalize=%s", nz))
  }

  # An empty mask: dense path, all zeros, no error.
  expect_equal(sum(eng(arr, integer(0), 2L, 2, c(1, 1, 1), TRUE) != 0), 0L)

  # normalize=FALSE deliberately reads out-of-mask neighbours, so non-finite
  # values there DO propagate -- and must propagate identically on both paths.
  for (nonfinite in list(NaN, NA_real_, Inf, -Inf)) {
    a <- arr; a[-idx] <- nonfinite
    for (mi in list(idx, seq_len(nvox))) {
      got <- eng(a, mi, 2L, 2, c(1, 1, 1), FALSE)
      ref <- dense_blur(a, mi, 2L, 2, c(1, 1, 1), FALSE)
      lbl <- sprintf("legacy branch, outside = %s, n_idx=%d", format(nonfinite), length(mi))
      expect_identical(is.nan(got), is.nan(ref), info = lbl)
      expect_identical(is.finite(got), is.finite(ref), info = lbl)
      ok <- is.finite(got)
      if (any(ok)) expect_lt(max(abs(got[ok] - ref[ok])), 1e-12, label = lbl)
    }
  }

  # A 4-D array must not start erroring: the separable kernel requires exactly
  # three dimensions, so the engine has to keep such input on the dense path.
  a4 <- array(rnorm(2 * nvox), c(D, 2L))
  expect_identical(eng(a4, seq_len(2 * nvox), 1L, 2, c(1, 1, 1), TRUE),
                   dense_blur(a4, seq_len(2 * nvox), 1L, 2, c(1, 1, 1), TRUE))
})

test_that("gaussian_blur() validates its arguments and clamps absurd windows", {
  D <- c(10L, 12L, 8L)
  sp <- NeuroSpace(D, spacing = c(1, 1, 1))
  vol <- NeuroVol(array(rnorm(prod(D)), D), sp)
  mk <- LogicalNeuroVol(array(TRUE, D), sp)

  expect_error(gaussian_blur(vol, mk, sigma = 2, window = Inf), "finite")
  expect_error(gaussian_blur(vol, mk, sigma = 2, window = c(1, 2)), "single")
  # window = NULL is no longer an error: it is the default, and means "derive
  # the half-width from sigma and the voxel size".
  expect_silent(gaussian_blur(vol, mk, sigma = 2, window = NULL))
  expect_equal(as.numeric(gaussian_blur(vol, mk, sigma = 2, window = NULL)),
               as.numeric(gaussian_blur(vol, mk, sigma = 2, window = 8)))
  expect_error(gaussian_blur(vol, mk, sigma = 2, window = NA), "single")
  expect_error(gaussian_blur(vol, mk, sigma = Inf), "finite")
  expect_error(gaussian_blur(vol, mk, sigma = NA), "finite")
  expect_error(gaussian_blur(vol, mk, sigma = c(1, 2)), "single")

  # A mask of the wrong shape used to reach an unbounded write and segfault.
  big <- LogicalNeuroVol(array(TRUE, c(20L, 20L, 20L)),
                         NeuroSpace(c(20L, 20L, 20L), spacing = c(1, 1, 1)))
  expect_error(gaussian_blur(vol, big), "same spatial dimensions")

  # Any tap at least max(dim) away is out of bounds from every voxel, so
  # clamping the window there changes no result.
  expect_equal(as.numeric(gaussian_blur(vol, mk, sigma = 2, window = 400)),
               as.numeric(gaussian_blur(vol, mk, sigma = 2, window = 12)),
               tolerance = 1e-12)
  expect_equal(as.numeric(gaussian_blur(vol, mk, sigma = 2, window = 1e10)),
               as.numeric(gaussian_blur(vol, mk, sigma = 2, window = 12)),
               tolerance = 1e-12)
})

test_that("gaussian_weights() no longer overflows at extreme sigma", {
  # The 3-D kernel used to carry a per-axis sqrt(2*pi*sigma) that cancelled in
  # the normalization but cubed to Inf first, making every weight NaN.
  for (sg in c(1e-8, 1, 1e100, 1e203, 1e204, 1e250, 1e300)) {
    w <- neuroim2:::gaussian_weights(1L, sg, c(1, 1, 1))
    expect_false(any(is.nan(w)), label = sprintf("sigma=%g", sg))
    expect_equal(sum(w), 1, tolerance = 1e-12, label = sprintf("sigma=%g", sg))
  }
  # ... and the two kernels agree there, which they did not before.
  set.seed(46)
  D <- c(9L, 8L, 7L); arr <- array(rnorm(prod(D)), D); idx <- seq_len(prod(D))
  for (sg in c(1e203, 1e204, 1e250)) {
    expect_lt(max(abs(dense_blur(arr, idx, 2L, sg, c(1, 1, 1), TRUE) -
                      sep_blur(arr, idx, 2L, sg, c(1, 1, 1), TRUE))), 1e-12,
              label = sprintf("sigma=%g", sg))
  }
  expect_error(neuroim2:::gaussian_weights(-1L, 2, c(1, 1, 1)), "non-negative")
})

test_that("dense kernel rejects an out-of-bounds mask index", {
  # This was an unbounded std::vector write, i.e. a segfault, reachable from
  # gaussian_blur() whenever mask and volume disagreed in size.
  arr <- array(0, c(6L, 6L, 6L))
  expect_error(dense_blur(arr, 217L, 1L, 2, c(1, 1, 1), TRUE), "out of bounds")
  expect_error(dense_blur(arr, 0L, 1L, 2, c(1, 1, 1), TRUE), "out of bounds")
  expect_error(dense_blur(arr, -5L, 1L, 2, c(1, 1, 1), FALSE), "out of bounds")
  expect_error(dense_blur(arr, 1L, 1L, 2, c(1, 1), TRUE), "spacing")
  expect_error(dense_blur(1:10, 1L, 1L, 2, c(1, 1, 1), TRUE), "dim")
})

test_that("separable kernel rejects malformed arguments", {
  D <- c(6L, 6L, 6L)
  arr <- array(0, D)
  expect_error(sep_blur(arr, 1L, 1L, 2, c(1, 1)), "spacing")
  expect_error(sep_blur(arr, 1L, 1L, 0, c(1, 1, 1)), "sigma")
  expect_error(sep_blur(arr, 1L, 1L, NaN, c(1, 1, 1)), "sigma")
  expect_error(sep_blur(arr, 1L, -1L, 2, c(1, 1, 1)), "window")
  expect_error(sep_blur(arr, 0L, 1L, 2, c(1, 1, 1)), "out of bounds")
  expect_error(sep_blur(arr, 217L, 1L, 2, c(1, 1, 1)), "out of bounds")
  expect_error(sep_blur(arr, -5L, 1L, 2, c(1, 1, 1)), "out of bounds")
  # ... on the unnormalized branch too, which checks bounds in its own loop.
  expect_error(sep_blur(arr, 217L, 1L, 2, c(1, 1, 1), FALSE), "out of bounds")

  expect_error(sep_blur(1:10, 1L, 1L, 2, c(1, 1, 1)), "dim")
  bad <- array(0, c(6L, 6L)); expect_error(sep_blur(bad, 1L, 1L, 2, c(1, 1, 1)), "dim")
})

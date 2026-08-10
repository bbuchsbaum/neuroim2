context("processing defects found against nilearn")

library(neuroim2)

ns <- asNamespace("neuroim2")

# The point-spread function is the operator: smooth an impulse and read the
# FWHM back out. Comparing the `sigma` that went in tells you nothing about the
# smoothing that came out, which is exactly how the truncation went unnoticed.
psf_fwhm <- function(vol, spacing, centre) {
  a <- as.array(vol)
  a[a < 0] <- 0
  prof <- apply(a, 1, sum)
  off <- (seq_along(prof) - centre) * spacing
  m <- sum(prof)
  mu <- sum(prof * off) / m
  2 * sqrt(2 * log(2)) * sqrt(sum(prof * (off - mu)^2) / m)
}

impulse <- function(spacing, n = 41L) {
  sp <- NeuroSpace(rep(n, 3L), rep(spacing, 3))
  d <- array(0, rep(n, 3L))
  d[(n + 1L) %/% 2L, (n + 1L) %/% 2L, (n + 1L) %/% 2L] <- 1
  NeuroVol(d, sp)
}

# ---------------------------------------------------------------------------
# gaussian_blur(): the kernel is sized from sigma
# ---------------------------------------------------------------------------

test_that("the default kernel delivers the smoothing that was requested", {
  for (vox in c(1, 2, 3)) {
    v <- impulse(vox)
    for (sg in c(2, 3, 4)) {
      got <- psf_fwhm(gaussian_blur(v, sigma = sg, normalize = FALSE), vox, 21L)
      want <- 2 * sqrt(2 * log(2)) * sg
      # 4 sigma of a Gaussian is 99.99% of its mass; the residual truncation is
      # well under a percent of the FWHM.
      expect_equal(got, want, tolerance = 0.01,
                   info = paste("spacing", vox, "sigma", sg))
    }
  }
})

test_that("fwhm is accepted and means what it says", {
  v <- impulse(2)
  for (f in c(4, 6, 9)) {
    expect_equal(psf_fwhm(gaussian_blur(v, fwhm = f, normalize = FALSE), 2, 21L),
                 f, tolerance = 0.01, info = paste("fwhm", f))
  }
  # fwhm and sigma are two spellings of one argument
  expect_equal(as.numeric(gaussian_blur(v, fwhm = 6, normalize = FALSE)),
               as.numeric(gaussian_blur(v, sigma = 6 / (2 * sqrt(2 * log(2))),
                                        normalize = FALSE)))
  expect_error(gaussian_blur(v, sigma = 2, fwhm = 6), "not both")
  expect_error(gaussian_blur(v, fwhm = -1), "positive")
  expect_error(gaussian_blur(v, fwhm = c(1, 2)), "single")
})

test_that("an explicit window still truncates exactly as it used to", {
  v <- impulse(2)
  # These are the numbers the truncated kernel produced before the default
  # changed; passing `window` must keep reproducing them.
  expect_equal(psf_fwhm(gaussian_blur(v, sigma = 2, window = 1, normalize = FALSE), 2, 21L),
               3.49, tolerance = 0.01)
  expect_equal(psf_fwhm(gaussian_blur(v, sigma = 2, window = 2, normalize = FALSE), 2, 21L),
               4.53, tolerance = 0.01)
  expect_equal(psf_fwhm(gaussian_blur(v, sigma = 4, window = 1, normalize = FALSE), 2, 21L),
               3.76, tolerance = 0.01)
})

test_that("the derived window follows scipy's truncate rule", {
  expect_equal(ns$.gaussian_window(2, c(2, 2, 2)), 4L)
  expect_equal(ns$.gaussian_window(2, c(1, 1, 1)), 8L)
  # Anisotropic voxels take the finest axis, so no axis is clipped early.
  expect_equal(ns$.gaussian_window(3, c(1, 1, 4)), 12L)
  expect_equal(ns$.gaussian_window(2, c(2, 2, 2), truncate = 2), 2L)
  expect_gte(ns$.gaussian_window(0.01, c(2, 2, 2)), 1L)
  # A degenerate spacing must not produce a degenerate window.
  expect_gte(ns$.gaussian_window(2, c(0, 0, 0)), 1L)
  expect_error(gaussian_blur(impulse(2), sigma = 2, truncate = 0), "positive")
})

test_that("smoothing a 4-D image still matches smoothing its volumes", {
  sp <- NeuroSpace(c(9L, 10L, 8L, 3L), c(2, 2, 3))
  set.seed(8)
  x <- DenseNeuroVec(array(rnorm(prod(dim(sp))), dim(sp)), sp)
  got <- gaussian_blur(x, sigma = 2)
  want <- vapply(seq_len(3L),
                 function(i) as.numeric(gaussian_blur(x[[i]], sigma = 2)),
                 numeric(prod(dim(sp)[1:3])))
  expect_equal(as.numeric(got@.Data), as.numeric(want))
})

# ---------------------------------------------------------------------------
# as_canonical(): a permutation, not a resample
# ---------------------------------------------------------------------------

odd_space <- function(dm) {
  NeuroSpace(dm, c(3, 4, 2),
             trans = rbind(c(0, 0, 2, -10), c(3, 0, 0, -20),
                           c(0, -4, 0, 30), c(0, 0, 0, 1)))
}

test_that("reorientation permutes the grid instead of resampling into it", {
  dm <- c(3L, 4L, 5L)
  v <- NeuroVol(array(as.numeric(seq_len(prod(dm))), dm), odd_space(dm))
  expect_equal(axcodes(space(v)), c("A", "I", "R"))

  cn <- as_canonical(v)
  expect_equal(axcodes(space(cn)), c("R", "A", "S"))
  # The axes permute, so the shape does too. Keeping the source shape is what
  # made the old version write into a grid that did not hold the data.
  expect_equal(dim(cn), c(5L, 3L, 4L))
  # Exact: every voxel survives, and no value is interpolated.
  expect_identical(sort(as.numeric(as.array(v))), sort(as.numeric(as.array(cn))))
  # Volumes below RNiftyReg's four-voxel minimum used to error outright.
  expect_s4_class(cn, "DenseNeuroVol")
})

test_that("the reoriented affine still maps voxels to the same world points", {
  dm <- c(3L, 4L, 5L)
  v <- NeuroVol(array(as.numeric(seq_len(prod(dm))), dm), odd_space(dm))
  cn <- as_canonical(v)
  # Pick a voxel by value and check it sits at the same world coordinate.
  for (val in c(1, 17, 60)) {
    a <- which(as.array(v) == val, arr.ind = TRUE)
    b <- which(as.array(cn) == val, arr.ind = TRUE)
    expect_equal(as.numeric(grid_to_coord(v, a)),
                 as.numeric(grid_to_coord(cn, b)),
                 tolerance = 1e-6, info = paste("value", val))
  }
})

test_that("reorientation round-trips and is a no-op when already canonical", {
  dm <- c(3L, 4L, 5L)
  v <- NeuroVol(array(as.numeric(seq_len(prod(dm))), dm), odd_space(dm))
  cn <- as_canonical(v)
  back <- as_canonical(cn, axcodes(space(v)))
  expect_identical(as.numeric(as.array(back)), as.numeric(as.array(v)))
  expect_equal(dim(back), dim(v))
  expect_identical(as_canonical(cn), cn)
})

test_that("a 4-D image reorients volume by volume", {
  dm <- c(3L, 4L, 5L, 3L)
  sp <- add_dim(odd_space(dm[1:3]), 3L)
  x <- DenseNeuroVec(array(as.numeric(seq_len(prod(dm))), dm), sp)
  cx <- as_canonical(x)
  expect_equal(dim(cx), c(5L, 3L, 4L, 3L))
  expect_identical(sort(as.numeric(as.array(x))), sort(as.numeric(as.array(cx))))
  for (i in 1:3) {
    expect_identical(as.numeric(as.array(cx[[i]])),
                     as.numeric(as.array(as_canonical(x[[i]]))), info = i)
  }
})

test_that("reorient() permutes a NeuroSpace's dimensions and spacing", {
  sp <- odd_space(c(3L, 4L, 5L))
  r <- reorient(sp, c("R", "A", "S"))
  expect_equal(dim(r), c(5L, 3L, 4L))
  expect_equal(spacing(r), spacing(sp)[c(3, 1, 2)])
  expect_equal(affine_to_axcodes(trans(r)), c("R", "A", "S"))
  # Asking for the orientation a space already has changes nothing.
  same <- reorient(r, c("R", "A", "S"))
  expect_equal(dim(same), dim(r))
  expect_equal(trans(same), trans(r))
})

# ---------------------------------------------------------------------------
# searchlight iterators
# ---------------------------------------------------------------------------

test_that("searchlight() and searchlight_coords() use the same centres", {
  sp <- NeuroSpace(c(10L, 11L, 9L), c(2, 2, 2))
  m <- array(FALSE, c(10, 11, 9))
  m[3:8, 3:9, 2:7] <- TRUE
  mask <- LogicalNeuroVol(m, sp)
  n_in <- sum(m)
  expect_lt(n_in, prod(dim(sp)))          # the two rules must be distinguishable

  for (nz in c(FALSE, TRUE)) {
    a <- searchlight(mask, radius = 5, nonzero = nz)
    b <- searchlight_coords(mask, radius = 5, nonzero = nz)
    expect_equal(length(a), n_in, info = paste("nonzero", nz))
    expect_equal(length(b), n_in, info = paste("nonzero", nz))
    for (i in c(1L, n_in %/% 2L, n_in)) {
      expect_identical(unname(coords(a[[i]])), unname(b[[i]]),
                       info = paste("nonzero", nz, "element", i))
    }
  }
})

test_that("every searchlight iterator agrees on how many centres there are", {
  sp <- NeuroSpace(c(9L, 9L, 9L), c(2, 2, 2))
  m <- array(FALSE, c(9, 9, 9)); m[2:7, 2:7, 2:7] <- TRUE
  mask <- LogicalNeuroVol(m, sp)
  n_in <- sum(m)
  expect_equal(length(searchlight(mask, radius = 4)), n_in)
  expect_equal(length(searchlight_coords(mask, radius = 4)), n_in)
  expect_equal(length(random_searchlight(mask, radius = 4)) > 0, TRUE)
})

test_that("searchlight_coords() validates its arguments", {
  sp <- NeuroSpace(c(6L, 6L, 6L), c(2, 2, 2))
  mask <- LogicalNeuroVol(array(TRUE, c(6, 6, 6)), sp)
  expect_error(searchlight_coords(as.array(mask), radius = 4), "NeuroVol")
  expect_error(searchlight_coords(mask, radius = 0), "positive")
  expect_error(searchlight_coords(mask, radius = 4, nonzero = 1), "TRUE or FALSE")
  expect_error(searchlight_coords(mask, radius = 4, cores = -1), "non-negative")
})

# ---------------------------------------------------------------------------
# resample target overlap
# ---------------------------------------------------------------------------

test_that("a target grid that misses the source is flagged", {
  sp <- NeuroSpace(c(12L, 12L, 12L), c(2, 2, 2),
                   trans = rbind(c(-2, 0, 0, 12), c(0, 2, 0, -12),
                                 c(0, 0, 2, -12), c(0, 0, 0, 1)))
  v <- NeuroVol(array(1, c(12, 12, 12)), sp)

  # Built from dims/spacing/origin: positive diagonal, so it mirrors the source.
  naive <- NeuroSpace(c(8L, 8L, 8L), c(3, 3, 3), origin = origin(v))
  expect_warning(resample(v, naive), "covers little of the source")

  # A target that keeps the source's axis directions does not warn.
  tx <- trans(v); tx[1:3, 1:3] <- tx[1:3, 1:3] %*% diag(rep(1.5, 3))
  good <- NeuroSpace(c(8L, 8L, 8L), c(3, 3, 3), trans = tx, origin = tx[1:3, 4])
  expect_silent(resample(v, good))
  expect_silent(resample(v, space(v)))
})

# ---------------------------------------------------------------------------
# reshapes that used to copy the payload
# ---------------------------------------------------------------------------

test_that("as.matrix() and mean() do not disturb the source image", {
  sp <- NeuroSpace(c(6L, 7L, 5L, 4L), c(2, 2, 2))
  vals <- as.numeric(seq_len(prod(dim(sp))))
  x <- DenseNeuroVec(array(vals, dim(sp)), sp)

  m <- as.matrix(x)
  expect_equal(dim(m), c(prod(dim(sp)[1:3]), 4L))
  expect_identical(as.numeric(m), vals)
  m[1, 1] <- -999                       # writing to the view must copy first
  expect_identical(as.numeric(as.array(x)), vals)
  expect_equal(dim(x), dim(sp))

  mu <- mean(x)
  expect_s4_class(mu, "DenseNeuroVol")
  expect_equal(dim(mu), dim(sp)[1:3])
  expect_equal(as.numeric(as.array(mu)),
               rowMeans(matrix(vals, nrow = prod(dim(sp)[1:3]), ncol = 4L)))
  expect_identical(as.numeric(as.array(x)), vals)

  # NA handling is unchanged
  vals2 <- vals; vals2[1] <- NA
  y <- DenseNeuroVec(array(vals2, dim(sp)), sp)
  expect_true(is.na(as.array(mean(y))[1]))
  expect_false(is.na(as.array(mean(y, na.rm = TRUE))[1]))
})

test_that("scale_series still centres and scales without touching the source", {
  sp <- NeuroSpace(c(5L, 5L, 4L, 6L), c(2, 2, 2))
  set.seed(12)
  vals <- rnorm(prod(dim(sp)))
  x <- DenseNeuroVec(array(vals, dim(sp)), sp)
  s <- scale_series(x, TRUE, TRUE)
  M <- matrix(vals, nrow = prod(dim(sp)[1:3]), ncol = 6L)
  M <- M - rowMeans(M)
  rsd <- sqrt(rowSums(M * M) / 5)
  rsd[rsd == 0] <- 1
  expect_equal(as.numeric(s@.Data), as.numeric(M / rsd))
  expect_identical(as.numeric(as.array(x)), vals)
})

test_that("the compiled 4-D driver is bit-identical to smoothing volume by volume", {
  # Both engines: a nearly-full mask takes the separable path, a sparse one the
  # dense kernel, which has no 4-D form and still loops.
  set.seed(19)
  for (dm in list(c(9L, 11L, 8L), c(13L, 7L, 10L))) {
    sp3 <- NeuroSpace(dm, c(2, 2, 2.5))
    for (nt in c(1L, 4L)) {
      a <- array(rnorm(prod(dm) * nt), c(dm, nt))
      v4 <- DenseNeuroVec(a, add_dim(sp3, nt))
      for (frac in c(0.6, 0.03)) {
        mk <- LogicalNeuroVol(array(runif(prod(dm)) < frac, dm), sp3)
        for (norm in c(TRUE, FALSE)) {
          for (w in c(1, 3)) {
            for (masked in c(TRUE, FALSE)) {
              got <- as.array(if (masked) gaussian_blur(v4, mk, sigma = 2, window = w, normalize = norm)
                              else gaussian_blur(v4, sigma = 2, window = w, normalize = norm))
              want <- vapply(seq_len(nt), function(i) {
                vi <- NeuroVol(array(a[, , , i], dm), sp3)
                as.array(if (masked) gaussian_blur(vi, mk, sigma = 2, window = w, normalize = norm)
                         else gaussian_blur(vi, sigma = 2, window = w, normalize = norm))
              }, array(0, dm))
              dim(want) <- c(dm, nt)
              expect_identical(got, want,
                               info = paste(paste(dm, collapse = "x"), nt, frac, norm, w, masked))
            }
          }
        }
      }
    }
  }
})

test_that("the 4-D driver rejects what it cannot smooth", {
  ns <- asNamespace("neuroim2")
  a <- array(1, c(3L, 3L, 3L))
  expect_error(ns$gaussian_blur_sep_4d_cpp(a, integer(0), 1L, 1, c(1, 1, 1), TRUE, TRUE),
               "4-D dim attribute")
  a4 <- array(1, c(3L, 3L, 3L, 2L))
  expect_error(ns$gaussian_blur_sep_4d_cpp(a4, 99L, 1L, 1, c(1, 1, 1), TRUE, FALSE),
               "out of bounds")
  expect_error(ns$gaussian_blur_sep_4d_cpp(a4, 1L, 1L, 0, c(1, 1, 1), TRUE, FALSE),
               "positive")
})

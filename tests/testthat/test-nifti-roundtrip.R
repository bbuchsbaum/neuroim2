test_that("NIfTI round-trip preserves sparse voxel ordering, geometry and labels", {
  sp <- NeuroSpace(c(2L, 2L, 2L), spacing = c(2, 3, 4),
                   origin = c(-5, 7, 11))
  ix <- c(1L, 3L, 6L)
  mask_data <- array(FALSE, dim(sp))
  mask_data[ix] <- TRUE
  mask <- LogicalNeuroVol(mask_data, sp)
  paths <- character()
  on.exit(unlink(paths), add = TRUE)

  # Non-square time x masked-voxel matrices also exercise singleton time and
  # the 64-volume chunk boundary in sparse expansion.
  for (nt in c(1L, 4L, 65L)) {
    values <- matrix(seq_len(nt * length(ix)) / 10, nt, length(ix))
    expected <- array(0, c(dim(sp), nt))
    for (t in seq_len(nt)) {
      for (k in seq_along(ix)) {
        expected[ix[k] + (t - 1L) * prod(dim(sp))] <- values[t, k]
      }
    }
    for (labels in list(character(), paste("trial", seq_len(nt)))) {
      vec <- NeuroVec(values, add_dim(sp, nt), mask = mask,
                      volume_labels = labels)
      expect_true(methods::validObject(vec))
      expect_equal(series(vec, ix), values)
      expect_equal(as.array(methods::as(vec, "DenseNeuroVec")), expected)

      for (suffix in c(".nii", ".nii.gz")) {
        path <- tempfile(fileext = suffix)
        paths <- c(paths, path)
        write_vec(vec, path)
        back <- read_vec(path)

        expect_equal(dim(back), dim(vec))
        expect_equal(as.array(back), expected, tolerance = 1e-6)
        expect_equal(as.numeric(series(back, ix)), as.numeric(values),
                     tolerance = 1e-6)
        expect_equal(as.numeric(series(back, setdiff(seq_len(8L), ix))),
                     rep(0, 5L * nt))
        expect_equal(spacing(back), spacing(vec), tolerance = 1e-6)
        expect_equal(origin(back), origin(vec), tolerance = 1e-6)
        expect_equal(trans(space(back)), trans(space(vec)), tolerance = 1e-6)
        expect_equal(volume_labels(back), labels)
      }
    }
  }
})

test_that("NIfTI round-trip preserves DenseNeuroVol data and affine", {
  vol <- make_vol(c(10, 10, 10), spacing = c(2, 2, 2), origin = c(-10, -10, -10))
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)

  expect_equal(as.numeric(vol@.Data), as.numeric(vol2@.Data), tolerance = 1e-5)
  expect_equal(dim(vol), dim(vol2))
  expect_equal(spacing(space(vol)), spacing(space(vol2)), tolerance = 1e-5)
  expect_equal(trans(space(vol)), trans(space(vol2)), tolerance = 1e-5)
})

test_that("NIfTI round-trip preserves DenseNeuroVol with non-unit spacing", {
  vol <- make_vol(c(8, 8, 8), spacing = c(3, 3, 4), origin = c(-12, -12, -16))
  tmp <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp))

  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)

  expect_equal(as.numeric(vol@.Data), as.numeric(vol2@.Data), tolerance = 1e-5)
  expect_equal(spacing(space(vol)), spacing(space(vol2)), tolerance = 1e-5)
})

test_that("NIfTI round-trip preserves sform affine", {
  # Create a volume with oblique affine (sform should be preferred)
  aff <- matrix(c(1.5, 0.1, 0, -80,
                   0.1, 2.0, 0, -120,
                   0, 0, 2.5, -60,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(10, 10, 10), trans = aff)
  vol <- DenseNeuroVol(rnorm(1000), sp)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)

  expect_equal(trans(space(vol)), trans(space(vol2)), tolerance = 1e-4)
})

test_that("NIfTI round-trip works for gzipped vs uncompressed", {
  vol <- make_vol(c(8, 8, 8))

  tmp_gz <- tempfile(fileext = ".nii.gz")
  tmp_nii <- tempfile(fileext = ".nii")
  on.exit(unlink(c(tmp_gz, tmp_nii)))

  write_vol(vol, tmp_gz)
  write_vol(vol, tmp_nii)

  vol_gz <- read_vol(tmp_gz)
  vol_nii <- read_vol(tmp_nii)

  expect_equal(as.numeric(vol_gz@.Data), as.numeric(vol_nii@.Data), tolerance = 1e-6)
  expect_equal(as.numeric(vol@.Data), as.numeric(vol_gz@.Data), tolerance = 1e-5)
})

test_that("NIfTI round-trip preserves DenseNeuroVec data", {
  skip_on_cran()
  vec <- make_vec(c(8, 8, 8), ntime = 5)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vec(vec, tmp)
  vec2 <- read_vec(tmp)

  expect_equal(dim(vec), dim(vec2))
  # Compare raw numeric data, ignoring attributes
  expect_equal(as.numeric(vec@.Data), as.numeric(vec2@.Data), tolerance = 1e-5)
})

test_that("NIfTI round-trip preserves 4D spatial affine", {
  skip_on_cran()
  sp4 <- NeuroSpace(c(8L, 8L, 8L, 5L), spacing = c(2, 2, 2))
  vec <- DenseNeuroVec(array(rnorm(8*8*8*5), c(8,8,8,5)), sp4)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vec(vec, tmp)
  vec2 <- read_vec(tmp)

  # 3D spatial affine should be preserved
  tx1 <- trans(space(vec))
  tx2 <- trans(space(vec2))
  expect_equal(tx1[1:3, 1:3], tx2[1:3, 1:3], tolerance = 1e-5)
})

test_that("NIfTI round-trip preserves LogicalNeuroVol as binary mask", {
  mask <- make_mask(c(10, 10, 10), frac = 0.4)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vol(mask, tmp)
  vol2 <- read_vol(tmp)

  # Logical -> numeric on disk; compare binary values
  expect_equal(as.numeric(mask@.Data), as.numeric(vol2@.Data), tolerance = 1e-6)
})

test_that("NIfTI round-trip preserves data range (no accidental scaling)", {
  # Use a known range to catch slope/intercept issues
  vol <- make_vol(c(8, 8, 8), data = seq_len(512) * 0.1)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))

  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)

  expect_equal(range(vol@.Data), range(vol2@.Data), tolerance = 1e-4)
})

test_that("read_vol with extdata file works", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(f == "", message = "extdata not available")
  vol <- read_vol(f)
  expect_s4_class(vol, "NeuroVol")
  expect_equal(length(dim(vol)), 3)
})

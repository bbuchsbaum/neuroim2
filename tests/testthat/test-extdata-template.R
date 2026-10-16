library(neuroim2)

# Regression tests for #50: inst/extdata/mni_downsampled.nii.gz is the MNI152
# 1 mm template (193 x 229 x 193) downsampled to 48 x 57 x 48, so its voxels
# are ~4.02 mm. Its sform used to be the identity plus a translation, which
# takes precedence over the (correct) qform and made the image read as 1 mm.

template_file <- function() {
  system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2")
}

test_that("bundled template reports ~4.02 mm voxels (#50)", {
  v <- read_vol(template_file())
  expect_equal(dim(v), c(48L, 57L, 48L))
  expect_equal(spacing(v), c(4.02083, 4.01754, 4.02083), tolerance = 1e-5)
  # Field of view matches the 1 mm source grid: 193 x 229 x 193 mm.
  expect_equal(dim(v) * spacing(v), c(193, 229, 193), tolerance = 1e-4)
})

test_that("bundled template's sform agrees with its qform and pixdim (#50)", {
  h <- read_header(template_file())@header
  expect_gt(h$sform_code, 0)
  expect_gt(h$qform_code, 0)
  expect_equal(unname(h$sform), unname(h$qform), tolerance = 1e-5)
  expect_equal(abs(diag(h$sform)[1:3]), h$pixdim[2:4], tolerance = 1e-5)
})

test_that("bundled template's world coordinates are anatomical (#50)", {
  v <- read_vol(template_file())
  # The far corner of the grid sits at the far corner of MNI space, not ~50 mm
  # from the origin as it did with 1 mm voxels.
  far <- grid_to_coord(space(v), matrix(dim(v), nrow = 1))
  expect_equal(as.vector(far), c(-96, -132, -78) + (dim(v) - 1) * spacing(v),
               tolerance = 1e-4)
  expect_true(all(as.vector(far) > c(90, 90, 105)))
  # The MNI origin (anterior commissure) lies inside the grid, and the brain
  # is there: the voxel at (0, 0, 0) mm is well above background.
  ac <- as.vector(coord_to_grid(space(v), matrix(c(0, 0, 0), nrow = 1)))
  expect_true(all(ac >= 1 & ac <= dim(v)))
  expect_gt(v[round(ac[1]), round(ac[2]), round(ac[3])],
            stats::quantile(as.array(v), 0.5))
})

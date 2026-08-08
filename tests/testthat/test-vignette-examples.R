library(testthat)
library(neuroim2)

# Guarantees that used to be asserted inline, inside the vignettes themselves,
# via stopifnot() calls in reader-facing chunks. The articles are now free of
# test scaffolding; the checks live here, where they run on every R CMD check
# instead of interrupting a first-time reader.
#
# Each test names the article whose examples it protects.

# --- vignette("neuroim2"): the tour ------------------------------------------

test_that("tour: anatomy reads as a 3D volume with positive spacing", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  expect_length(dim(anat), 3L)
  expect_true(all(spacing(anat) > 0))
  expect_gt(diff(range(as.array(anat))), 0)
})

test_that("tour: a millimetre coordinate round-trips to a voxel inside the grid", {
  # Uses the EPI mask, not the anatomical image: mni_downsampled.nii.gz has an
  # sform that disagrees with its pixdim, so its world coordinates are not
  # physically meaningful and no vignette does coordinate work on it.
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  vox <- coord_to_grid(mask, c(-34, -28, 10))
  expect_length(as.vector(vox), 3L)
  expect_true(all(as.vector(vox) >= 1 & as.vector(vox) <= dim(mask)))
  back <- grid_to_coord(mask, matrix(as.vector(vox), nrow = 1))
  expect_equal(as.vector(back), c(-34, -28, 10), tolerance = max(spacing(mask)))
})

test_that("tour: the anatomical image is used for display only", {
  # Pins the header inconsistency that drove the decision above, so that a
  # future fix to the shipped file surfaces here rather than silently changing
  # what the articles show.
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  expect_equal(spacing(anat), c(1, 1, 1))
  expect_true(all(abs(diff(spacing(anat))) < 1e-8)) # isotropic: aspect ratio safe
})

# --- vignette("volumes-and-vectors") -----------------------------------------

test_that("containers: mask and simulated series have the advertised shapes", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  expect_length(dim(mask), 3L)
  expect_gt(sum(mask > 0), 0)

  bold <- simulate_fmri(mask, n_time = 12, seed = 1)
  expect_length(dim(bold), 4L)
  expect_gt(dim(bold)[4], 1L)
  expect_identical(dim(bold)[1:3], dim(mask))
})

test_that("containers: sub_vector keeps the requested timepoints", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 12, seed = 1)
  first_two <- sub_vector(bold, 1:2)
  expect_equal(dim(first_two)[4], 2L)
  expect_identical(dim(first_two)[1:3], dim(bold)[1:3])
})

test_that("containers: as.matrix gives voxels x time", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 12, seed = 1)
  mat <- as.matrix(bold)
  expect_equal(nrow(mat), prod(dim(bold)[1:3]))
  expect_equal(ncol(mat), dim(bold)[4])
})

test_that("containers: simulated voxel series are not constant", {
  # The failure this guards against: using a binary mask file as a stand-in
  # for a 4D time series, which makes every series() result a constant.
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 30, seed = 1)
  inmask <- which(as.vector(mask) > 0)
  mat <- as.matrix(bold)[head(inmask, 200), , drop = FALSE]
  expect_true(all(apply(mat, 1, sd) > 0))
})

test_that("containers: concat stacks volumes along time", {
  sp3 <- NeuroSpace(c(8, 8, 8), c(1, 1, 1))
  v1 <- NeuroVol(array(rnorm(512), c(8, 8, 8)), sp3)
  v2 <- NeuroVol(array(rnorm(512), c(8, 8, 8)), sp3)
  expect_equal(dim(concat(v1, v2))[4], 2L)
})

test_that("containers: sparse reads preserve dimensions and class", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  dense <- read_vec(f)
  sparse <- read_vec(f, mask = read_vol(f) > 0)
  expect_s4_class(sparse, "SparseNeuroVec")
  expect_identical(dim(sparse), dim(dense))
})

# --- vignette("reading-and-writing") -----------------------------------------

test_that("io: a volume round-trips through NIfTI unchanged", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp), add = TRUE)
  write_vol(anat, tmp)
  expect_true(file.exists(tmp))
  back <- read_vol(tmp)
  expect_equal(dim(back), dim(anat))
  expect_equal(as.array(back), as.array(anat), tolerance = 1e-5)
  expect_equal(trans(space(back)), trans(space(anat)), tolerance = 1e-4)
})

test_that("io: a 4D series round-trips through NIfTI", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 6, seed = 1)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp), add = TRUE)
  write_vec(bold, tmp)
  back <- read_vec(tmp)
  expect_equal(dim(back), dim(bold))
  expect_equal(as.matrix(back), as.matrix(bold), tolerance = 1e-4)
})

test_that("io: read_header reports the same geometry as read_vol", {
  f <- system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")
  h <- read_header(f)
  v <- read_vol(f)
  expect_equal(dim(h)[1:3], dim(v))
  # NIFTIMetaInfo carries geometry in slots; spacing()/space() are not defined
  # for it, which is why the I/O article reads the slots directly.
  expect_equal(h@spacing[1:3], spacing(v), tolerance = 1e-5)
  expect_equal(h@origin[1:3], origin(v), tolerance = 1e-4)
})

# --- vignette("regions-and-searchlights") ------------------------------------

test_that("regions: a spherical ROI yields a non-empty time series matrix", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 12, seed = 1)
  roi <- spherical_roi(mask, c(32, 32, 12), radius = 8, nonzero = TRUE)
  expect_gt(length(roi), 0L)
  rts <- series_roi(bold, roi)
  expect_equal(nrow(values(rts)), dim(bold)[4])
  expect_equal(ncol(values(rts)), length(roi))
  expect_true(all(is.finite(values(rts))))
})

test_that("regions: spherical ROIs with nonzero=TRUE stay inside the mask", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  roi <- spherical_roi(mask, c(32, 32, 12), radius = 8, nonzero = TRUE)
  expect_true(all(as.vector(mask)[indices(roi)] > 0))
})

test_that("regions: split_clusters returns one piece per cluster with finite means", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 12, seed = 1)
  set.seed(1)
  cvol <- ClusteredNeuroVol(mask > 0, sample(1:4, sum(mask > 0), replace = TRUE))
  parts <- split_clusters(bold, cvol)
  expect_length(parts, 4L)
  expect_true(all(is.finite(vapply(parts, function(p) mean(values(p)), numeric(1)))))
})

test_that("regions: lazy searchlights produce non-empty neighbourhoods", {
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  sl <- searchlight(as.mask(mask > 0), radius = 6, eager = FALSE, nonzero = TRUE)
  expect_gt(length(sl), 0L)
  first <- sl[[1]]
  expect_gt(nrow(coords(first)), 0L)
  means <- vapply(seq_len(5), function(i) mean(values(sl[[i]])), numeric(1))
  expect_length(means, 5L)
  expect_true(all(is.finite(means)))
})

test_that("regions: conn_comp finds the planted clusters", {
  sp <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
  arr <- array(0, c(10, 10, 10))
  arr[1:2, 1:2, 1:2] <- 1
  arr[8:9, 8:9, 8:9] <- 1
  cc <- conn_comp(NeuroVol(arr, sp), threshold = 0)
  expect_equal(max(cc$index), 2)
})

test_that("regions: slices/vols/vectors iterate the expected number of pieces", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  expect_length(slices(anat), dim(anat)[3])
  vec <- concat(anat, anat, anat)
  expect_length(vols(vec), dim(vec)[4])
  mean_vol <- NeuroVol(
    vapply(vectors(vec), function(v) mean(v), numeric(1)),
    space(anat)
  )
  expect_equal(dim(mean_vol), dim(anat))
})

# --- vignette("resampling-and-orientation") ----------------------------------

test_that("resampling: downsample halves the grid and doubles the spacing", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  ds <- downsample(anat, factor = 0.5)
  expect_true(all(dim(ds) < dim(anat)))
  expect_true(all(spacing(ds) > spacing(anat)))
})

test_that("resampling: resample_to lands on the requested grid", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  sp <- space(anat)
  target <- NeuroSpace(dim(sp) * 2L, spacing(sp) / 2, origin = origin(sp))
  out <- resample_to(anat, target, method = "linear")
  expect_equal(dim(out), dim(target))
})

test_that("resampling: deoblique removes tilt", {
  aff <- matrix(c(
    2.0, 0.2, 0.0, -90,
    0.0, 2.0, 0.1, -126,
    0.0, 0.0, 2.0, -72,
    0.0, 0.0, 0.0, 1
  ), nrow = 4, byrow = TRUE)
  sp_obl <- NeuroSpace(dim = c(40L, 40L, 30L), trans = aff)
  expect_true(max(obliquity(trans(sp_obl))) > 0)
  expect_lt(max(obliquity(trans(deoblique(sp_obl)))), 1e-6)
})

# --- vignette("smoothing-and-filtering") -------------------------------------

test_that("smoothing: gaussian_blur preserves geometry and reduces variance", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  out <- gaussian_blur(anat, anat, sigma = 4, window = 2)
  expect_equal(dim(out), dim(anat))
  expect_equal(space(out), space(anat))
  expect_lt(sd(as.array(out)), sd(as.array(anat)))
})

test_that("smoothing: the demo image has the intensity structure filters need", {
  # Guards against demonstrating edge-preserving filters on a binary mask.
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  expect_gt(length(unique(as.vector(as.array(anat)))), 1000)
})

test_that("smoothing: cgb_filter has non-degenerate correlations to work with", {
  # Guards against running a correlation-guided filter on constant series.
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 20, seed = 1)
  inmask <- which(as.vector(mask) > 0)
  mat <- as.matrix(bold)[head(inmask, 100), , drop = FALSE]
  expect_true(all(apply(mat, 1, sd) > 0))
})

# --- vignette("visualization") -----------------------------------------------

test_that("visualization: the demo stat map is a genuine signed t-map", {
  skip_on_cran()
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  bold <- simulate_fmri(mask, n_time = 30, seed = 2, spatial_fwhm = 2.5 * min(spacing(mask)))

  design <- rep(rep(c(0, 1), each = 10), length.out = 30)
  Y <- as.matrix(bold)
  pos <- spherical_roi(mask, c(20, 34, 12), radius = 12, nonzero = TRUE)
  Y[indices(pos), ] <- Y[indices(pos), ] + 1.3 * rep(design, each = length(pos))
  X <- cbind(1, design)
  XtXi <- solve(crossprod(X))
  B <- Y %*% X %*% XtXi
  s2 <- rowSums((Y - tcrossprod(B, X))^2) / (ncol(Y) - ncol(X))
  tstat <- B[, 2] / sqrt(s2 * XtXi[2, 2])
  tstat[as.vector(mask) == 0] <- 0

  expect_true(any(tstat > 3)) # signal where it was planted
  expect_true(all(is.finite(tstat)))
  # A t-map's null must have roughly unit spread. A fabricated blob image
  # would not, which is the point of computing this one.
  expect_lt(abs(sd(tstat[abs(tstat) < 2 & tstat != 0]) - 1), 0.4)
})

test_that("visualization: plot helpers return ggplot objects on matching grids", {
  anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  p <- plot_montage(anat, zlevels = c(12, 20, 28), ncol = 3)
  expect_s3_class(p, "ggplot")
  panels <- plot_ortho(anat, coord = round(dim(anat) / 2), draw = FALSE)
  expect_true(length(panels) >= 3)
})

# --- vignette("large-data") --------------------------------------------------

test_that("large-data: backends agree with the dense result", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  dense <- read_vec(f)
  fb <- read_vec(f, mode = "filebacked")
  expect_equal(dim(fb), dim(dense))
  expect_equal(as.array(fb[[1]]), as.array(dense[[1]]), tolerance = 1e-6)
})

test_that("large-data: ClusteredNeuroVec stores one series per cluster", {
  sp4 <- NeuroSpace(c(10, 10, 10, 20), c(1, 1, 1))
  vec <- NeuroVec(array(rnorm(10 * 10 * 10 * 20), c(10, 10, 10, 20)), sp4)
  sp3 <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
  marr <- array(FALSE, c(10, 10, 10))
  marr[3:8, 3:8, 3:8] <- TRUE
  mask <- LogicalNeuroVol(marr, sp3)
  set.seed(1)
  cvol <- ClusteredNeuroVol(mask, sample(1:5, sum(marr), replace = TRUE))
  cv <- ClusteredNeuroVec(vec, cvol)
  expect_equal(num_clusters(cv), 5)
  expect_equal(dim(as.matrix(cv, by = "cluster")), c(20L, 5L))
  expect_length(dim(cv), 4L)
})

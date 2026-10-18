library(testthat)
library(neuroim2)

# mmap and bigstatsr are Suggests: neuroim2 must work without them (e.g. webR).

test_that("masked reads fall back to streaming when mmap is unavailable", {
  fname <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")) > 0
  ref <- read_vec(fname, mask = mask)

  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  streamed <- read_vec(fname, mask = mask)
  expect_s4_class(streamed, "SparseNeuroVec")
  expect_equal(as.matrix(streamed), as.matrix(ref))
})

test_that("memory-mapped modes fail informatively when mmap is unavailable", {
  fname <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  expect_error(read_vec(fname, mode = "mmap"), "Package 'mmap' is required")
  expect_error(MappedNeuroVec(fname), "Package 'mmap' is required")
})

test_that("BigNeuroVec validity requires an FBM in the data slot", {
  sp <- NeuroSpace(c(2, 2, 1, 2))
  mask <- LogicalNeuroVol(array(TRUE, c(2, 2, 1)), drop_dim(sp))
  expect_error(
    new("BigNeuroVec", space = sp, mask = mask,
        map = IndexLookupVol(space(mask), 1:4), data = matrix(0, 2, 4)),
    "FBM"
  )
})

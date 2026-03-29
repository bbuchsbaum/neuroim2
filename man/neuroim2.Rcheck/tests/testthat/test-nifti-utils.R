context("nifti utils")

library(neuroim2)

# verify helper functions for nifti metadata

test_that("internal nifti helpers return expected values", {
  expect_equal(neuroim2:::.getDataCode("FLOAT"), 16)
  expect_equal(neuroim2:::.getDataSize("DOUBLE"), 8)
  expect_equal(neuroim2:::.getDataStorage(8), "INT")
  expect_equal(neuroim2:::.niftiExt("nifti-gz"), list(header="nii.gz", data="nii.gz"))
})

# check error handling

test_that("nifti helper functions throw errors on invalid input", {
  expect_error(neuroim2:::.getDataCode("FOO"))
  expect_error(neuroim2:::.niftiExt("bogus"))
})

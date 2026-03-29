context("file_format coverage")

library(neuroim2)

# ---- FileFormat construction ----

test_that("NIFTIFormat object exists and has correct slots", {
  expect_s4_class(neuroim2:::NIFTI, "NIFTIFormat")
  expect_equal(neuroim2:::NIFTI@header_extension, "nii")
  expect_equal(neuroim2:::NIFTI@data_extension,   "nii")
  expect_equal(neuroim2:::NIFTI@file_format,      "NIFTI")
})

test_that("NIFTI_GZ object has gzip encoding", {
  expect_s4_class(neuroim2:::NIFTI_GZ, "NIFTIFormat")
  expect_equal(neuroim2:::NIFTI_GZ@header_extension, "nii.gz")
  expect_equal(neuroim2:::NIFTI_GZ@data_extension,   "nii.gz")
})

test_that("NIFTI_PAIR object has hdr/img extensions", {
  expect_s4_class(neuroim2:::NIFTI_PAIR, "NIFTIFormat")
  expect_equal(neuroim2:::NIFTI_PAIR@header_extension, "hdr")
  expect_equal(neuroim2:::NIFTI_PAIR@data_extension,   "img")
})

test_that("AFNI object has correct extensions", {
  expect_s4_class(neuroim2:::AFNI, "AFNIFormat")
  expect_equal(neuroim2:::AFNI@header_extension, "HEAD")
  expect_equal(neuroim2:::AFNI@data_extension,   "BRIK")
})

# ---- header_file_matches ----

test_that("header_file_matches is TRUE for matching extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_true(header_file_matches(fmt, "brain.hdr"))
})

test_that("header_file_matches is FALSE for wrong extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_false(header_file_matches(fmt, "brain.img"))
  expect_false(header_file_matches(fmt, "brain.hdr.gz"))
})

test_that("header_file_matches is TRUE for .nii file with NIFTIFormat", {
  fmt <- neuroim2:::NIFTI
  expect_true(header_file_matches(fmt, "brain.nii"))
})

# ---- data_file_matches ----

test_that("data_file_matches is TRUE for matching extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_true(data_file_matches(fmt, "brain.img"))
})

test_that("data_file_matches is FALSE for wrong extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_false(data_file_matches(fmt, "brain.hdr"))
})

# ---- header_file / data_file ----

test_that("header_file returns same name when already a header file", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(header_file(fmt, "brain.hdr"), "brain.hdr")
})

test_that("header_file derives header name from data file", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(header_file(fmt, "brain.img"), "brain.hdr")
})

test_that("data_file returns same name when already a data file", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(data_file(fmt, "brain.img"), "brain.img")
})

test_that("data_file derives data name from header file", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(data_file(fmt, "brain.hdr"), "brain.img")
})

test_that("header_file errors on unrecognized extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_error(header_file(fmt, "brain.xyz"))
})

test_that("data_file errors on unrecognized extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_error(data_file(fmt, "brain.xyz"))
})

# ---- strip_extension ----

test_that("strip_extension removes header extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(strip_extension(fmt, "brain.hdr"), "brain")
})

test_that("strip_extension removes data extension", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_equal(strip_extension(fmt, "brain.img"), "brain")
})

test_that("strip_extension errors on unrecognized file", {
  fmt <- neuroim2:::NIFTI_PAIR
  expect_error(strip_extension(fmt, "brain.xyz"))
})

# ---- find_descriptor ----

test_that("find_descriptor returns NIFTI for a real .nii file", {
  tmp <- tempfile(fileext = ".nii")
  file.create(tmp)
  on.exit(unlink(tmp))
  result <- neuroim2:::find_descriptor(tmp)
  expect_s4_class(result, "NIFTIFormat")
})

test_that("find_descriptor returns NULL for unknown extension", {
  tmp <- tempfile(fileext = ".xyz")
  file.create(tmp)
  on.exit(unlink(tmp))
  result <- neuroim2:::find_descriptor(tmp)
  expect_null(result)
})

test_that("find_descriptor returns NIFTI_GZ for .nii.gz file", {
  tmp <- tempfile(fileext = ".nii.gz")
  file.create(tmp)
  on.exit(unlink(tmp))
  result <- neuroim2:::find_descriptor(tmp)
  expect_s4_class(result, "NIFTIFormat")
})

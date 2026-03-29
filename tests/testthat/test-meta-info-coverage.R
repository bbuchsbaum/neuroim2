library(testthat)
library(neuroim2)

nii_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")

# ---------------------------------------------------------------------------
# read_header
# ---------------------------------------------------------------------------

test_that("read_header returns a FileMetaInfo object", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  expect_s4_class(hdr, "FileMetaInfo")
  expect_s4_class(hdr, "NIFTIMetaInfo")
})

test_that("read_header dim() returns numeric vector of length >= 3", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  d <- dim(hdr)
  expect_true(is.numeric(d))
  expect_true(length(d) >= 3)
})

test_that("read_header errors on nonexistent file", {
  expect_error(read_header("/nonexistent/path/file.nii"))
})

# ---------------------------------------------------------------------------
# show,FileMetaInfo
# ---------------------------------------------------------------------------

test_that("show method on FileMetaInfo produces output", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  out <- capture.output(show(hdr))
  expect_true(length(out) > 0)
  expect_true(any(grepl("header_file|NIFTIMetaInfo", out)))
})

# ---------------------------------------------------------------------------
# trans() on MetaInfo and NIFTIMetaInfo
# ---------------------------------------------------------------------------

test_that("trans() on NIFTIMetaInfo returns 4x4 matrix", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  tx <- trans(hdr)
  expect_true(is.matrix(tx))
  expect_equal(dim(tx), c(4L, 4L))
})

test_that("trans() on MetaInfo returns correct-sized matrix", {
  mi <- MetaInfo(Dim = c(10L, 10L, 10L), spacing = c(2, 2, 2), origin = c(1, 2, 3))
  tx <- trans(mi)
  expect_true(is.matrix(tx))
  # 3 spatial dims + 1 homogeneous row → 4x4
  expect_equal(dim(tx), c(4L, 4L))
  # origin should appear in the last column
  expect_equal(tx[1:3, 4], c(1, 2, 3))
})

# ---------------------------------------------------------------------------
# MetaInfo constructor
# ---------------------------------------------------------------------------

test_that("MetaInfo constructor creates valid object", {
  mi <- MetaInfo(Dim = c(64L, 64L, 32L), spacing = c(1, 1, 2))
  expect_s4_class(mi, "MetaInfo")
  expect_equal(mi@dims, c(64L, 64L, 32L))
})

test_that("MetaInfo rejects non-positive dims", {
  expect_error(MetaInfo(Dim = c(-1L, 10L, 10L), spacing = c(1, 1, 1)))
})

test_that("MetaInfo rejects non-positive spacing", {
  expect_error(MetaInfo(Dim = c(10L, 10L, 10L), spacing = c(-1, 1, 1)))
})

test_that("MetaInfo rejects non-finite origin", {
  expect_error(MetaInfo(Dim = c(10L, 10L, 10L), spacing = c(1, 1, 1),
                        origin = c(Inf, 0, 0)))
})

test_that("MetaInfo rejects unsupported data_type", {
  expect_error(MetaInfo(Dim = c(10L, 10L, 10L), spacing = c(1, 1, 1),
                        data_type = "UINT32"))
})

# ---------------------------------------------------------------------------
# NIFTIMetaInfo constructor
# ---------------------------------------------------------------------------

test_that("NIFTIMetaInfo constructor errors on bad descriptor type", {
  expect_error(NIFTIMetaInfo("not_a_descriptor", list(file_type = "nifti")))
})

test_that("NIFTIMetaInfo constructor errors on non-list header", {
  fmt <- new("NIFTIFormat",
             file_format = "NIFTI",
             header_encoding = "raw",
             header_extension = "nii",
             data_encoding = "raw",
             data_extension = "nii")
  expect_error(NIFTIMetaInfo(fmt, "not_a_list"))
})

# ---------------------------------------------------------------------------
# meta_info(character) — file path
# ---------------------------------------------------------------------------

test_that("meta_info from file path returns list with expected fields", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_true(is.list(mi))
  expected_fields <- c("dim", "spacing", "origin", "trans", "path",
                       "filename", "format", "dtype", "bytes_per_element",
                       "nvox", "nvol", "size_bytes", "time_step")
  expect_true(all(expected_fields %in% names(mi)))
})

test_that("meta_info dim field is integer vector of length >= 3", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_true(is.integer(mi$dim))
  expect_true(length(mi$dim) >= 3)
})

test_that("meta_info trans field is 4x4 matrix", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_true(is.matrix(mi$trans))
  expect_equal(dim(mi$trans), c(4L, 4L))
})

test_that("meta_info nvox equals product of first 3 dims", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_equal(mi$nvox, as.integer(prod(mi$dim[1:3])))
})

test_that("meta_info nvol matches 4th dim", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  if (length(mi$dim) >= 4) {
    expect_equal(mi$nvol, as.integer(mi$dim[4]))
  } else {
    expect_equal(mi$nvol, 1L)
  }
})

test_that("meta_info size_bytes is nvox * nvol * bytes_per_element", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expected <- as.numeric(mi$nvox) * as.numeric(mi$nvol) * as.numeric(mi$bytes_per_element)
  expect_equal(mi$size_bytes, expected)
})

test_that("meta_info filename is basename of path", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_equal(mi$filename, basename(mi$path))
})

test_that("meta_info format field is non-empty string", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_true(is.character(mi$format))
  expect_true(nchar(mi$format) > 0)
})

test_that("meta_info time_step is numeric for NIfTI", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  mi <- meta_info(nii_path)
  expect_true(is.numeric(mi$time_step))
})

test_that("meta_info vectorised over multiple paths returns a list", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  result <- meta_info(c(nii_path, nii_path))
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(all(sapply(result, is.list)))
})

# ---------------------------------------------------------------------------
# meta_info(FileMetaInfo)
# ---------------------------------------------------------------------------

test_that("meta_info from FileMetaInfo returns identical result to path version", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  mi_from_obj  <- meta_info(hdr)
  mi_from_path <- meta_info(nii_path)
  expect_equal(mi_from_obj$dim,     mi_from_path$dim)
  expect_equal(mi_from_obj$spacing, mi_from_path$spacing)
  expect_equal(mi_from_obj$origin,  mi_from_path$origin)
  expect_equal(mi_from_obj$nvox,    mi_from_path$nvox)
  expect_equal(mi_from_obj$nvol,    mi_from_path$nvol)
})

# ---------------------------------------------------------------------------
# .data_scale_params
# ---------------------------------------------------------------------------

test_that(".data_scale_params returns slope and intercept fields", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  pars <- neuroim2:::.data_scale_params(hdr, 1L)
  expect_true(is.list(pars))
  expect_true("slope" %in% names(pars))
  expect_true("intercept" %in% names(pars))
  expect_true(is.numeric(pars$slope))
  expect_true(is.numeric(pars$intercept))
})

test_that(".data_scale_params slope is never exactly zero (NIfTI identity rule)", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  pars <- neuroim2:::.data_scale_params(hdr, 1L)
  # NIfTI spec: slope==0 is treated as 1 (identity)
  expect_true(pars$slope != 0)
})

test_that(".data_scale_params errors on bad index", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  expect_error(neuroim2:::.data_scale_params(hdr, 0L))
  expect_error(neuroim2:::.data_scale_params(hdr, NA_integer_))
})

test_that(".data_scale_params on MetaInfo (no slope slot) returns defaults", {
  mi <- MetaInfo(Dim = c(10L, 10L, 10L), spacing = c(1, 1, 1))
  # MetaInfo has no slope/intercept slots; function should fall back to 1/0
  pars <- neuroim2:::.data_scale_params(mi, 1L)
  expect_equal(pars$slope, 1)
  expect_equal(pars$intercept, 0)
})

# ---------------------------------------------------------------------------
# .apply_data_scaling
# ---------------------------------------------------------------------------

test_that(".apply_data_scaling scales data correctly", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  raw <- c(1, 2, 3, 4, 5)
  pars <- neuroim2:::.data_scale_params(hdr, 1L)
  scaled <- neuroim2:::.apply_data_scaling(raw, hdr, 1L)
  expected <- as.numeric(raw) * pars$slope + pars$intercept
  expect_equal(scaled, expected)
})

test_that(".apply_data_scaling returns numeric", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  result <- neuroim2:::.apply_data_scaling(1:10, hdr, 1L)
  expect_true(is.numeric(result))
  expect_equal(length(result), 10)
})

# ---------------------------------------------------------------------------
# NIFTIMetaInfo specific fields via read_header
# ---------------------------------------------------------------------------

test_that("NIFTIMetaInfo slope and intercept slots exist", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  expect_true(.hasSlot(hdr, "slope"))
  expect_true(.hasSlot(hdr, "intercept"))
})

test_that("NIFTIMetaInfo header slot is a list", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  expect_true(is.list(hdr@header))
})

test_that("NIFTIMetaInfo spacing is positive", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  expect_true(all(hdr@spacing > 0))
})

test_that("data_reader for NIFTIMetaInfo returns BinaryReader", {
  skip_if(nii_path == "", "extdata NIfTI file not available")
  hdr <- read_header(nii_path)
  reader <- data_reader(hdr, offset = 0)
  expect_s4_class(reader, "BinaryReader")
  close(reader)
})

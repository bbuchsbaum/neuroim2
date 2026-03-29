context("header API")

# ---- header(character) ----

test_that("header(character) returns a NeuroHeader", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_s3_class(h, "NeuroHeader")
})

test_that("header(character) dim field is an integer vector", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.integer(h$dim))
  expect_true(length(h$dim) >= 3)
})

test_that("header(character) spacing is numeric with length >= 3", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.numeric(h$spacing))
  expect_true(length(h$spacing) >= 3)
  expect_true(all(h$spacing > 0))
})

test_that("header(character) origin is numeric", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.numeric(h$origin))
})

test_that("header(character) trans is a 4x4 matrix", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.matrix(h$trans))
  expect_equal(dim(h$trans), c(4, 4))
})

test_that("header(character) sform is a list with matrix and code", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.list(h$sform))
  expect_true("matrix" %in% names(h$sform))
  expect_true("code" %in% names(h$sform))
})

test_that("header(character) qform is a list with matrix and code", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.list(h$qform))
  expect_true("matrix" %in% names(h$qform))
  expect_true("code" %in% names(h$qform))
})

test_that("header(character) data_type is character", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.character(h$data_type) || is.na(h$data_type))
})

test_that("header(character) raw is a list", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.list(h$raw))
})

test_that("header(character) TR field exists", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true("TR" %in% names(h))
})

test_that("header(character) has all required fields", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  required <- c("dim", "pixdim", "spacing", "origin", "trans",
                "qform", "sform", "intent_code", "intent_name",
                "descrip", "data_type", "bitpix", "scl_slope",
                "scl_inter", "cal_min", "cal_max", "TR", "raw")
  expect_true(all(required %in% names(h)))
})

# ---- header(FileMetaInfo) ----

test_that("header(FileMetaInfo) returns a NeuroHeader", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  mi <- read_header(f)
  h <- header(mi)
  expect_s3_class(h, "NeuroHeader")
})

test_that("header(FileMetaInfo) and header(character) produce same dim", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h_char <- header(f)
  mi <- read_header(f)
  h_mi <- header(mi)
  expect_equal(h_char$dim, h_mi$dim)
})

test_that("header(FileMetaInfo) and header(character) produce same spacing", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h_char <- header(f)
  mi <- read_header(f)
  h_mi <- header(mi)
  expect_equal(h_char$spacing, h_mi$spacing)
})

# ---- print.NeuroHeader ----

test_that("print.NeuroHeader produces output without error", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  out <- capture.output(print(h))
  expect_true(length(out) > 0)
})

test_that("print.NeuroHeader output contains NIfTI", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  out <- capture.output(print(h))
  expect_true(any(grepl("NIfTI", out)))
})

test_that("print.NeuroHeader returns invisibly", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  out <- capture.output(result <- print(h))
  expect_s3_class(result, "NeuroHeader")
})

# ---- 4D file TR ----

test_that("header TR is NA or numeric for 3D file", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(is.numeric(h$TR) || is.na(h$TR))
})

context("nifti extensions additional coverage")

library(neuroim2)

# ---- NiftiExtension with raw data ----

test_that("NiftiExtension errors on non-character non-raw data", {
  expect_error(NiftiExtension(ecode = 6L, data = 42), "character string or raw vector")
})

test_that("NiftiExtension with raw vector input works", {
  rd <- as.raw(c(0x41, 0x42, 0x43, 0x00))  # "ABC\0"
  ext <- NiftiExtension(ecode = 6L, data = rd)
  expect_s4_class(ext, "NiftiExtension")
  expect_equal(ext@esize %% 16L, 0L)
})

# ---- parse_extension unknown code ----

test_that("parse_extension returns raw and warns for unknown ecode", {
  ext <- NiftiExtension(ecode = 99L, data = "some data")
  expect_warning(result <- parse_extension(ext), "Unknown extension code")
  expect_true(is.raw(result))
})

test_that("parse_extension returns raw data unchanged for known non-text code", {
  # ecode 32 = CIFTI — not 4 or 6, should return raw without warning
  ext <- NiftiExtension(ecode = 32L, data = charToRaw("cifti"))
  result <- parse_extension(ext)
  expect_true(is.raw(result))
})

# ---- parse_afni_extension ----

test_that("parse_afni_extension errors when ecode != 4", {
  ext <- NiftiExtension(ecode = 6L, data = "comment")
  expect_error(parse_afni_extension(ext), "ecode != 4")
})

test_that("parse_afni_extension with as_xml=FALSE returns character", {
  afni_xml <- '<?xml version="1.0" ?><AFNI_attributes></AFNI_attributes>'
  ext <- NiftiExtension(ecode = 4L, data = afni_xml)
  result <- parse_afni_extension(ext, as_xml = FALSE)
  expect_type(result, "character")
  expect_true(nchar(result) > 0)
})

# ---- show methods ----

test_that("show(NiftiExtension) produces output", {
  ext <- NiftiExtension(ecode = 6L, data = "hello world")
  out <- capture.output(show(ext))
  expect_true(any(grepl("NIfTI Extension", out)))
  expect_true(any(grepl("ecode", out)))
})

test_that("show(NiftiExtension) for AFNI ecode shows preview", {
  afni_xml <- '<?xml version="1.0" ?><AFNI_attributes></AFNI_attributes>'
  ext <- NiftiExtension(ecode = 4L, data = afni_xml)
  out <- capture.output(show(ext))
  expect_true(any(grepl("preview", out)))
})

test_that("show(NiftiExtension) truncates long previews", {
  long_str <- paste(rep("x", 100), collapse = "")
  ext <- NiftiExtension(ecode = 6L, data = long_str)
  out <- capture.output(show(ext))
  expect_true(any(grepl("\\.\\.\\.", out)))
})

test_that("show(NiftiExtensionList) empty list produces output", {
  el <- new("NiftiExtensionList")
  out <- capture.output(show(el))
  expect_true(any(grepl("0 extension", out)))
})

test_that("show(NiftiExtensionList) with extensions shows count and bytes", {
  ext1 <- NiftiExtension(ecode = 6L, data = "test")
  ext2 <- NiftiExtension(ecode = 4L, data = "more")
  el <- new("NiftiExtensionList", list(ext1, ext2))
  out <- capture.output(show(el))
  expect_true(any(grepl("2 extension", out)))
  expect_true(any(grepl("Total size", out)))
})

# ---- has_extensions ----

test_that("has_extensions on list with NULL extensions returns FALSE", {
  lst <- list(other = "data")
  expect_false(has_extensions(lst))
})

test_that("has_extensions on list with non-empty extensions returns TRUE", {
  ext <- NiftiExtension(ecode = 6L, data = "x")
  lst <- list(extensions = new("NiftiExtensionList", list(ext)))
  expect_true(has_extensions(lst))
})

# ---- extension filter on empty list ----

test_that("extension filter on empty NiftiExtensionList returns empty", {
  el <- new("NiftiExtensionList")
  result <- extension(el, 6L)
  expect_equal(length(result), 0L)
})

# ---- total_extension_size ----

test_that("total_extension_size with empty list returns 4L", {
  el <- new("NiftiExtensionList")
  expect_equal(neuroim2:::total_extension_size(el), 4L)
})

test_that("total_extension_size with NULL returns 4L", {
  expect_equal(neuroim2:::total_extension_size(NULL), 4L)
})

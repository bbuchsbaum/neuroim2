# Tests for NIfTI extension support

test_that("NiftiExtension constructor works correctly", {
  # Test with character data
  ext <- NiftiExtension(ecode = 6L, data = "test comment")
  expect_s4_class(ext, "NiftiExtension")
  expect_equal(ext@ecode, 6L)
  expect_true(ext@esize >= 16)
  expect_equal(ext@esize %% 16, 0)  # Must be multiple of 16

  # Test with raw data
  raw_data <- charToRaw("raw test")
  ext2 <- NiftiExtension(ecode = 6L, data = raw_data)
  expect_s4_class(ext2, "NiftiExtension")
})

test_that("NiftiExtension padding is correct", {
  # Test various string lengths to verify padding
  for (len in c(1, 5, 10, 15, 20, 30)) {
    data <- paste(rep("x", len), collapse = "")
    ext <- NiftiExtension(ecode = 6L, data = data)
    expect_equal(ext@esize %% 16, 0, info = paste("Length", len))
    expect_equal(length(ext@edata), ext@esize - 8L)
  }
})

test_that("NiftiExtensionList works correctly", {
  # Empty list
  ext_list <- new("NiftiExtensionList")
  expect_s4_class(ext_list, "NiftiExtensionList")
  expect_equal(length(ext_list), 0)

  # List with extensions
  ext1 <- NiftiExtension(ecode = 6L, data = "comment 1")
  ext2 <- NiftiExtension(ecode = 6L, data = "comment 2")
  ext_list <- new("NiftiExtensionList", list(ext1, ext2))
  expect_equal(length(ext_list), 2)
})

test_that("NiftiExtensionList validation works", {
  # Should fail with non-NiftiExtension objects
  expect_error(
    new("NiftiExtensionList", list("not an extension")),
    "NiftiExtension"
  )
})

test_that("ecode_name returns correct names", {
  expect_equal(ecode_name(4L), "AFNI")
  expect_equal(ecode_name(6L), "comment")
  expect_equal(ecode_name(32L), "CIFTI")
  expect_equal(ecode_name(999L), "unknown")
})

test_that("parse_extension works for comments", {
  ext <- NiftiExtension(ecode = 6L, data = "This is a test")
  parsed <- parse_extension(ext)
  expect_equal(parsed, "This is a test")
})

test_that("parse_extension works for AFNI XML", {
  afni_xml <- '<?xml version="1.0" ?><AFNI_attributes></AFNI_attributes>'
  ext <- NiftiExtension(ecode = 4L, data = afni_xml)
  parsed <- parse_extension(ext)

  if (requireNamespace("xml2", quietly = TRUE)) {
    expect_s3_class(parsed, "xml_document")
  } else {
    expect_type(parsed, "character")
  }
})

test_that("total_extension_size calculates correctly", {
  # Empty list: just extender (4 bytes)
  expect_equal(neuroim2:::total_extension_size(NULL), 4L)
  expect_equal(neuroim2:::total_extension_size(new("NiftiExtensionList")), 4L)

  # With extensions
  ext1 <- NiftiExtension(ecode = 6L, data = "test")
  ext_list <- new("NiftiExtensionList", list(ext1))
  expected_size <- 4L + ext1@esize
  expect_equal(neuroim2:::total_extension_size(ext_list), expected_size)
})

test_that("extension method filters by ecode", {
  ext1 <- NiftiExtension(ecode = 6L, data = "comment")
  ext2 <- NiftiExtension(ecode = 4L, data = "<AFNI/>")
  ext3 <- NiftiExtension(ecode = 6L, data = "another comment")

  ext_list <- new("NiftiExtensionList", list(ext1, ext2, ext3))

  # Filter for comments
  comments <- extension(ext_list, 6L)
  expect_equal(length(comments), 2)

  # Filter for AFNI
  afni <- extension(ext_list, 4L)
  expect_equal(length(afni), 1)

  # Filter for non-existent
  none <- extension(ext_list, 999L)
  expect_equal(length(none), 0)
})

test_that("has_extensions works correctly", {
  # Empty list
  ext_list <- new("NiftiExtensionList")
  expect_false(has_extensions(ext_list))

  # Non-empty list
  ext <- NiftiExtension(ecode = 6L, data = "test")
  ext_list <- new("NiftiExtensionList", list(ext))
  expect_true(has_extensions(ext_list))

  # List with extensions field
  header <- list(extensions = ext_list)
  expect_true(has_extensions(header))

  # List without extensions
  header <- list(other = "data")
  expect_false(has_extensions(header))
})

test_that("as_nifti_header includes extensions correctly", {
  # Create a simple volume
  vol_data <- array(rnorm(8 * 8 * 8), dim = c(8, 8, 8))
  space <- NeuroSpace(dim = c(8L, 8L, 8L), origin = c(0, 0, 0), spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(vol_data, space)

  # Without extensions
  hdr1 <- as_nifti_header(vol, file_name = "test.nii")
  expect_equal(hdr1$vox_offset, 352)  # 348 + 4 byte extender
  expect_equal(length(hdr1$extensions), 0)

  # With extensions
  ext <- NiftiExtension(ecode = 6L, data = "test")
  ext_list <- new("NiftiExtensionList", list(ext))
  hdr2 <- as_nifti_header(vol, file_name = "test.nii", extensions = ext_list)

  expect_gt(hdr2$vox_offset, 352)
  expect_equal(length(hdr2$extensions), 1)
  expect_equal(hdr2$vox_offset, 348 + neuroim2:::total_extension_size(ext_list))
})

test_that("extensions round-trip through file I/O", {
  skip_on_cran()

  # Create volume
  vol_data <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  space <- NeuroSpace(dim = c(4L, 4L, 4L), origin = c(0, 0, 0), spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(vol_data, space)

  # Create extensions
  ext1 <- NiftiExtension(ecode = 6L, data = "Round trip test")
  ext2 <- NiftiExtension(ecode = 6L, data = "Second extension")
  ext_list <- new("NiftiExtensionList", list(ext1, ext2))

  # Create header with extensions
  tmp_file <- tempfile(fileext = ".nii")
  hdr <- as_nifti_header(vol, file_name = tmp_file, extensions = ext_list)

  # Write file
  conn <- file(tmp_file, open = "wb")
  write_nifti_header(hdr, conn, close = FALSE)
  writer <- BinaryWriter(conn, hdr$vox_offset, "FLOAT", 4L, .Platform$endian)
  write_elements(writer, as.numeric(vol))
  close(writer)

  # Read back
  hdr_read <- read_nifti_header(tmp_file)

  # Verify extensions
  expect_equal(length(hdr_read$extensions), 2)
  expect_equal(hdr_read$extensions[[1]]@ecode, 6L)
  expect_equal(hdr_read$extensions[[2]]@ecode, 6L)

  # Parse and verify content
  content1 <- parse_extension(hdr_read$extensions[[1]])
  content2 <- parse_extension(hdr_read$extensions[[2]])
  expect_equal(content1, "Round trip test")
  expect_equal(content2, "Second extension")

  # Clean up
  unlink(tmp_file)
})

test_that("AFNI extension attributes can be extracted", {
  skip_if_not_installed("xml2")

  afni_xml <- paste0(
    '<?xml version="1.0" ?>',
    '<AFNI_attributes self_idcode="TEST123">',
    '<AFNI_atr atr_name="HISTORY_NOTE" ni_type="String">Test history</AFNI_atr>',
    '<AFNI_atr atr_name="BRICK_LABS" ni_type="String">Vol1~Vol2</AFNI_atr>',
    '</AFNI_attributes>'
  )

  ext <- NiftiExtension(ecode = 4L, data = afni_xml)

  # List attributes
  attrs <- list_afni_attributes(ext)
  expect_true("HISTORY_NOTE" %in% attrs)
  expect_true("BRICK_LABS" %in% attrs)

  # Get specific attribute
  history <- get_afni_attribute(ext, "HISTORY_NOTE")
  expect_equal(trimws(history), "Test history")
})

test_that("NiftiExtensionCodes contains expected codes", {
  expect_true("AFNI" %in% names(NiftiExtensionCodes))
  expect_true("comment" %in% names(NiftiExtensionCodes))
  expect_true("CIFTI" %in% names(NiftiExtensionCodes))
  expect_true("DICOM" %in% names(NiftiExtensionCodes))

  expect_equal(NiftiExtensionCodes["AFNI"], c(AFNI = 4L))
  expect_equal(NiftiExtensionCodes["comment"], c(comment = 6L))
})

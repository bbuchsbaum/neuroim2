# AFNI stores geometry in DICOM (RAI) coordinates: +x = Left, +y = Posterior.
# The voxel-to-world transform must be converted to RAS (negate x and y), as
# nibabel's AFNI reader does; previously the DICOM origin was used as-is and
# axis flips were dropped, mirroring datasets in x and y.

write_afni_pair <- function(dir, name, attrs, dims = c(4L, 5L, 6L)) {
  block <- function(type, nm, vals) {
    vals <- if (type == "string") paste0("'", vals, "~") else paste(vals, collapse = " ")
    count <- if (type == "string") nchar(vals) - 1L else length(strsplit(vals, " ")[[1]])
    c("", paste0("type = ", type, "-attribute"), paste0("name = ", nm),
      paste0("count = ", count), vals)
  }
  lines <- c(
    block("integer", "DATASET_DIMENSIONS", c(dims, 0L, 0L)),
    block("integer", "DATASET_RANK", c(3L, 1L)),
    block("integer", "BRICK_TYPES", 3L),
    block("float", "BRICK_FLOAT_FACS", 0),
    block("string", "BYTEORDER_STRING", "LSB_FIRST")
  )
  for (nm in names(attrs)) {
    type <- if (nm == "ORIENT_SPECIFIC") "integer" else "float"
    lines <- c(lines, block(type, nm, attrs[[nm]]))
  }
  head <- file.path(dir, paste0(name, "+orig.HEAD"))
  writeLines(lines, head)
  con <- file(file.path(dir, paste0(name, "+orig.BRIK")), "wb")
  writeBin(as.numeric(seq_len(prod(dims))), con, size = 4, endian = "little")
  close(con)
  head
}

test_that("AFNI IJK_TO_DICOM_REAL is converted from DICOM to RAS", {
  dir <- tempfile("afni")
  dir.create(dir)
  # Geometry attributes of nibabel's example4d+orig test dataset: an RAI
  # (DICOM) grid, i.e. LPI in RAS terms
  ijk <- c(3, 0, 0, -49.5, 0, 3, 0, -82.312, 0, 0, 3, -52.3511)
  head <- write_afni_pair(dir, "lpi", list(
    ORIENT_SPECIFIC = c(0L, 3L, 4L),
    ORIGIN = c(-49.5, -82.312, -52.3511),
    DELTA = c(3, 3, 3),
    IJK_TO_DICOM = ijk,
    IJK_TO_DICOM_REAL = ijk
  ))
  expected <- rbind(
    c(-3, 0, 0, 49.5),
    c(0, -3, 0, 82.312),
    c(0, 0, 3, -52.3511),
    c(0, 0, 0, 1)
  )
  meta <- read_header(head)
  expect_equal(trans(meta), expected)
  expect_equal(meta@spacing, c(3, 3, 3))
  expect_equal(meta@origin, c(49.5, 82.312, -52.3511))

  vol <- read_vol(head)
  expect_equal(trans(space(vol)), expected)
  # Voxel (0,0,0) sits at the RAS origin, not at the raw DICOM origin
  expect_equal(as.numeric(grid_to_coord(space(vol), matrix(1, 1, 3))),
               c(49.5, 82.312, -52.3511))
})

test_that("AFNI geometry falls back to ORIENT_SPECIFIC/DELTA/ORIGIN", {
  # L2R, P2A, I2S axes: an RAS-ordered grid stored in DICOM coordinates
  hdr <- list(
    ORIENT_SPECIFIC = list(content = c(1L, 2L, 4L)),
    DELTA = list(content = c(-2, -3, 4)),
    ORIGIN = list(content = c(10, 20, -30))
  )
  expect_equal(
    .afni_affine(hdr),
    rbind(c(2, 0, 0, -10), c(0, 3, 0, -20), c(0, 0, 4, -30), c(0, 0, 0, 1))
  )
  # The same fallback on the RAI grid above gives the same affine as its
  # IJK_TO_DICOM matrix
  rai <- list(
    ORIENT_SPECIFIC = list(content = c(0L, 3L, 4L)),
    DELTA = list(content = c(3, 3, 3)),
    ORIGIN = list(content = c(-49.5, -82.312, -52.3511))
  )
  expect_equal(.afni_affine(rai)[1:3, 4], c(49.5, 82.312, -52.3511))
  expect_equal(diag(.afni_affine(rai))[1:3], c(-3, -3, 3))
  expect_error(.afni_affine(list()), "ORIENT_SPECIFIC")
})

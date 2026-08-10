context("NIfTI conformance and the compiled I/O path")

library(neuroim2)

ns <- asNamespace("neuroim2")

# ---------------------------------------------------------------------------
# Helpers: build files by hand so the tests do not depend on the writer they
# are checking, and patch raw header fields the writer would not produce.
# ---------------------------------------------------------------------------

patch_bytes <- function(path, offset, raw_bytes) {
  con <- file(path, open = "r+b")
  on.exit(close(con))
  seek(con, offset, origin = "start", rw = "write")
  writeBin(raw_bytes, con)
  invisible(path)
}

f32 <- function(x) writeBin(as.double(x), raw(), size = 4, endian = "little")
i16 <- function(x) writeBin(as.integer(x), raw(), size = 2, endian = "little")

read_scl <- function(path) {
  con <- file(path, open = "rb"); on.exit(close(con))
  seek(con, 112)
  readBin(con, double(), n = 2, size = 4, endian = "little")
}

read_field_i16 <- function(path, offset) {
  con <- file(path, open = "rb"); on.exit(close(con))
  seek(con, offset)
  readBin(con, integer(), n = 1, size = 2, endian = "little")
}

demo_space3 <- function(d = c(5L, 6L, 4L)) NeuroSpace(d, c(2, 2, 2.5))

# ---------------------------------------------------------------------------
# Datatype coverage
# ---------------------------------------------------------------------------

test_that("every fixed-width NIfTI datatype round-trips", {
  sp <- demo_space3()
  n <- prod(dim(sp))
  cases <- list(
    UBYTE  = c(0, 255, 7),
    BYTE   = c(-128, 127, -3),
    SHORT  = c(-32768, 32767, 12),
    USHORT = c(0, 65535, 1000),
    INT    = c(-2147483648, 2147483647, -7),
    UINT   = c(0, 4294967295, 5),
    LONG   = c(-1e15, 1e15, 3),
    ULONG  = c(0, 1e15, 9),
    FLOAT  = c(-1.5, 2.25, 0),
    DOUBLE = c(-1 / 3, pi, 0)
  )
  for (dt in names(cases)) {
    vals <- rep_len(cases[[dt]], n)
    vol <- NeuroVol(array(vals, dim(sp)), sp)
    f <- tempfile(fileext = ".nii")
    write_vol(vol, f, data_type = dt)
    back <- read_vol(f)
    expect_equal(as.numeric(back@.Data), vals,
                 tolerance = if (dt == "FLOAT") 1e-6 else 0,
                 info = dt)
    expect_equal(header(back)$data_type, dt, info = dt)
    unlink(f)
  }
})

test_that("int32 INT_MIN survives instead of becoming NA", {
  # R's integer NA *is* INT_MIN, so routing an int32 payload through an R
  # integer vector loses that value silently. The compiled reader lands in
  # double and never sees an R integer.
  sp <- demo_space3(c(2L, 2L, 2L))
  vals <- c(-2147483648, 2147483647, 0, 1, -1, 5, -5, 42)
  vol <- NeuroVol(array(vals, dim(sp)), sp)
  f <- tempfile(fileext = ".nii")
  write_vol(vol, f, data_type = "INT")
  back <- as.numeric(read_vol(f)@.Data)
  expect_false(anyNA(back))
  expect_equal(back, vals)
})

test_that("datatypes a real-valued volume cannot hold are refused by name", {
  for (code in c(32L, 128L, 1792L, 2304L, 1L)) {
    expect_error(ns$.getDataStorage(code), "neuroim2 images hold real-valued",
                 info = as.character(code))
  }
  # ... and an unknown code still says so plainly
  expect_error(ns$.getDataStorage(9999L), "Unsupported NIfTI data-type code")
})

# ---------------------------------------------------------------------------
# Scaling conventions
# ---------------------------------------------------------------------------

test_that("all three spellings of \"no scaling\" are honoured on read", {
  sp <- demo_space3(c(3L, 3L, 2L))
  vals <- as.numeric(seq_len(prod(dim(sp))))
  vol <- NeuroVol(array(vals, dim(sp)), sp)

  for (spelling in list(c(0, 0), c(NaN, NaN), c(NA_real_, NA_real_),
                        c(Inf, 0), c(0, 7))) {
    f <- tempfile(fileext = ".nii")
    write_vol(vol, f, data_type = "SHORT")
    patch_bytes(f, 112, c(f32(spelling[1]), f32(spelling[2])))
    # slope == 0 / NaN / NA / Inf all mean identity, and the intercept that
    # accompanies a dead slope is ignored with it.
    expect_equal(as.numeric(read_vol(f)@.Data), vals,
                 info = paste(spelling, collapse = "/"))
    unlink(f)
  }
})

test_that("a real slope and intercept are applied on read", {
  sp <- demo_space3(c(3L, 3L, 2L))
  stored <- as.numeric(seq_len(prod(dim(sp))))
  vol <- NeuroVol(array(stored, dim(sp)), sp)
  f <- tempfile(fileext = ".nii")
  write_vol(vol, f, data_type = "SHORT")
  patch_bytes(f, 112, c(f32(2.5), f32(-4)))
  expect_equal(as.numeric(read_vol(f)@.Data), stored * 2.5 - 4)
})

# ---------------------------------------------------------------------------
# Writing: datatype choice and scaling
# ---------------------------------------------------------------------------

test_that("integer output is scaled to fit rather than truncated", {
  sp <- demo_space3(c(4L, 4L, 4L))
  d <- array(seq(-3.7, 3.7, length.out = prod(dim(sp))), dim(sp))
  vol <- NeuroVol(d, sp)
  f <- tempfile(fileext = ".nii")
  write_vol(vol, f, data_type = "SHORT")

  back <- as.numeric(read_vol(f)@.Data)
  # Truncation toward zero -- what writeBin(as.integer(x)) did -- costs ~0.5;
  # fitting the range to int16 costs ~range/65535.
  expect_lt(max(abs(back - as.numeric(d))), diff(range(d)) / 60000)
  expect_equal(read_field_i16(f, 70), 4L)   # still int16 on disk
  expect_true(read_scl(f)[1] != 1)          # ... via a real slope
})

test_that("data already representable in the target type is written unscaled", {
  sp <- demo_space3(c(4L, 4L, 4L))
  for (dt in c("UBYTE", "SHORT", "INT")) {
    vals <- rep_len(c(0, 1, 0, 1, 2, 3), prod(dim(sp)))
    vol <- NeuroVol(array(vals, dim(sp)), sp)
    f <- tempfile(fileext = ".nii")
    write_vol(vol, f, data_type = dt)
    expect_equal(read_scl(f), c(1, 0), info = dt)
    expect_identical(as.numeric(read_vol(f)@.Data), vals, info = dt)
    unlink(f)
  }
})

test_that("with no data_type the source datatype is kept when it stays exact", {
  sp <- demo_space3()
  vals <- as.numeric(seq_len(prod(dim(sp))))
  src <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(vals, dim(sp)), sp), src, data_type = "SHORT")

  v <- read_vol(src)
  out <- tempfile(fileext = ".nii")
  write_vol(v, out)
  expect_equal(read_field_i16(out, 70), 4L)          # SHORT, not promoted
  expect_identical(as.numeric(read_vol(out)@.Data), vals)

  # ... and promoted to float once the values no longer fit it
  v2 <- NeuroVol(array(vals + 0.5, dim(sp)), sp)
  v2@header <- v@header
  out2 <- tempfile(fileext = ".nii")
  write_vol(v2, out2)
  expect_equal(read_field_i16(out2, 70), 16L)        # FLOAT
  expect_equal(as.numeric(read_vol(out2)@.Data), vals + 0.5, tolerance = 1e-6)
})

test_that("a read-write round trip reproduces the payload byte for byte", {
  sp <- demo_space3()
  vals <- as.numeric(seq_len(prod(dim(sp))))
  src <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(vals, dim(sp)), sp), src, data_type = "SHORT")
  patch_bytes(src, 112, c(f32(0.5), f32(10)))        # a genuinely scaled file

  v <- read_vol(src)
  out <- tempfile(fileext = ".nii")
  write_vol(v, out)

  raw_in <- readBin(src, "raw", n = file.size(src))
  raw_out <- readBin(out, "raw", n = file.size(out))
  # Same datatype, same slope/intercept, same stored integers.
  expect_equal(read_field_i16(out, 70), read_field_i16(src, 70))
  expect_equal(read_scl(out), read_scl(src))
  expect_identical(raw_out[353:length(raw_out)], raw_in[353:length(raw_in)])
})

# ---------------------------------------------------------------------------
# Header preservation
# ---------------------------------------------------------------------------

test_that("acquisition metadata survives a read-write round trip", {
  sp <- NeuroSpace(c(4L, 5L, 3L, 6L), c(2, 2, 3))
  vec <- DenseNeuroVec(array(as.numeric(1:360), dim(sp)), sp)
  src <- tempfile(fileext = ".nii")
  write_vec(vec, src, data_type = "SHORT")

  # Fields that describe the acquisition, not the array. Offsets are NIfTI-1.
  patch_bytes(src, 76 + 4 * 4, f32(2.5))                     # pixdim[4] = TR
  patch_bytes(src, 123, as.raw(10L))                         # xyzt_units mm+sec
  patch_bytes(src, 68, i16(5L))                              # intent_code
  patch_bytes(src, 124, c(f32(3), f32(-3)))                  # cal_max, cal_min
  patch_bytes(src, 132, f32(0.05))                           # slice_duration
  patch_bytes(src, 136, f32(1.5))                            # toffset
  patch_bytes(src, 148, charToRaw("a probe"))                # descrip
  patch_bytes(src, 228, charToRaw("aux.txt"))                # aux_file
  patch_bytes(src, 252, c(i16(1L), i16(4L)))                 # qform/sform code
  patch_bytes(src, 328, charToRaw("zscore"))                 # intent_name

  v <- read_vec(src)
  out <- tempfile(fileext = ".nii")
  write_vec(v, out)

  a <- ns$read_nifti_header(src)
  b <- ns$read_nifti_header(out)
  for (fld in c("pixdim", "xyzt_units", "intent_code", "cal_max", "cal_min",
                "slice_duration", "toffset", "description", "auxfile",
                "qform_code", "sform_code", "intent_name", "dimensions",
                "datatype", "bitpix")) {
    expect_equal(b[[fld]], a[[fld]], info = fld)
  }
  # And the accessor sees them too.
  expect_equal(header(v)$TR, 2.5)
  expect_equal(header(v)$sform$code, 4L)
  expect_equal(header(v)$descrip, "a probe")
})

test_that("an image built in memory writes a neutral header", {
  sp <- demo_space3()
  vol <- NeuroVol(array(1.5, dim(sp)), sp)
  expect_length(vol@header, 0L)
  f <- tempfile(fileext = ".nii")
  write_vol(vol, f)
  h <- ns$read_nifti_header(f)
  expect_equal(h$qform_code, 1L)
  expect_equal(h$sform_code, 1L)
  expect_equal(as.numeric(h$pixdim[5:8]), rep(0, 4))
})

# ---------------------------------------------------------------------------
# Dimension handling
# ---------------------------------------------------------------------------

test_that("singleton axes declared by dim[0] are preserved", {
  # A single-slice volume used to come back 2-D, and a 1-volume run 3-D,
  # because the axis count was inferred by stopping at the first dim == 1.
  single_slice <- NeuroSpace(c(6L, 7L, 1L), c(2, 2, 2))
  f <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(1, dim(single_slice)), single_slice), f)
  expect_equal(dim(read_header(f)), c(6L, 7L, 1L))
  expect_equal(dim(read_vol(f)), c(6L, 7L, 1L))

  one_volume <- NeuroSpace(c(4L, 4L, 3L, 1L), c(2, 2, 2))
  g <- tempfile(fileext = ".nii")
  write_vec(DenseNeuroVec(array(1, dim(one_volume)), one_volume), g)
  expect_equal(dim(read_header(g)), c(4L, 4L, 3L, 1L))
  expect_equal(dim(read_vec(g)), c(4L, 4L, 3L, 1L))
})

test_that("a 5-D file keeps its shape in the header and is not silently folded", {
  # NeuroHyperVec stores [features x trials x in-mask voxels]; the space says
  # what that unpacks to.
  sp5 <- NeuroSpace(c(3L, 4L, 2L, 1L, 3L), c(2, 2, 2))
  mask5 <- LogicalNeuroVol(array(TRUE, dim(sp5)[1:3]),
                           NeuroSpace(dim(sp5)[1:3], c(2, 2, 2)))
  nvox <- prod(dim(sp5)[1:3])
  hv <- NeuroHyperVec(array(as.numeric(seq_len(3L * 1L * nvox)),
                            c(3L, 1L, nvox)), sp5, mask5)
  f <- tempfile(fileext = ".nii")
  write_vec(hv, f)

  expect_equal(dim(read_header(f)), c(3L, 4L, 2L, 1L, 3L))
  h <- read_hyper_vec(f)
  expect_equal(dim(h), c(3L, 4L, 2L, 1L, 3L))
  expect_equal(as.numeric(ns$dense_array_5d(h)),
               as.numeric(ns$dense_array_5d(hv)))

  # read_vec folds the degenerate 4th axis, which is a documented convenience
  # of that entry point rather than a property of the file.
  v <- read_vec(f)
  expect_equal(dim(v), c(3L, 4L, 2L, 3L))

  # A genuinely 5-D file is refused with a pointer to the right reader.
  sp5b <- NeuroSpace(c(3L, 4L, 2L, 2L, 3L), c(2, 2, 2))
  hv2 <- NeuroHyperVec(array(1, c(3L, 2L, nvox)), sp5b, mask5)
  g <- tempfile(fileext = ".nii")
  write_vec(hv2, g)
  expect_error(read_vec(g), "read_hyper_vec")
})

# ---------------------------------------------------------------------------
# Affine provenance
# ---------------------------------------------------------------------------

test_that("qform_code = sform_code = 0 falls back to the spec's METHOD 1", {
  sp <- NeuroSpace(c(8L, 9L, 7L), c(2, 2, 2), origin = c(90, -126, -72))
  f <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(1, dim(sp)), sp), f)
  patch_bytes(f, 252, c(i16(0L), i16(0L)))
  patch_bytes(f, 280, rep(as.raw(0L), 48))          # zero the srow rows

  got <- trans(read_header(f))
  # Voxel size only, image centred on the origin, first axis flipped.
  want <- rbind(c(-2, 0, 0, 7), c(0, 2, 0, -8), c(0, 0, 2, -6), c(0, 0, 0, 1))
  expect_equal(unname(got), want, tolerance = 1e-6)
})

test_that("METHOD 1 is not used when either transform is declared", {
  sp <- NeuroSpace(c(8L, 9L, 7L), c(2, 2, 2), origin = c(90, -126, -72))
  f <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(1, dim(sp)), sp), f)
  base <- trans(read_header(f))

  g <- tempfile(fileext = ".nii"); file.copy(f, g)
  patch_bytes(g, 252, c(i16(1L), i16(0L)))          # qform only
  patch_bytes(g, 280, rep(as.raw(0L), 48))
  expect_equal(unname(trans(read_header(g))), unname(base), tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# NIfTI-2
# ---------------------------------------------------------------------------

test_that("NIfTI-2 round-trips, plain and gzipped", {
  sp <- NeuroSpace(c(5L, 6L, 4L, 3L), c(2, 2, 2.5))
  vals <- as.numeric(seq_len(prod(dim(sp))))
  vec <- DenseNeuroVec(array(vals, dim(sp)), sp)

  for (ext in c(".nii", ".nii.gz")) {
    f <- tempfile(fileext = ext)
    write_vec(vec, f, version = 2)
    h <- read_header(f)
    expect_equal(h@header$version, 2L, info = ext)
    expect_equal(h@header$magic, "n+2", info = ext)
    expect_equal(h@data_offset, 544, info = ext)
    back <- read_vec(f)
    expect_identical(as.numeric(back@.Data), vals, info = ext)
    expect_equal(trans(back), trans(vec), tolerance = 1e-9, info = ext)
    unlink(f)
  }

  # format = strings select it too
  f <- tempfile(fileext = ".nii")
  write_vec(vec, f, "NIFTI2")
  expect_equal(read_header(f)@header$version, 2L)
})

test_that("NIfTI-1 is used unless the image cannot fit it", {
  expect_equal(ns$.nifti_output_version(c(3, 10, 10, 10, 1, 1, 1, 1)), 1L)
  expect_equal(ns$.nifti_output_version(c(3, 40000, 10, 10, 1, 1, 1, 1)), 2L)
  expect_equal(ns$.nifti_output_version(c(3, 10, 10, 10, 1, 1, 1, 1), 2), 2L)
  expect_error(ns$.nifti_output_version(c(3, 40000, 10, 10, 1, 1, 1, 1), 1),
               "cannot store a dimension larger than 32767")
  expect_error(ns$.nifti_output_version(c(3, 10, 10, 10, 1, 1, 1, 1), 3),
               "must be 1 or 2")
})

# ---------------------------------------------------------------------------
# Byte order and compression
# ---------------------------------------------------------------------------

test_that("big-endian files read identically to little-endian ones", {
  sp <- demo_space3()
  vals <- as.numeric(seq_len(prod(dim(sp))))
  vol <- NeuroVol(array(vals, dim(sp)), sp)

  le <- tempfile(fileext = ".nii")
  write_vol(vol, le, data_type = "SHORT")
  hdr <- ns$as_nifti_header(vol, "x", data_type = "SHORT",
                            values = as.numeric(vals))
  hdr$endian <- "big"
  be <- tempfile(fileext = ".nii")
  ns$.write_nifti_file(hdr, be, as.numeric(vals), "SHORT")

  expect_equal(read_header(be)@endian, "big")
  expect_identical(as.numeric(read_vol(be)@.Data), vals)
  expect_identical(as.numeric(read_vol(be)@.Data),
                   as.numeric(read_vol(le)@.Data))
})

test_that("gzip round-trips for both volumes and vectors", {
  sp <- NeuroSpace(c(5L, 6L, 4L, 3L), c(2, 2, 2))
  vals <- as.numeric(seq_len(prod(dim(sp))))
  vec <- DenseNeuroVec(array(vals, dim(sp)), sp)
  f <- tempfile(fileext = ".nii.gz")
  write_vec(vec, f)
  expect_identical(as.numeric(read_vec(f)@.Data), vals)
  expect_identical(as.numeric(read_vol(f, index = 2)@.Data),
                   vals[prod(dim(sp)[1:3]) + seq_len(prod(dim(sp)[1:3]))])
})

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

test_that("a truncated file says what was expected and what arrived", {
  sp <- demo_space3()
  f <- tempfile(fileext = ".nii")
  write_vol(NeuroVol(array(1, dim(sp)), sp), f, data_type = "SHORT")
  b <- readBin(f, "raw", n = file.size(f))
  short <- tempfile(fileext = ".nii")
  writeBin(b[seq_len(length(b) - 40L)], short)

  expect_error(read_vol(short), "truncated")
  expect_error(read_vol(short), "supplies only")
})

test_that("a file that is not an image is rejected with the reason", {
  f <- tempfile(fileext = ".nii")
  writeBin(as.raw(rep(0L, 400)), f)
  expect_error(read_header(f), "does not start with a NIfTI")

  g <- tempfile(fileext = ".nii")
  writeBin(as.raw(rep(0L, 3)), g)
  expect_error(read_header(g), "3 bytes long")
})

# ---------------------------------------------------------------------------
# The compiled path itself
# ---------------------------------------------------------------------------

test_that("in-place dense constructors are identical to new()", {
  sp4 <- NeuroSpace(c(3L, 4L, 2L, 5L), c(2, 2, 2))
  d4 <- as.numeric(seq_len(prod(dim(sp4))))
  a <- new("DenseNeuroVec", .Data = array(d4, dim(sp4)), space = sp4,
           label = "L", volume_labels = rep("", 5L), header = list(x = 1))
  b <- ns$.dense_neurovec_inplace(d4, sp4, "L", rep("", 5L), list(x = 1))
  expect_identical(a, b)
  expect_true(isTRUE(validObject(b)))

  sp3 <- demo_space3()
  d3 <- as.numeric(seq_len(prod(dim(sp3))))
  c1 <- new("DenseNeuroVol", .Data = array(d3, dim(sp3)), space = sp3,
            header = list(y = 2))
  c2 <- ns$.dense_neurovol_inplace(d3, sp3, list(y = 2))
  expect_identical(c1, c2)
  expect_true(isTRUE(validObject(c2)))

  expect_error(ns$.dense_neurovol_inplace(seq_len(4), sp3), "plain double")
  expect_error(ns$.dense_neurovol_inplace(as.numeric(1:3), sp3), "does not match")
})

test_that("read_mapped_vols returns voxels-by-volumes and applies scaling", {
  sp <- NeuroSpace(c(4L, 3L, 2L, 4L), c(2, 2, 2))
  nvox <- prod(dim(sp)[1:3])
  vals <- as.numeric(seq_len(prod(dim(sp))))
  f <- tempfile(fileext = ".nii")
  write_vec(DenseNeuroVec(array(vals, dim(sp)), sp), f, data_type = "SHORT")
  patch_bytes(f, 112, c(f32(3), f32(1)))

  meta <- read_header(f)
  m <- ns$read_mapped_vols(meta, c(3L, 1L))
  expect_equal(dim(m), c(nvox, 2L))
  expect_equal(m[, 1], vals[(2 * nvox + 1):(3 * nvox)] * 3 + 1)
  expect_equal(m[, 2], vals[1:nvox] * 3 + 1)

  expect_error(ns$read_mapped_vols(meta, 0L), "must be in range")
  expect_error(ns$read_mapped_vols(meta, 99L), "must be in range")
  expect_equal(dim(ns$read_mapped_vols(meta, integer(0))), c(nvox, 0L))
})

test_that("the compiled reader rejects hostile arguments", {
  f <- tempfile(fileext = ".nii")
  sp <- demo_space3()
  write_vol(NeuroVol(array(1, dim(sp)), sp), f, data_type = "SHORT")
  expect_error(ns$nifti_read_data_cpp(f, 352, -1, 4L, FALSE, FALSE),
               "non-negative")
  expect_error(ns$nifti_read_data_cpp(f, -1, 10, 4L, FALSE, FALSE),
               "non-negative")
  expect_error(ns$nifti_read_data_cpp(f, 352, 10, 12345L, FALSE, FALSE),
               "unhandled datatype")
  expect_error(ns$nifti_read_data_cpp(tempfile(), 352, 10, 4L, FALSE, FALSE),
               "cannot open")
  expect_length(ns$nifti_read_data_cpp(f, 352, 0, 4L, FALSE, FALSE), 0L)
})

test_that("write and read agree on non-finite payloads", {
  sp <- demo_space3(c(3L, 3L, 2L))
  vals <- c(NaN, Inf, -Inf, rep(1.5, prod(dim(sp)) - 3L))
  vol <- NeuroVol(array(vals, dim(sp)), sp)
  for (dt in c("FLOAT", "DOUBLE")) {
    f <- tempfile(fileext = ".nii")
    write_vol(vol, f, data_type = dt)
    back <- as.numeric(read_vol(f)@.Data)
    expect_true(is.nan(back[1]), info = dt)
    expect_identical(back[2:3], c(Inf, -Inf), info = dt)
    unlink(f)
  }
  # An integer type has no representation for a NaN or an infinity. They are
  # stored as 0, which comes back as the scaling intercept -- finite, and
  # indistinguishable from a real 0. What must not happen is a garbage value or
  # a saturated extreme leaking into the data.
  spread <- c(NaN, Inf, -Inf, seq(-4, 4, length.out = prod(dim(sp)) - 3L))
  volx <- NeuroVol(array(spread, dim(sp)), sp)
  f <- tempfile(fileext = ".nii")
  write_vol(volx, f, data_type = "SHORT")
  back <- as.numeric(read_vol(f)@.Data)
  expect_true(all(is.finite(back)))
  expect_equal(back[-(1:3)], spread[-(1:3)], tolerance = 1e-3)
  expect_true(all(abs(back[1:3]) <= max(abs(spread[-(1:3)]))))
})

# ---------------------------------------------------------------------------
# 4-D Gaussian blur
# ---------------------------------------------------------------------------

test_that("gaussian_blur on a NeuroVec matches per-volume blurring", {
  sp <- NeuroSpace(c(8L, 9L, 6L, 4L), c(2, 2, 3))
  set.seed(4)
  vec <- DenseNeuroVec(array(rnorm(prod(dim(sp))), dim(sp)), sp)
  mk <- LogicalNeuroVol(array(runif(prod(dim(sp)[1:3])) < 0.4, dim(sp)[1:3]),
                        NeuroSpace(dim(sp)[1:3], c(2, 2, 3)))

  for (w in 1:2) for (nm in c(TRUE, FALSE)) {
    got <- gaussian_blur(vec, sigma = 1.5, window = w, normalize = nm)
    want <- vapply(seq_len(dim(sp)[4]), function(i)
      as.numeric(gaussian_blur(vec[[i]], sigma = 1.5, window = w,
                               normalize = nm)),
      numeric(prod(dim(sp)[1:3])))
    expect_equal(as.numeric(got@.Data), as.numeric(want),
                 info = paste(w, nm))

    got_m <- gaussian_blur(vec, mk, sigma = 1.5, window = w, normalize = nm)
    want_m <- vapply(seq_len(dim(sp)[4]), function(i)
      as.numeric(gaussian_blur(vec[[i]], mk, sigma = 1.5, window = w,
                               normalize = nm)),
      numeric(prod(dim(sp)[1:3])))
    expect_equal(as.numeric(got_m@.Data), as.numeric(want_m),
                 info = paste("masked", w, nm))
  }

  expect_s4_class(gaussian_blur(vec, sigma = 1.5), "DenseNeuroVec")
  expect_equal(dim(gaussian_blur(vec, sigma = 1.5)), dim(sp))
})

test_that("the full-mask shortcut agrees with an explicit index set", {
  set.seed(5)
  for (d in list(c(9L, 11L, 7L), c(5L, 4L, 12L), c(1L, 7L, 7L))) {
    arr <- array(rnorm(prod(d)), d)
    full <- seq_len(prod(d))
    for (w in 1:3) for (sg in c(0.5, 2.5)) for (nm in c(TRUE, FALSE)) {
      short <- ns$gaussian_blur_sep_cpp(arr, integer(0), w, sg, c(1, 1, 1), nm, TRUE)
      long  <- ns$gaussian_blur_sep_cpp(arr, full, w, sg, c(1, 1, 1), nm, FALSE)
      dense <- ns$gaussian_blur_cpp(arr, full, w, sg, c(1, 1, 1), nm)
      expect_identical(as.numeric(short), as.numeric(long),
                       info = paste(d, w, sg, nm))
      expect_equal(as.numeric(short), as.numeric(dense), tolerance = 1e-12,
                   info = paste(d, w, sg, nm))
    }
  }
})

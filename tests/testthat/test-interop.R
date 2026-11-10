# as_neuro() / as_nifti(): moving images between neuroim2 and RNifti (and,
# through RNifti, oro.nifti) without writing a file.

interop_file <- function(name = "mni_downsampled.nii.gz") {
  system.file("extdata", name, package = "neuroim2")
}

plain <- function(m) {
  m <- unclass(m)
  attributes(m) <- list(dim = dim(m))
  m
}

test_that("an RNifti image becomes the volume read_vol() gives", {
  f <- interop_file()
  img <- RNifti::readNifti(f)
  vol <- as_neuro(img)
  ref <- read_vol(f)

  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(dim(vol), dim(ref))
  expect_equal(trans(vol), trans(ref))
  expect_equal(spacing(vol), spacing(ref))
  expect_equal(as.array(vol), as.array(ref))
  expect_equal(trans(vol), plain(RNifti::xform(img)), tolerance = 1e-6)
})

test_that("a 4-D RNifti image becomes the vector read_vec() gives", {
  f <- interop_file("global_mask_v4.nii")
  vec <- as_neuro(RNifti::readNifti(f))
  ref <- read_vec(f)

  expect_s4_class(vec, "DenseNeuroVec")
  expect_equal(dim(vec), dim(ref))
  expect_equal(trans(vec), trans(ref))
  expect_equal(as.array(vec), as.array(ref))
})

test_that("an internal (pointer-only) RNifti image converts", {
  f <- interop_file()
  img <- RNifti::readNifti(f, internal = TRUE)
  expect_s3_class(img, "internalImage")
  vol <- as_neuro(img)
  expect_equal(as.array(vol), as.array(read_vol(f)))
  expect_equal(trans(vol), trans(read_vol(f)))
})

test_that("a trailing singleton dimension is dropped", {
  img <- RNifti::asNifti(array(1:24, c(2, 3, 4, 1)))
  vol <- as_neuro(img)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(dim(vol), c(2L, 3L, 4L))
})

test_that("as_nifti() carries values and geometry", {
  ref <- read_vol(interop_file())
  img <- as_nifti(ref)

  expect_s3_class(img, "niftiImage")
  expect_equal(dim(img), dim(ref))
  expect_equal(plain(img), as.array(ref))
  expect_equal(plain(RNifti::xform(img)), trans(ref), tolerance = 1e-6)
  expect_equal(RNifti::pixdim(img), spacing(ref), tolerance = 1e-5)
  hdr <- RNifti::niftiHeader(img)
  expect_equal(hdr$sform_code, 1L)
  expect_equal(hdr$qform_code, 1L)
  expect_equal(RNifti::niftiHeader(as_nifti(ref, xform_code = 4))$sform_code, 4L)
  expect_error(as_nifti(ref, xform_code = 9), "xform_code")
})

test_that("RNifti::asNifti() dispatches on neuroim2 objects", {
  ref <- read_vol(interop_file())
  a <- RNifti::asNifti(ref)
  b <- as_nifti(ref)
  expect_s3_class(a, "niftiImage")
  expect_equal(plain(a), plain(b))
  expect_equal(plain(RNifti::xform(a)), plain(RNifti::xform(b)))
})

test_that("neuroim2 -> RNifti -> neuroim2 is the identity", {
  vol <- read_vol(interop_file())
  back <- as_neuro(as_nifti(vol))
  expect_equal(trans(back), trans(vol))
  expect_equal(as.array(back), as.array(vol))

  vec <- read_vec(interop_file("global_mask_v4.nii"))
  back4 <- as_neuro(as_nifti(vec))
  expect_s4_class(back4, "DenseNeuroVec")
  expect_equal(dim(back4), dim(vec))
  expect_equal(trans(back4), trans(vec))
  expect_equal(as.array(back4), as.array(vec))
})

test_that("an oblique, sheared affine survives through the sform", {
  th <- pi / 7
  rot <- rbind(c(cos(th), -sin(th), 0), c(sin(th), cos(th), 0), c(0, 0, 1))
  lin <- rot %*% diag(c(2, 2.5, 3))
  lin[1, 3] <- 0.4   # shear: the qform cannot hold this, the sform can
  aff <- rbind(cbind(lin, c(-40, 12, 7)), c(0, 0, 0, 1))
  sp <- NeuroSpace(c(5L, 6L, 7L), trans = aff)
  vol <- NeuroVol(array(rnorm(5 * 6 * 7), c(5, 6, 7)), sp)

  img <- as_nifti(vol)
  expect_equal(plain(RNifti::xform(img)), trans(sp), tolerance = 1e-6)
  expect_equal(trans(as_neuro(img)), trans(sp), tolerance = 1e-6)
})

test_that("a file written by RNifti from as_nifti() reads back the same", {
  vol <- read_vol(interop_file())
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))
  RNifti::writeNifti(as_nifti(vol), tmp)
  back <- read_vol(tmp)
  expect_equal(trans(back), trans(vol))
  expect_equal(as.array(back), as.array(vol), tolerance = 1e-6)
})

test_that("sparse and logical volumes convert through as.array()", {
  vol <- read_vol(interop_file())
  mask <- as.mask(vol > 0)
  mimg <- as_nifti(mask)
  expect_type(plain(mimg), "integer")
  expect_equal(sum(mimg), sum(mask))

  sv <- as.sparse(vol, which(mask))
  simg <- as_nifti(sv)
  expect_equal(plain(simg), as.array(vol) * as.array(mask))
  expect_equal(plain(RNifti::xform(simg)), trans(vol), tolerance = 1e-6)
})

test_that("as_neuro() passes neuroim2 objects and file names through", {
  f <- interop_file()
  vol <- read_vol(f)
  expect_identical(as_neuro(vol), vol)
  expect_equal(as.array(as_neuro(f)), as.array(vol))
  expect_s4_class(as_neuro(interop_file("global_mask_v4.nii")), "NeuroVec")
})

test_that("unsupported inputs say what is supported", {
  expect_error(as_neuro(array(0, c(2, 2, 2))), "NeuroSpace")
  expect_error(as_neuro(RNifti::asNifti(matrix(0, 3, 3))), "3-D and 4-D")
  expect_error(as_nifti(array(0, c(2, 2, 2))), "asNifti")
})

test_that("oro.nifti objects convert in both directions", {
  skip_if_not_installed("oro.nifti")
  f <- interop_file()
  ref <- read_vol(f)
  o <- oro.nifti::readNIfTI(f, reorient = FALSE)

  vol <- as_neuro(o)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(trans(vol), trans(ref))
  expect_equal(as.array(vol), as.array(ref))

  o2 <- oro.nifti::nii2oro(as_nifti(ref))
  expect_s4_class(o2, "nifti")
  expect_equal(trans(as_neuro(o2)), trans(ref), tolerance = 1e-6)
  expect_equal(as.array(as_neuro(o2)), as.array(ref))
})

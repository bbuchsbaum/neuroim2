library(testthat)
library(neuroim2)

# Helper: build AxisSet objects via new() since constructors are not exported
make_ax1 <- function(ax = LEFT_RIGHT) new("AxisSet1D", ndim = 1L, i = ax)
make_ax2 <- function(i = LEFT_RIGHT, j = POST_ANT) new("AxisSet2D", ndim = 2L, i = i, j = j)
make_ax3 <- function(i = LEFT_RIGHT, j = POST_ANT, k = INF_SUP) new("AxisSet3D", ndim = 3L, i = i, j = j, k = k)

# ---- AxisSet constructors via new() ----

test_that("AxisSet1D created via new() has correct slots", {
  ax <- make_ax1(LEFT_RIGHT)
  expect_s4_class(ax, "AxisSet1D")
  expect_equal(ax@ndim, 1L)
  expect_identical(ax@i, LEFT_RIGHT)
})

test_that("AxisSet2D created via new() has correct slots", {
  ax <- make_ax2(LEFT_RIGHT, POST_ANT)
  expect_s4_class(ax, "AxisSet2D")
  expect_equal(ax@ndim, 2L)
  expect_identical(ax@i, LEFT_RIGHT)
  expect_identical(ax@j, POST_ANT)
})

test_that("AxisSet3D created via new() has correct slots", {
  ax <- make_ax3(LEFT_RIGHT, POST_ANT, INF_SUP)
  expect_s4_class(ax, "AxisSet3D")
  expect_equal(ax@ndim, 3L)
  expect_identical(ax@i, LEFT_RIGHT)
  expect_identical(ax@j, POST_ANT)
  expect_identical(ax@k, INF_SUP)
})

# ---- ndim ----

test_that("ndim returns correct dimension for AxisSet objects", {
  expect_equal(ndim(make_ax1(TIME)),             1L)
  expect_equal(ndim(make_ax2(LEFT_RIGHT, POST_ANT)), 2L)
  expect_equal(ndim(make_ax3()),                 3L)
})

# ---- perm_mat ----

test_that("perm_mat returns correct matrix for AxisSet2D", {
  ax <- make_ax2(LEFT_RIGHT, POST_ANT)
  pm <- perm_mat(ax)
  expect_true(is.matrix(pm))
  expect_equal(dim(pm), c(3L, 2L))
  expect_equal(pm[, 1], LEFT_RIGHT@direction)
  expect_equal(pm[, 2], POST_ANT@direction)
})

test_that("perm_mat returns 3x3 matrix for AxisSet3D", {
  ax <- make_ax3()
  pm <- perm_mat(ax)
  expect_true(is.matrix(pm))
  expect_equal(dim(pm), c(3L, 3L))
})

test_that("perm_mat works via NeuroSpace", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  pm <- perm_mat(sp)
  expect_true(is.matrix(pm))
})

# ---- drop_dim ----

test_that("drop_dim on AxisSet2D without dimnum returns AxisSet1D", {
  ax  <- make_ax2(LEFT_RIGHT, POST_ANT)
  out <- drop_dim(ax)
  expect_s4_class(out, "AxisSet1D")
  expect_identical(out@i, LEFT_RIGHT)
})

test_that("drop_dim on AxisSet2D with dimnum=1 keeps second axis", {
  ax  <- make_ax2(LEFT_RIGHT, POST_ANT)
  out <- drop_dim(ax, 1)
  expect_s4_class(out, "AxisSet1D")
  expect_identical(out@i, POST_ANT)
})

test_that("drop_dim on AxisSet2D with dimnum=2 keeps first axis", {
  ax  <- make_ax2(LEFT_RIGHT, POST_ANT)
  out <- drop_dim(ax, 2)
  expect_s4_class(out, "AxisSet1D")
  expect_identical(out@i, LEFT_RIGHT)
})

test_that("drop_dim on AxisSet3D without dimnum drops third axis", {
  ax  <- make_ax3(LEFT_RIGHT, POST_ANT, INF_SUP)
  out <- drop_dim(ax)
  expect_s4_class(out, "AxisSet2D")
  expect_identical(out@i, LEFT_RIGHT)
  expect_identical(out@j, POST_ANT)
})

test_that("drop_dim on AxisSet3D with dimnum=1 removes first axis", {
  ax  <- make_ax3(LEFT_RIGHT, POST_ANT, INF_SUP)
  out <- drop_dim(ax, 1)
  expect_s4_class(out, "AxisSet2D")
  expect_identical(out@i, POST_ANT)
  expect_identical(out@j, INF_SUP)
})

test_that("drop_dim on AxisSet3D with dimnum=2 removes second axis", {
  ax  <- make_ax3(LEFT_RIGHT, POST_ANT, INF_SUP)
  out <- drop_dim(ax, 2)
  expect_s4_class(out, "AxisSet2D")
  expect_identical(out@i, LEFT_RIGHT)
  expect_identical(out@j, INF_SUP)
})

test_that("drop_dim on AxisSet3D with dimnum=3 removes third axis", {
  ax  <- make_ax3(LEFT_RIGHT, POST_ANT, INF_SUP)
  out <- drop_dim(ax, 3)
  expect_s4_class(out, "AxisSet2D")
  expect_identical(out@i, LEFT_RIGHT)
  expect_identical(out@j, POST_ANT)
})

# ---- findAnatomy3D / matchAnatomy3D ----

test_that("findAnatomy3D with default LPI returns AxisSet3D", {
  ax <- findAnatomy3D()
  expect_s4_class(ax, "AxisSet3D")
})

test_that("findAnatomy3D with RAS codes returns correct axes", {
  ax <- findAnatomy3D("R", "A", "S")
  expect_s4_class(ax, "AxisSet3D")
  expect_identical(ax@i, RIGHT_LEFT)
  expect_identical(ax@j, ANT_POST)
  expect_identical(ax@k, SUP_INF)
})

test_that("findAnatomy3D with full-word codes works", {
  ax <- findAnatomy3D("LEFT", "ANTERIOR", "INFERIOR")
  expect_s4_class(ax, "AxisSet3D")
})

test_that("findAnatomy3D with invalid axis errors", {
  expect_error(findAnatomy3D("X", "A", "S"))
})

# ---- affine_to_orientation ----

test_that("affine_to_orientation on identity returns RAS orientation", {
  aff  <- diag(4)
  ornt <- affine_to_orientation(aff)
  expect_true(is.matrix(ornt))
  expect_equal(nrow(ornt), 3L)
  expect_equal(ncol(ornt), 2L)
  expect_equal(unname(ornt[, 1]), c(1, 2, 3))
  expect_equal(unname(ornt[, 2]), c(1, 1, 1))
})

test_that("affine_to_orientation detects axis flips on negative diagonal", {
  aff  <- diag(c(-2, 3, 4, 1))
  ornt <- affine_to_orientation(aff)
  # First column negative scaling => flip = -1
  expect_equal(unname(ornt[1, "flip"]), -1)
  expect_equal(unname(ornt[2, "flip"]),  1)
  expect_equal(unname(ornt[3, "flip"]),  1)
})

test_that("affine_to_orientation rejects non-matrix input", {
  expect_error(affine_to_orientation(list(1, 2, 3)))
})

test_that("affine_to_orientation rejects too-small matrix", {
  expect_error(affine_to_orientation(matrix(1, 1, 1)))
})

# ---- orientation_to_axcodes ----

test_that("orientation_to_axcodes on identity returns R A S", {
  ornt  <- affine_to_orientation(diag(4))
  codes <- orientation_to_axcodes(ornt)
  expect_equal(codes, c("R", "A", "S"))
})

test_that("orientation_to_axcodes detects L P I for negative identity", {
  aff   <- diag(c(-1, -1, -1, 1))
  ornt  <- affine_to_orientation(aff)
  codes <- orientation_to_axcodes(ornt)
  expect_equal(codes, c("L", "P", "I"))
})

# ---- axcodes_to_orientation ----

test_that("axcodes_to_orientation roundtrips with orientation_to_axcodes for RAS", {
  codes  <- c("R", "A", "S")
  ornt   <- axcodes_to_orientation(codes)
  expect_true(is.matrix(ornt))
  expect_equal(nrow(ornt), 3L)
  codes2 <- orientation_to_axcodes(ornt)
  expect_equal(codes, codes2)
})

test_that("axcodes_to_orientation roundtrips for LPI", {
  codes <- c("L", "P", "I")
  ornt  <- axcodes_to_orientation(codes)
  expect_equal(orientation_to_axcodes(ornt), codes)
})

test_that("axcodes_to_orientation sets correct axis indices and flip signs for RAS", {
  ornt <- axcodes_to_orientation(c("R", "A", "S"))
  # R = positive end of L-R axis (axis 1, flip +1)
  expect_equal(unname(ornt[1, "axis"]), 1)
  expect_equal(unname(ornt[1, "flip"]), 1)
  # A = positive end of P-A axis (axis 2, flip +1)
  expect_equal(unname(ornt[2, "axis"]), 2)
  expect_equal(unname(ornt[2, "flip"]), 1)
  # S = positive end of I-S axis (axis 3, flip +1)
  expect_equal(unname(ornt[3, "axis"]), 3)
  expect_equal(unname(ornt[3, "flip"]), 1)
})

test_that("axcodes_to_orientation errors on invalid code", {
  expect_error(axcodes_to_orientation(c("X", "A", "S")))
})

test_that("axcodes_to_orientation errors on duplicated axis codes", {
  expect_error(axcodes_to_orientation(c("R", "L", "S")))
})

# ---- affine_to_axcodes ----

test_that("affine_to_axcodes on identity affine gives R A S", {
  codes <- affine_to_axcodes(diag(4))
  expect_equal(codes, c("R", "A", "S"))
})

test_that("affine_to_axcodes with negative diagonal gives L P I", {
  codes <- affine_to_axcodes(diag(c(-1, -1, -1, 1)))
  expect_equal(codes, c("L", "P", "I"))
})

# ---- orientation_transform ----

test_that("orientation_transform identity when start == end", {
  ornt  <- axcodes_to_orientation(c("R", "A", "S"))
  xform <- orientation_transform(ornt, ornt)
  expect_equal(unname(xform[, "axis"]), c(1, 2, 3))
  expect_equal(unname(xform[, "flip"]), c(1, 1, 1))
})

test_that("orientation_transform detects flip when sign changes", {
  start <- axcodes_to_orientation(c("R", "A", "S"))
  end   <- axcodes_to_orientation(c("L", "A", "S"))
  xform <- orientation_transform(start, end)
  expect_equal(unname(xform[1, "flip"]), -1)
  expect_equal(unname(xform[2, "flip"]),  1)
  expect_equal(unname(xform[3, "flip"]),  1)
})

test_that("orientation_transform detects axis permutation", {
  start <- axcodes_to_orientation(c("R", "A", "S"))
  end   <- axcodes_to_orientation(c("A", "R", "S"))
  xform <- orientation_transform(start, end)
  expect_true(is.matrix(xform))
  expect_equal(nrow(xform), 3L)
})

test_that("orientation_transform errors when shapes differ", {
  s2 <- axcodes_to_orientation(c("R", "A"))
  s3 <- axcodes_to_orientation(c("R", "A", "S"))
  expect_error(orientation_transform(s2, s3))
})

# ---- apply_orientation ----

test_that("apply_orientation identity transform leaves array unchanged", {
  arr   <- array(seq_len(24), c(2, 3, 4))
  ornt  <- axcodes_to_orientation(c("R", "A", "S"))
  xform <- orientation_transform(ornt, ornt)
  out   <- apply_orientation(arr, xform)
  expect_equal(dim(out), dim(arr))
  expect_equal(as.numeric(out), as.numeric(arr))
})

test_that("apply_orientation flip on axis 1 reverses first dimension", {
  arr  <- array(1:8, c(2, 2, 2))
  ornt <- matrix(c(1, -1, 2, 1, 3, 1), nrow = 3, ncol = 2, byrow = TRUE,
                 dimnames = list(NULL, c("axis", "flip")))
  out  <- apply_orientation(arr, ornt)
  expect_equal(dim(out), c(2L, 2L, 2L))
  expect_equal(out[1, , ], arr[2, , ])
  expect_equal(out[2, , ], arr[1, , ])
})

test_that("apply_orientation axis permutation changes dims", {
  arr   <- array(seq_len(24), c(2, 3, 4))
  start <- axcodes_to_orientation(c("R", "A", "S"))
  end   <- axcodes_to_orientation(c("A", "R", "S"))
  xform <- orientation_transform(start, end)
  out   <- apply_orientation(arr, xform)
  expect_equal(length(out), length(arr))
})

test_that("apply_orientation errors when array has fewer dims than ornt", {
  arr  <- array(1:6, c(2, 3))
  ornt <- axcodes_to_orientation(c("R", "A", "S"))
  expect_error(apply_orientation(arr, ornt))
})

# ---- orientation_inverse_affine ----

test_that("orientation_inverse_affine returns square homogeneous matrix", {
  ornt  <- axcodes_to_orientation(c("R", "A", "S"))
  shape <- c(10, 12, 8)
  aff   <- orientation_inverse_affine(ornt, shape)
  expect_true(is.matrix(aff))
  expect_equal(dim(aff), c(4L, 4L))
})

test_that("orientation_inverse_affine for no-flip no-permutation has identity rotation", {
  ornt  <- axcodes_to_orientation(c("R", "A", "S"))
  shape <- c(10, 10, 10)
  aff   <- orientation_inverse_affine(ornt, shape)
  expect_equal(aff[1:3, 1:3], diag(3), tolerance = 1e-10, ignore_attr = TRUE)
})

# ---- show methods (smoke tests) ----

test_that("show method for NamedAxis does not error", {
  expect_output(show(LEFT_RIGHT))
})

test_that("show method for AxisSet1D does not error", {
  expect_output(show(TimeAxis))
})

test_that("show method for AxisSet2D does not error", {
  ax <- make_ax2(LEFT_RIGHT, POST_ANT)
  expect_output(show(ax))
})

test_that("show method for AxisSet3D does not error", {
  ax <- make_ax3()
  expect_output(show(ax))
})

# ---- OrientationList objects ----

test_that("OrientationList2D is a named list of AxisSet2D", {
  expect_true(is.list(OrientationList2D))
  expect_true(all(vapply(OrientationList2D, is, logical(1), "AxisSet2D")))
})

test_that("OrientationList3D is a named list of AxisSet3D", {
  expect_true(is.list(OrientationList3D))
  expect_true(all(vapply(OrientationList3D, is, logical(1), "AxisSet3D")))
})

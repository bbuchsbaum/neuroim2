## Tests for Arith, Compare, Logic, and ! operations on NeuroVol types

sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))

make_dense <- function(vals = rnorm(125)) {
  DenseNeuroVol(array(vals, c(5, 5, 5)), sp)
}

make_logical <- function() {
  LogicalNeuroVol(array(sample(c(TRUE, FALSE), 125, replace = TRUE), c(5, 5, 5)), sp)
}

make_sparse <- function() {
  SparseNeuroVol(c(1.0, 2.0, 3.0), sp, indices = c(1L, 10L, 50L))
}

make_clustered <- function() {
  mask <- LogicalNeuroVol(array(c(rep(TRUE, 50), rep(FALSE, 75)), c(5, 5, 5)), sp)
  clusters <- rep(c(1L, 2L), length.out = 50)
  ClusteredNeuroVol(mask, clusters)
}

# ---------------------------------------------------------------------------
# as.dense identity
# ---------------------------------------------------------------------------
test_that("as.dense is identity for DenseNeuroVol", {
  vol <- make_dense()
  expect_identical(as.dense(vol), vol)
})

test_that("as.dense is identity for LogicalNeuroVol", {
  lvol <- make_logical()
  expect_identical(as.dense(lvol), lvol)
})

# ---------------------------------------------------------------------------
# Scalar Arith — DenseNeuroVol
# ---------------------------------------------------------------------------
test_that("DenseNeuroVol + numeric returns DenseNeuroVol with correct space", {
  vol <- make_dense()
  r <- vol + 5
  expect_s4_class(r, "DenseNeuroVol")
  expect_equal(r@.Data, vol@.Data + 5)
  expect_identical(space(r), space(vol))
})

test_that("numeric * DenseNeuroVol returns DenseNeuroVol", {
  vol <- make_dense()
  r <- 2 * vol
  expect_s4_class(r, "DenseNeuroVol")
  expect_equal(r@.Data, 2 * vol@.Data)
})

# ---------------------------------------------------------------------------
# Scalar Arith — SparseNeuroVol
# ---------------------------------------------------------------------------
test_that("SparseNeuroVol + numeric returns DenseNeuroVol", {
  svol <- make_sparse()
  r <- svol + 1
  expect_s4_class(r, "DenseNeuroVol")
  expect_equal(r[1, 1, 1], 2)   # index 1 had value 1.0
})

test_that("numeric * SparseNeuroVol returns DenseNeuroVol", {
  svol <- make_sparse()
  r <- 3 * svol
  expect_s4_class(r, "DenseNeuroVol")
})

# ---------------------------------------------------------------------------
# ClusteredNeuroVol Arith — warns and returns DenseNeuroVol
# ---------------------------------------------------------------------------
test_that("ClusteredNeuroVol * numeric warns and returns DenseNeuroVol", {
  cvol <- make_clustered()
  expect_warning(r <- cvol * 2, "cluster structure")
  expect_s4_class(r, "DenseNeuroVol")
})

test_that("ClusteredNeuroVol + ClusteredNeuroVol warns", {
  cvol <- make_clustered()
  expect_warning(r <- cvol + cvol, "cluster structure")
  expect_s4_class(r, "DenseNeuroVol")
})

test_that("ClusteredNeuroVol + DenseNeuroVol warns", {
  cvol <- make_clustered()
  dvol <- make_dense()
  expect_warning(r <- cvol + dvol, "cluster structure")
  expect_s4_class(r, "DenseNeuroVol")
})

test_that("DenseNeuroVol + ClusteredNeuroVol warns", {
  cvol <- make_clustered()
  dvol <- make_dense()
  expect_warning(r <- dvol + cvol, "cluster structure")
  expect_s4_class(r, "DenseNeuroVol")
})

# ---------------------------------------------------------------------------
# Compare — returns LogicalNeuroVol
# ---------------------------------------------------------------------------
test_that("DenseNeuroVol > numeric returns LogicalNeuroVol", {
  vol <- make_dense()
  r <- vol > 0
  expect_s4_class(r, "LogicalNeuroVol")
  expect_identical(space(r), space(vol))
})

test_that("numeric < DenseNeuroVol returns LogicalNeuroVol", {
  vol <- make_dense()
  r <- 0 < vol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("DenseNeuroVol == DenseNeuroVol returns LogicalNeuroVol", {
  vol <- make_dense()
  r <- vol == vol
  expect_s4_class(r, "LogicalNeuroVol")
  expect_true(all(r@.Data))
})

test_that("SparseNeuroVol > numeric returns LogicalNeuroVol", {
  svol <- make_sparse()
  r <- svol > 1
  expect_s4_class(r, "LogicalNeuroVol")
  expect_identical(space(r), space(svol))
})

test_that("numeric >= SparseNeuroVol returns LogicalNeuroVol", {
  svol <- make_sparse()
  r <- 2 >= svol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("ClusteredNeuroVol > numeric returns LogicalNeuroVol", {
  cvol <- make_clustered()
  r <- cvol > 1
  expect_s4_class(r, "LogicalNeuroVol")
})

# ---------------------------------------------------------------------------
# Logic — & and |
# ---------------------------------------------------------------------------
test_that("DenseNeuroVol & DenseNeuroVol returns LogicalNeuroVol", {
  v1 <- make_dense(sample(0:1, 125, replace = TRUE))
  v2 <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- v1 & v2
  expect_s4_class(r, "LogicalNeuroVol")
  expect_identical(space(r), space(v1))
})

test_that("DenseNeuroVol | DenseNeuroVol returns LogicalNeuroVol", {
  v1 <- make_dense(sample(0:1, 125, replace = TRUE))
  v2 <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- v1 | v2
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("LogicalNeuroVol & LogicalNeuroVol returns LogicalNeuroVol", {
  l1 <- make_logical()
  l2 <- make_logical()
  r <- l1 & l2
  expect_s4_class(r, "LogicalNeuroVol")
  expect_equal(r@.Data, l1@.Data & l2@.Data)
})

test_that("SparseNeuroVol & DenseNeuroVol returns LogicalNeuroVol", {
  svol <- make_sparse()
  dvol <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- svol & dvol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("ClusteredNeuroVol & DenseNeuroVol returns LogicalNeuroVol", {
  cvol <- make_clustered()
  dvol <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- cvol & dvol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("ClusteredNeuroVol | SparseNeuroVol returns LogicalNeuroVol", {
  cvol <- make_clustered()
  svol <- make_sparse()
  r <- cvol | svol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("NeuroVol & logical scalar returns LogicalNeuroVol", {
  vol <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- vol & TRUE
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("logical scalar | NeuroVol returns LogicalNeuroVol", {
  vol <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- FALSE | vol
  expect_s4_class(r, "LogicalNeuroVol")
})

# ---------------------------------------------------------------------------
# Logical NOT (!)
# ---------------------------------------------------------------------------
test_that("!LogicalNeuroVol returns LogicalNeuroVol with inverted values", {
  lvol <- make_logical()
  r <- !lvol
  expect_s4_class(r, "LogicalNeuroVol")
  expect_equal(r@.Data, !lvol@.Data)
})

test_that("!DenseNeuroVol returns LogicalNeuroVol", {
  vol <- make_dense(sample(0:1, 125, replace = TRUE))
  r <- !vol
  expect_s4_class(r, "LogicalNeuroVol")
})

test_that("!SparseNeuroVol returns LogicalNeuroVol", {
  svol <- make_sparse()
  r <- !svol
  expect_s4_class(r, "LogicalNeuroVol")
  # non-zero voxels should become FALSE
  expect_false(r[1, 1, 1])
  # zero voxels should become TRUE
  expect_true(r[2, 1, 1])
})

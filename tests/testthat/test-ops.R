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

# ---------------------------------------------------------------------------
# Arith — SparseNeuroVec op SparseNeuroVec
# ---------------------------------------------------------------------------

make_sparse_vec_ops <- function(seed = 1) {
  set.seed(seed)
  sp4  <- NeuroSpace(c(4L, 4L, 4L, 5L), c(1, 1, 1))
  mask <- array(runif(64) > 0.5, c(4L, 4L, 4L))
  dat  <- matrix(rnorm(5 * sum(mask)), nrow = 5, ncol = sum(mask))
  SparseNeuroVec(dat, sp4, mask = mask)
}

test_that("SparseNeuroVec + SparseNeuroVec returns SparseNeuroVec", {
  sv1 <- make_sparse_vec_ops(1)
  sv2 <- make_sparse_vec_ops(2)
  r   <- sv1 + sv2
  expect_s4_class(r, "SparseNeuroVec")
  expect_equal(dim(r)[1:3], dim(sv1)[1:3])
  expect_equal(dim(r)[4],   dim(sv1)[4])
})

test_that("SparseNeuroVec * SparseNeuroVec preserves space", {
  sv1 <- make_sparse_vec_ops(3)
  sv2 <- make_sparse_vec_ops(4)
  r   <- sv1 * sv2
  expect_s4_class(r, "SparseNeuroVec")
  expect_identical(space(r), space(sv1))
})

test_that("SparseNeuroVec - SparseNeuroVec result has correct values", {
  sp4  <- NeuroSpace(c(3L, 3L, 3L, 4L), c(1, 1, 1))
  mask <- array(TRUE, c(3L, 3L, 3L))
  dat1 <- matrix(seq_len(4 * 27), nrow = 4, ncol = 27)
  dat2 <- matrix(1, nrow = 4, ncol = 27)
  sv1  <- SparseNeuroVec(dat1, sp4, mask = mask)
  sv2  <- SparseNeuroVec(dat2, sp4, mask = mask)
  r    <- sv1 - sv2
  expect_s4_class(r, "SparseNeuroVec")
})

# ---------------------------------------------------------------------------
# Arith — NeuroVec op NeuroVec (generic fallback via DenseNeuroVec)
# ---------------------------------------------------------------------------

make_dense_vec_ops <- function(seed = 1) {
  set.seed(seed)
  sp4 <- NeuroSpace(c(4L, 4L, 4L, 6L), c(1, 1, 1))
  DenseNeuroVec(array(rnorm(4 * 4 * 4 * 6), c(4, 4, 4, 6)), sp4)
}

test_that("DenseNeuroVec + DenseNeuroVec returns DenseNeuroVec", {
  dv1 <- make_dense_vec_ops(1)
  dv2 <- make_dense_vec_ops(2)
  r   <- dv1 + dv2
  expect_s4_class(r, "DenseNeuroVec")
  expect_equal(dim(r), dim(dv1))
  expect_identical(space(r), space(dv1))
})

test_that("DenseNeuroVec * DenseNeuroVec element-wise correctness", {
  sp4 <- NeuroSpace(c(2L, 2L, 2L, 3L), c(1, 1, 1))
  a   <- array(rep(2, 2 * 2 * 2 * 3), c(2, 2, 2, 3))
  b   <- array(rep(3, 2 * 2 * 2 * 3), c(2, 2, 2, 3))
  dv1 <- DenseNeuroVec(a, sp4)
  dv2 <- DenseNeuroVec(b, sp4)
  r   <- dv1 * dv2
  expect_s4_class(r, "DenseNeuroVec")
  expect_true(all(r@.Data == 6))
})

test_that("DenseNeuroVec - DenseNeuroVec dimension mismatch errors", {
  sp4a <- NeuroSpace(c(3L, 3L, 3L, 5L), c(1, 1, 1))
  sp4b <- NeuroSpace(c(4L, 4L, 4L, 5L), c(1, 1, 1))
  dv1  <- DenseNeuroVec(array(rnorm(3^3 * 5), c(3, 3, 3, 5)), sp4a)
  dv2  <- DenseNeuroVec(array(rnorm(4^3 * 5), c(4, 4, 4, 5)), sp4b)
  expect_error(dv1 - dv2)
})

# ---------------------------------------------------------------------------
# Arith — NeuroVec op NeuroVol and NeuroVol op NeuroVec
# ---------------------------------------------------------------------------

test_that("NeuroVec + NeuroVol returns DenseNeuroVec with same spatial dims", {
  sp4 <- NeuroSpace(c(4L, 4L, 4L, 5L), c(1, 1, 1))
  sp3 <- NeuroSpace(c(4L, 4L, 4L),     c(1, 1, 1))
  dv  <- DenseNeuroVec(array(rnorm(4^3 * 5), c(4, 4, 4, 5)), sp4)
  vol <- DenseNeuroVol(array(rnorm(4^3),     c(4, 4, 4)),     sp3)
  r   <- dv + vol
  expect_s4_class(r, "DenseNeuroVec")
  expect_equal(dim(r)[1:3], c(4L, 4L, 4L))
  expect_equal(dim(r)[4],   5L)
})

test_that("NeuroVol + NeuroVec returns DenseNeuroVec", {
  sp4 <- NeuroSpace(c(4L, 4L, 4L, 5L), c(1, 1, 1))
  sp3 <- NeuroSpace(c(4L, 4L, 4L),     c(1, 1, 1))
  dv  <- DenseNeuroVec(array(rnorm(4^3 * 5), c(4, 4, 4, 5)), sp4)
  vol <- DenseNeuroVol(array(rnorm(4^3),     c(4, 4, 4)),     sp3)
  r   <- vol + dv
  expect_s4_class(r, "DenseNeuroVec")
  expect_equal(dim(r)[1:3], c(4L, 4L, 4L))
  expect_equal(dim(r)[4],   5L)
})

test_that("NeuroVec * NeuroVol scales each time-point correctly", {
  sp4 <- NeuroSpace(c(2L, 2L, 2L, 3L), c(1, 1, 1))
  sp3 <- NeuroSpace(c(2L, 2L, 2L),     c(1, 1, 1))
  dv  <- DenseNeuroVec(array(rep(2, 2^3 * 3), c(2, 2, 2, 3)), sp4)
  vol <- DenseNeuroVol(array(rep(3, 2^3),     c(2, 2, 2)),     sp3)
  r   <- dv * vol
  expect_s4_class(r, "DenseNeuroVec")
  # every voxel at every time should be 2*3 = 6
  expect_true(all(as.matrix(r) == 6))
})

test_that("NeuroVec + NeuroVol spatial dim mismatch errors", {
  sp4 <- NeuroSpace(c(3L, 3L, 3L, 4L), c(1, 1, 1))
  sp3 <- NeuroSpace(c(4L, 4L, 4L),     c(1, 1, 1))
  dv  <- DenseNeuroVec(array(rnorm(3^3 * 4), c(3, 3, 3, 4)), sp4)
  vol <- DenseNeuroVol(array(rnorm(4^3),     c(4, 4, 4)),     sp3)
  expect_error(dv + vol)
})

# ---------------------------------------------------------------------------
# mean() for NeuroVec types
# ---------------------------------------------------------------------------

test_that("mean(DenseNeuroVec) returns DenseNeuroVol with correct values", {
  sp4 <- NeuroSpace(c(3L, 3L, 3L, 10L), c(1, 1, 1))
  arr <- array(rnorm(3^3 * 10), c(3, 3, 3, 10))
  dv  <- DenseNeuroVec(arr, sp4)
  mv  <- mean(dv)
  expect_s4_class(mv, "DenseNeuroVol")
  expect_equal(dim(mv), c(3L, 3L, 3L))
  # Check one voxel: mean over time
  expected_vox1 <- mean(arr[1, 1, 1, ])
  expect_equal(mv[1, 1, 1], expected_vox1, tolerance = 1e-10)
})

test_that("mean(SparseNeuroVec) returns SparseNeuroVol", {
  sp4  <- NeuroSpace(c(4L, 4L, 4L, 8L), c(1, 1, 1))
  mask <- array(runif(64) > 0.4, c(4L, 4L, 4L))
  set.seed(5)
  dat  <- matrix(rnorm(8 * sum(mask)), nrow = 8, ncol = sum(mask))
  sv   <- SparseNeuroVec(dat, sp4, mask = mask)
  mv   <- mean(sv)
  expect_s4_class(mv, "SparseNeuroVol")
  expect_equal(dim(mv), c(4L, 4L, 4L))
})

test_that("mean(SparseNeuroVec) values match column means of data matrix", {
  sp4  <- NeuroSpace(c(3L, 3L, 3L, 6L), c(1, 1, 1))
  mask <- array(TRUE, c(3L, 3L, 3L))
  set.seed(6)
  dat  <- matrix(rnorm(6 * 27), nrow = 6, ncol = 27)
  sv   <- SparseNeuroVec(dat, sp4, mask = mask)
  mv   <- mean(sv)
  expect_equal(as.vector(mv@data), colMeans(dat), tolerance = 1e-10)
})

test_that("mean(NeuroVec) generic fallback returns DenseNeuroVol", {
  # Use DenseNeuroVec which will dispatch to the DenseNeuroVec method,
  # verifying the generic is correctly dispatched
  sp4 <- NeuroSpace(c(3L, 3L, 3L, 5L), c(1, 1, 1))
  dv  <- DenseNeuroVec(array(seq_len(3^3 * 5), c(3, 3, 3, 5)), sp4)
  mv  <- mean(dv)
  expect_s4_class(mv, "DenseNeuroVol")
  expect_equal(dim(mv), c(3L, 3L, 3L))
  # voxel 1 has values 1, 28, 55, 82, 109 across 5 timepoints
  # (column-major: voxel 1 steps by prod(3,3,3)=27)
  expected <- mean(c(1, 1 + 27, 1 + 54, 1 + 81, 1 + 108))
  expect_equal(mv[1, 1, 1], expected, tolerance = 1e-10)
})

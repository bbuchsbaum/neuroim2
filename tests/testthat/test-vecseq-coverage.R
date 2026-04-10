library(testthat)
library(neuroim2)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_test_vec <- function(ntime = 6L) {
  sp <- NeuroSpace(c(5L, 5L, 5L, ntime), spacing = c(1, 1, 1))
  DenseNeuroVec(array(rnorm(5 * 5 * 5 * ntime), c(5, 5, 5, ntime)), sp)
}

# ---------------------------------------------------------------------------
# NeuroVecSeq constructor
# ---------------------------------------------------------------------------

test_that("NeuroVecSeq from varargs of DenseNeuroVec", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(6L)
  vs <- NeuroVecSeq(v1, v2)
  expect_s4_class(vs, "NeuroVecSeq")
  expect_equal(dim(vs)[4], 10L)
  expect_equal(dim(vs)[1:3], c(5L, 5L, 5L))
})

test_that("NeuroVecSeq from list via do.call", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs <- do.call(NeuroVecSeq, list(v1, v2))
  expect_s4_class(vs, "NeuroVecSeq")
  expect_equal(dim(vs)[4], 8L)
})

test_that("NeuroVecSeq errors if non-NeuroVec objects passed", {
  v1 <- make_test_vec(4L)
  expect_error(NeuroVecSeq(v1, "not_a_vec"))
})

test_that("NeuroVecSeq errors on spatially incompatible vecs", {
  sp_a <- NeuroSpace(c(5L, 5L, 5L, 4L))
  sp_b <- NeuroSpace(c(6L, 6L, 6L, 4L))
  va <- DenseNeuroVec(array(0, c(5, 5, 5, 4)), sp_a)
  vb <- DenseNeuroVec(array(0, c(6, 6, 6, 4)), sp_b)
  expect_error(NeuroVecSeq(va, vb))
})

test_that("NeuroVecSeq length equals sum of time dims", {
  v1 <- make_test_vec(3L)
  v2 <- make_test_vec(5L)
  v3 <- make_test_vec(2L)
  vs <- NeuroVecSeq(v1, v2, v3)
  expect_equal(dim(vs)[4], 10L)
})

test_that("as.matrix on NeuroVecSeq matches concatenated NeuroVec", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(6L)
  vs <- NeuroVecSeq(v1, v2)
  ref <- concat(v1, v2)

  expect_equal(as.matrix(vs), as.matrix(ref))
  expect_equal(as(vs, "matrix"), as.matrix(ref))
})

test_that("NeuroVecSeq can be materialized as DenseNeuroVec", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(6L)
  vs <- NeuroVecSeq(v1, v2)
  ref <- concat(v1, v2)

  dense <- as.dense(vs)
  coerced <- as(vs, "DenseNeuroVec")
  plain <- as(vs, "NeuroVec")

  expect_s4_class(dense, "DenseNeuroVec")
  expect_s4_class(coerced, "DenseNeuroVec")
  expect_s4_class(plain, "DenseNeuroVec")
  expect_false(inherits(plain, "NeuroVecSeq"))
  expect_equal(dim(dense), dim(ref))
  expect_equal(space(dense), space(ref))
  expect_equal(as.matrix(dense), as.matrix(ref))
  expect_equal(as.matrix(coerced), as.matrix(ref))
  expect_equal(as.matrix(plain), as.matrix(ref))
})

# ---------------------------------------------------------------------------
# [[ extraction
# ---------------------------------------------------------------------------

test_that("[[ extracts the correct volume from NeuroVecSeq", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs <- NeuroVecSeq(v1, v2)
  # Volume 1 should match v1[[1]]
  vol1_seq  <- vs[[1]]
  vol1_ref  <- v1[[1]]
  expect_equal(as.array(vol1_seq), as.array(vol1_ref))
})

test_that("[[ crosses segment boundary correctly", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs <- NeuroVecSeq(v1, v2)
  # Volume 5 is the first volume of v2
  vol5_seq <- vs[[5]]
  vol5_ref <- v2[[1]]
  expect_equal(as.array(vol5_seq), as.array(vol5_ref))
})

test_that("[[ errors on out-of-bounds index", {
  vs <- NeuroVecSeq(make_test_vec(4L))
  expect_error(vs[[5]])
  expect_error(vs[[0]])
})

# ---------------------------------------------------------------------------
# sub_vector
# ---------------------------------------------------------------------------

test_that("sub_vector on NeuroVecSeq returns correct number of volumes", {
  v1 <- make_test_vec(6L)
  v2 <- make_test_vec(6L)
  vs <- NeuroVecSeq(v1, v2)
  sv <- sub_vector(vs, 1:4)
  expect_s4_class(sv, "NeuroVecSeq")
  expect_equal(dim(sv)[4], 4L)
})

test_that("sub_vector spanning two segments returns correct data", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  ref <- concat(v1, v2)
  vs  <- NeuroVecSeq(v1, v2)
  # Indices 3:6 span both segments
  sv  <- sub_vector(vs, 3:6)
  sv_ref <- sub_vector(ref, 3:6)
  expect_equal(dim(sv)[4], 4L)
  expect_equal(as.array(sv[[1]]), as.array(sv_ref[[1]]))
})

test_that("sub_vector errors when index exceeds 4th dim", {
  vs <- NeuroVecSeq(make_test_vec(4L))
  expect_error(sub_vector(vs, 1:10))
})

# ---------------------------------------------------------------------------
# series extraction
# ---------------------------------------------------------------------------

test_that("series on NeuroVecSeq matches concat reference for voxel indices", {
  v1 <- make_test_vec(6L)
  v2 <- make_test_vec(6L)
  ref <- concat(v1, v2)
  vs  <- NeuroVecSeq(v1, v2)
  idx <- c(1L, 10L, 25L, 50L, 125L)
  mat_seq <- series(vs, idx)
  mat_ref <- series(ref, idx)
  expect_equal(mat_seq, mat_ref)
})

test_that("series on NeuroVecSeq with single index and drop=TRUE returns vector", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs  <- NeuroVecSeq(v1, v2)
  result <- series(vs, 1L, drop = TRUE)
  expect_true(is.numeric(result))
  expect_true(is.vector(result))
  expect_equal(length(result), 8L)
})

test_that("series on NeuroVecSeq with i,j,k coords returns correct length", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs  <- NeuroVecSeq(v1, v2)
  result <- series(vs, 2L, 3L, 4L)
  expect_equal(length(result), 8L)
})

test_that("series on NeuroVecSeq matches concat for i,j,k coords", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  ref <- concat(v1, v2)
  vs  <- NeuroVecSeq(v1, v2)
  result_seq <- series(vs, 2L, 3L, 4L)
  result_ref <- series(ref, 2L, 3L, 4L)
  expect_equal(result_seq, result_ref)
})

test_that("series on NeuroVecSeq with matrix coords", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  ref <- concat(v1, v2)
  vs  <- NeuroVecSeq(v1, v2)
  coords <- matrix(c(1L, 2L, 3L,
                     2L, 3L, 4L), nrow = 2, byrow = TRUE)
  mat_seq <- series(vs, coords)
  mat_ref <- series(ref, coords)
  expect_equal(dim(mat_seq), dim(mat_ref))
  expect_equal(mat_seq, mat_ref)
})

# ---------------------------------------------------------------------------
# linear_access
# ---------------------------------------------------------------------------

test_that("linear_access on NeuroVecSeq matches concat reference", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  ref <- concat(v1, v2)
  vs  <- NeuroVecSeq(v1, v2)
  set.seed(42)
  idx <- sample(seq_len(prod(dim(vs))), 30)
  expect_equal(linear_access(vs, idx), linear_access(ref, idx))
})

# ---------------------------------------------------------------------------
# vectors()
# ---------------------------------------------------------------------------

test_that("vectors() on NeuroVecSeq returns deflist of correct length", {
  v1 <- make_test_vec(4L)
  vs <- NeuroVecSeq(v1)
  vl <- vectors(vs)
  expect_equal(length(vl), prod(dim(vs)[1:3]))
})

test_that("vectors() on NeuroVecSeq with numeric subset", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(4L)
  vs <- NeuroVecSeq(v1, v2)
  ref <- concat(v1, v2)
  subset_idx <- c(1L, 10L, 25L)
  vl_seq <- vectors(vs, subset = subset_idx)
  vl_ref <- vectors(ref, subset = subset_idx)
  expect_equal(length(vl_seq), 3L)
  for (i in seq_along(subset_idx)) {
    expect_equal(vl_seq[[i]], vl_ref[[i]])
  }
})

test_that("vectors() on NeuroVecSeq with logical subset", {
  v1 <- make_test_vec(4L)
  vs <- NeuroVecSeq(v1)
  nvox <- prod(dim(vs)[1:3])
  mask_lgl <- rep(FALSE, nvox)
  mask_lgl[c(1L, 5L, 10L)] <- TRUE
  vl <- vectors(vs, subset = mask_lgl)
  expect_equal(length(vl), 3L)
})

# ---------------------------------------------------------------------------
# concat produces equivalent result to NeuroVecSeq
# ---------------------------------------------------------------------------

test_that("NeuroVecSeq space matches concat space", {
  v1 <- make_test_vec(4L)
  v2 <- make_test_vec(6L)
  vs  <- NeuroVecSeq(v1, v2)
  ref <- concat(v1, v2)
  expect_equal(dim(vs), dim(ref))
  expect_equal(space(vs), space(ref))
})

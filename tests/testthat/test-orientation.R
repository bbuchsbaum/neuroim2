context("orientation")

library(neuroim2)

# Test findAnatomy3D with explicit axis abbreviations

test_that("findAnatomy3D returns the expected orientation", {
  orient <- findAnatomy3D("L", "P", "I")
  expect_true(inherits(orient, "AxisSet3D"))
  expect_identical(orient@i, LEFT_RIGHT)
  expect_identical(orient@j, POST_ANT)
  expect_identical(orient@k, INF_SUP)
})

# Test findAnatomy on a simple permutation matrix

test_that("findAnatomy recovers orientation from identity matrix", {
  pmat <- diag(3)
  orient <- neuroim2:::findAnatomy(pmat)
  expect_identical(orient@i, LEFT_RIGHT)
  expect_identical(orient@j, POST_ANT)
  expect_identical(orient@k, INF_SUP)
})

# ---- axcodes on NeuroSpace ----

test_that("axcodes(NeuroSpace) returns character vector of length 3", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
  codes <- axcodes(sp)
  expect_type(codes, "character")
  expect_length(codes, 3)
})

test_that("axcodes(NeuroSpace) values are valid anatomical codes", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
  codes <- axcodes(sp)
  valid <- c("R", "L", "A", "P", "S", "I")
  expect_true(all(codes %in% valid))
})

test_that("axcodes(NeuroSpace) RAS affine gives R A S codes", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(10L, 10L, 10L), trans = aff)
  codes <- axcodes(sp)
  expect_equal(codes, c("R", "A", "S"))
})

test_that("axcodes(NeuroSpace) LPS affine gives L P S codes", {
  aff <- diag(4)
  aff[1, 1] <- -2; aff[2, 2] <- -2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(10L, 10L, 10L), trans = aff)
  codes <- axcodes(sp)
  expect_equal(codes, c("L", "P", "S"))
})

# ---- axcodes on NeuroObj (DenseNeuroVol) ----

test_that("axcodes(DenseNeuroVol) returns character vector of length 3", {
  vol <- make_vol(dim = c(8L, 8L, 8L))
  codes <- axcodes(vol)
  expect_type(codes, "character")
  expect_length(codes, 3)
})

test_that("axcodes(DenseNeuroVol) agrees with axcodes(space(vol))", {
  vol <- make_vol(dim = c(8L, 8L, 8L))
  expect_equal(axcodes(vol), axcodes(space(vol)))
})

# ---- axcodes on matrix ----

test_that("axcodes(matrix) works on a 4x4 RAS affine", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  codes <- axcodes(aff)
  expect_type(codes, "character")
  expect_length(codes, 3)
  expect_equal(codes, c("R", "A", "S"))
})

test_that("axcodes(matrix) result matches axcodes(NeuroSpace) built from same affine", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(10L, 10L, 10L), trans = aff)
  expect_equal(axcodes(aff), axcodes(sp))
})

# ---- as_canonical ----

test_that("as_canonical returns same object when already RAS", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(6L, 6L, 6L), trans = aff)
  vol <- DenseNeuroVol(array(rnorm(216), c(6L, 6L, 6L)), sp)
  result <- as_canonical(vol)
  expect_identical(result, vol)
})

test_that("as_canonical actually reorients a non-RAS volume", {
  # Affine (-2,-2,2) => LPS; as_canonical flips L->R and P->A
  aff <- diag(4)
  aff[1, 1] <- -2; aff[2, 2] <- -2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(6L, 6L, 6L), trans = aff)
  vol <- DenseNeuroVol(array(rnorm(216), c(6L, 6L, 6L)), sp)
  original_codes <- axcodes(vol)
  expect_false(identical(original_codes, c("R", "A", "S")))
  result <- as_canonical(vol)
  # The result must differ from original orientation
  expect_false(identical(axcodes(result), original_codes))
  # R and A axes should be canonical (x=R, y=A)
  result_codes <- axcodes(result)
  expect_equal(result_codes[1], "R")
  expect_equal(result_codes[2], "A")
})

test_that("as_canonical result is a DenseNeuroVol", {
  aff <- diag(4)
  aff[1, 1] <- -2; aff[2, 2] <- -2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(6L, 6L, 6L), trans = aff)
  vol <- DenseNeuroVol(array(rnorm(216), c(6L, 6L, 6L)), sp)
  result <- as_canonical(vol)
  expect_s4_class(result, "DenseNeuroVol")
})

test_that("as_canonical preserves total voxel count", {
  aff <- diag(4)
  aff[1, 1] <- -2; aff[2, 2] <- -2; aff[3, 3] <- 2
  sp <- NeuroSpace(c(6L, 7L, 8L), trans = aff)
  vol <- DenseNeuroVol(array(rnorm(6 * 7 * 8), c(6L, 7L, 8L)), sp)
  result <- as_canonical(vol)
  expect_equal(prod(dim(result)), 6L * 7L * 8L)
})

test_that("as_canonical on already-RAS NeuroVec returns same object", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  sp4 <- NeuroSpace(c(4L, 4L, 4L, 5L), trans = aff)
  vec <- DenseNeuroVec(array(rnorm(4 * 4 * 4 * 5), c(4L, 4L, 4L, 5L)), sp4)
  result <- as_canonical(vec)
  expect_identical(result, vec)
})

context("orientation additional coverage")

library(neuroim2)

# ---- as_canonical error path ----

test_that("as_canonical errors on non-NeuroVol/NeuroVec object", {
  expect_error(as_canonical("not an image"), class = "error")
})

test_that("as_canonical errors on a plain numeric vector", {
  expect_error(as_canonical(1:10), class = "error")
})

# ---- as_canonical already-RAS NeuroVec (identity branch) ----

test_that("as_canonical already-RAS NeuroVec returns the same object", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  sp4 <- NeuroSpace(c(4L, 4L, 4L, 3L), trans = aff)
  vec <- DenseNeuroVec(array(rnorm(4 * 4 * 4 * 3), c(4L, 4L, 4L, 3L)), sp4)
  result <- as_canonical(vec)
  expect_identical(result, vec)
})

# ---- axcodes on matrix via as_canonical / axcodes paths ----

test_that("affine_to_axcodes round-trips through axcodes(matrix)", {
  aff <- diag(4)
  aff[1, 1] <- 2; aff[2, 2] <- 2; aff[3, 3] <- 2
  codes <- axcodes(aff)
  expect_equal(codes, c("R", "A", "S"))
})

test_that("axcodes on negative-diagonal affine gives L P I", {
  aff <- diag(4)
  aff[1, 1] <- -2; aff[2, 2] <- -2; aff[3, 3] <- -2
  codes <- axcodes(aff)
  expect_equal(codes, c("L", "P", "I"))
})

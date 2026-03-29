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

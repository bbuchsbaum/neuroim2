context("neurovol regressions")

library(neuroim2)

# Helpers
make_space <- function(d) NeuroSpace(d, spacing=rep(1,length(d)))

test_that("SparseNeuroVol as.array preserves indices", {
  sp <- make_space(c(3,3,3))
  idx <- c(2, 5, 10)
  vals <- c(1, 2, 3)
  sv <- SparseNeuroVol(vals, sp, indices = idx)
  arr <- as(sv, "array")
  expect_equal(arr[idx], vals)
  expect_true(all(arr[-idx] == 0))
})

test_that("SparseNeuroVol rejects out-of-range indices", {
  sp <- make_space(c(2,2,2))
  expect_error(SparseNeuroVol(1, sp, indices = 9), "within 1..prod")
})

test_that("map_values with numeric list keys works and character keys fail", {
  sp <- make_space(c(2,2,1))
  v <- NeuroVol(array(c(1,2,1,2), dim=c(2,2,1)), sp)
  lookup <- list(`1`=10, `2`=20)
  out <- map_values(v, lookup)
  expect_s4_class(out, "NeuroVol")
  expect_equal(unique(as.vector(out)), c(10,20))
  bad_lookup <- list(a=1, b=2)
  expect_error(map_values(v, bad_lookup))
})

test_that("DenseNeuroVol as.matrix exists and matches data", {
  sp <- make_space(c(2,2,1))
  arr <- array(1:4, dim=c(2,2,1))
  v <- DenseNeuroVol(arr, sp)
  m <- as.matrix(v)
  expect_equal(m, matrix(arr, ncol=1))
})

test_that("mapf produces expected size and respects mask", {
  sp <- make_space(c(5,5,5))
  arr <- array(1, dim=c(5,5,5))
  v <- NeuroVol(arr, sp)
  kdim <- c(3,3,3)
  ker <- Kernel(kerndim=kdim, vdim=c(1,1,1), FUN=function(d) as.numeric(d==0)) # center =1, else 0
  msk <- LogicalNeuroVol(array(FALSE, dim=c(5,5,5)), sp)
  msk@.Data[3,3,3] <- TRUE
  out <- mapf(v, ker, mask = msk)
  expect_equal(dim(out), dim(v))
  expect_equal(out[3,3,3], 1)
  expect_true(all(out[-c(3+5*(3-1)+25*(3-1))] == 0))
})

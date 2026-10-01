test_that("HyperVec access preserves independent trial and feature coordinates", {
  d <- c(2L, 2L, 1L, 3L, 2L)
  mask <- LogicalNeuroVol(array(c(TRUE, FALSE, TRUE, TRUE), d[1:3]),
                          NeuroSpace(d[1:3]))
  ref <- array(0, d)
  dat <- array(0, c(d[5], d[4], 3L))
  for (v in seq_len(3)) for (t in seq_len(d[4])) for (f in seq_len(d[5])) {
    value <- 100 * f + 10 * t + which(mask)[v]
    dat[f, t, v] <- value
    ref[which(mask)[v] + (t - 1) * 4 + (f - 1) * 12] <- value
  }
  x <- NeuroHyperVec(dat, NeuroSpace(d), mask)
  ind <- c(24L, 2L, 1L, 13L, 1L)
  expect_equal(linear_access(x, ind), as.vector(ref)[ind])
  expect_equal(x[c(2, 1), , , c(3, 1, 3), c(2, 1), drop = FALSE],
               ref[c(2, 1), , , c(3, 1, 3), c(2, 1), drop = FALSE])
  expect_equal(series(x, 1, 2, 1), t(ref[1, 2, 1, , ]))
  expect_equal(series(x, 2, 1, 1), matrix(0, 2, 3))
  expect_error(series(x, 1L), "Must provide")
  expect_error(as.matrix(x))
  expect_error(split_reduce(x, factor(1:4)))
})

test_that("clustered access broadcasts cluster values with outside-mask NA", {
  sp <- NeuroSpace(c(2, 2, 1))
  mask <- LogicalNeuroVol(array(c(TRUE, FALSE, TRUE, TRUE), dim(sp)), sp)
  clusters <- ClusteredNeuroVol(mask, c(1L, 1L, 2L))
  arr <- array(as.numeric(1:12), c(2, 2, 1, 3))
  x <- ClusteredNeuroVec(DenseNeuroVec(arr, NeuroSpace(dim(arr))), clusters)
  expected <- cbind(c(2, 6, 10), c(4, 8, 12))
  expect_equal(unname(as.matrix(x)), expected)
  expect_equal(series(x, c(1, 2, 1), c(1, 1, 2), c(1, 1, 1)),
               cbind(expected[, 1], rep(NA_real_, 3), expected[, 1]))
  expect_equal(x[c(2, 1, 1), c(1, 2, 1), c(1, 1, 1), c(3, 1), drop = FALSE],
               cbind(rep(NA_real_, 2), expected[c(3, 1), 1], expected[c(3, 1), 1]))
  expect_error(linear_access(x, 1L))
  expect_error(series(x, 1L))
  expect_error(as.matrix(x, by = "voxel"), "not yet implemented")
  expect_error(split_reduce(x, factor(1:4)))
})

test_that("NeuroBucket inheritance does not advertise working data access", {
  sp <- NeuroSpace(c(2, 2, 1))
  vol <- NeuroVol(array(1:4, dim(sp)), sp)
  x <- new("NeuroBucket", data = list(vol, vol), labels = c("a", "b"),
           space = add_dim(sp, 2))
  expect_error(linear_access(x, 1L))
  expect_error(x[1L])
  expect_error(series(x, 1L))
  expect_error(as.matrix(x))
  expect_error(split_reduce(x, factor(1:4)))
})

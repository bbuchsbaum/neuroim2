library(testthat)
library(neuroim2)

context("lazy sparse accessors")

register_mock_lazy_sparse_neurovec <- function() {
  if (!methods::isClass("MockLazySparseNeuroVec")) {
    methods::setClass(
      "MockLazySparseNeuroVec",
      slots = c(data = "matrix", backing = "matrix"),
      contains = c("AbstractSparseNeuroVec")
    )
  }

  if (!methods::hasMethod("matricized_access", c("MockLazySparseNeuroVec", "integer"))) {
    methods::setMethod(
      "matricized_access",
      signature(x = "MockLazySparseNeuroVec", i = "integer"),
      function(x, i, ...) x@backing[, i, drop = FALSE]
    )
  }

  if (!methods::hasMethod("matricized_access", c("MockLazySparseNeuroVec", "numeric"))) {
    methods::setMethod(
      "matricized_access",
      signature(x = "MockLazySparseNeuroVec", i = "numeric"),
      function(x, i, ...) x@backing[, i, drop = FALSE]
    )
  }

  if (!methods::hasMethod("temporal_access", c("MockLazySparseNeuroVec", "integer"))) {
    methods::setMethod(
      "temporal_access",
      signature(x = "MockLazySparseNeuroVec", i = "integer"),
      function(x, i, ...) x@backing[i, , drop = FALSE]
    )
  }

  if (!methods::hasMethod("temporal_access", c("MockLazySparseNeuroVec", "numeric"))) {
    methods::setMethod(
      "temporal_access",
      signature(x = "MockLazySparseNeuroVec", i = "numeric"),
      function(x, i, ...) x@backing[i, , drop = FALSE]
    )
  }
}

make_mock_lazy_sparse_vec <- function() {
  register_mock_lazy_sparse_neurovec()

  dims <- c(2, 2, 1, 3)
  sp <- NeuroSpace(dims)
  mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = dims[1:3])
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(dims[1:3]))
  backing <- matrix(c(10, 20, 30,
                      40, 50, 60), nrow = dims[4], byrow = FALSE)
  placeholder <- matrix(-999, nrow = dims[4], ncol = sum(mask_arr))

  dense_arr <- array(0, dim = dims)
  dense_arr[1, 1, 1, ] <- backing[, 1]
  dense_arr[1, 2, 1, ] <- backing[, 2]

  list(
    lazy = methods::new(
      "MockLazySparseNeuroVec",
      space = sp,
      label = "mock-lazy",
      mask = mask,
      map = IndexLookupVol(space(mask), as.integer(which(mask_arr))),
      data = placeholder,
      backing = backing
    ),
    dense = DenseNeuroVec(dense_arr, sp)
  )
}

test_that("lazy sparse subclasses use accessor hooks for series and [", {
  objs <- make_mock_lazy_sparse_vec()
  x <- objs$lazy
  dense <- objs$dense

  expect_equal(series(x, 1L), c(10, 20, 30))
  expect_equal(series(x, 1L), series(dense, 1L))
  expect_equal(series(x, 1, 1, 1), series(dense, 1, 1, 1))
  expect_equal(x[1, 1, 1, 1], dense[1, 1, 1, 1])
  expect_equal(
    x[1:2, 1:2, 1, 1:3, drop = FALSE],
    dense[1:2, 1:2, 1, 1:3, drop = FALSE]
  )
})

test_that("lazy sparse subclasses use temporal access for [[, sub_vector, and as.matrix", {
  objs <- make_mock_lazy_sparse_vec()
  x <- objs$lazy
  dense <- objs$dense

  expect_equal(as.vector(x[[1]]), as.vector(dense[[1]]))
  expect_equal(as.matrix(sub_vector(x, 1:3)), as.matrix(sub_vector(dense, 1:3)))
  expect_equal(as.matrix(x), as.matrix(dense))
  expect_equal(as(x, "matrix"), as.matrix(dense))
})

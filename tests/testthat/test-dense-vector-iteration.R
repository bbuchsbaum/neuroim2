test_that("dense vector iteration agrees with independent array time courses", {
  shapes <- list(c(3L,4L,2L,5L), c(1L,1L,1L,7L), c(3L,2L,1L,1L),
                 c(1L,1L,1L,1L))
  for (d in shapes) {
    a <- array(seq_len(prod(d)) / 7, d)
    x <- DenseNeuroVec(a, NeuroSpace(d, spacing=c(2,3,4)))
    it <- vectors(x)
    expect_s3_class(it, "deflist")
    expect_length(it, prod(d[1:3]))
    grid <- arrayInd(seq_len(prod(d[1:3])), d[1:3])
    expected <- lapply(seq_len(nrow(grid)), function(i) {
      unname(a[grid[i,1], grid[i,2], grid[i,3], ])
    })
    # Access out of order and repeatedly before realising the full list.
    for (i in c(length(it), 1L, length(it))) {
      expect_identical(it[[i]], expected[[i]])
    }
    expect_equal(as.list(it), expected)
    expect_equal(vapply(it, mean, numeric(1)),
                 as.numeric(apply(a, 1:3, mean)))
    expect_identical(as.array(x), a)
  }
})

test_that("dense iterators retain missing values and the original value snapshot", {
  a <- array(c(NA_real_, NaN, Inf, -Inf, 0, 2, 3, 4, 5, 6, 7, 8),
             c(2L,1L,2L,3L))
  x <- DenseNeuroVec(a, NeuroSpace(dim(a)))
  it <- vectors(x)
  # Changing the caller's binding must not change the already-created iterator.
  x[1,1,1,1] <- 99
  expect_identical(it[[1]], c(NA_real_, 0, 5))
  expect_identical(it[[2]], c(NaN, 2, 6))
  expect_identical(it[[3]], c(Inf, 3, 7))
  expect_identical(it[[4]], c(-Inf, 4, 8))
})

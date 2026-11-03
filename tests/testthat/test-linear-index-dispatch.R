# vol[i] is linear indexing and vol[i, ] / vol[i, , k] index the first
# dimension, decided by the shape of the call rather than by the values of i.
# Before #60 a vector of linear indices that all fitted within dim(x)[1]
# returned whole slices on dense volumes, and sparse volumes never did linear
# indexing at all.

first_row_fixture <- function() {
  sp <- NeuroSpace(c(6, 6, 4))
  a <- array(0, c(6, 6, 4))
  a[6] <- 5   # first x-row: linear index <= dim(x)[1]
  a[8] <- 7
  list(
    sp = sp, a = a,
    vols = list(
      dense = DenseNeuroVol(a, sp),
      logical = LogicalNeuroVol(a > 0, sp),
      sparse = SparseNeuroVol(c(5, 7), sp, indices = c(6, 8)),
      clustered = ClusteredNeuroVol(LogicalNeuroVol(a > 0, sp), c(1L, 2L))
    ),
    expected = list(
      dense = c(5, 7), logical = c(TRUE, TRUE),
      sparse = c(5, 7), clustered = c(1L, 2L)
    )
  )
}

test_that("vol[i] is linear indexing even when every i fits the first dimension", {
  f <- first_row_fixture()
  for (nm in names(f$vols)) {
    v <- f$vols[[nm]]
    expect_equal(v[6], f$expected[[nm]][1], info = nm)
    expect_equal(v[6L], f$expected[[nm]][1], info = nm)
    expect_equal(v[c(6, 8)], f$expected[[nm]], info = nm)
    expect_equal(v[c(6L, 8L)], f$expected[[nm]], info = nm)
    expect_false(inherits(v[c(6, 8)], "sparseVector"), info = nm)
  }
})

test_that("vol[which(mask)] returns the masked values for a first-row mask", {
  f <- first_row_fixture()
  dv <- f$vols$dense
  m <- f$a > 0
  expect_equal(dv[which(m)], f$a[m])
  # Indices beyond the first dimension were already linear; they still are.
  expect_equal(dv[c(6, 8, 100)], f$a[c(6, 8, 100)])
})

test_that("vol[i, ] and vol[i, , k] still index the first dimension", {
  f <- first_row_fixture()
  for (nm in names(f$vols)) {
    v <- f$vols[[nm]]
    expect_equal(dim(v[6, , ]), c(6L, 4L), info = nm)
    expect_equal(dim(v[6, ]), c(6L, 4L), info = nm)
    expect_length(v[6, , 2], 6L)
    expect_equal(dim(v[c(1, 6), , 1]), c(2L, 6L), info = nm)
    expect_equal(dim(v[6, , 2, drop = FALSE]), c(1L, 6L, 1L), info = nm)
    # An explicit drop only makes sense for slicing, so it selects slicing.
    expect_equal(dim(v[6, drop = FALSE]), c(1L, 6L, 4L), info = nm)
  }
  dv <- f$vols$dense
  expect_equal(dv[6, , 2], f$a[6, , 2])
  expect_equal(dv[2:3, , ], f$a[2:3, , ])
})

test_that("conn_comp handles a cluster lying in the first x-row (#60)", {
  sp <- NeuroSpace(c(6, 6, 4))
  m <- array(FALSE, c(6, 6, 4))
  m[6] <- TRUE
  cc <- conn_comp(LogicalNeuroVol(m, sp))
  expect_equal(num_clusters(cc$index), 1L)
  expect_equal(which(as.array(cc$index) > 0), 6L)

  m[c(5, 6)] <- TRUE
  cc2 <- conn_comp(LogicalNeuroVol(m, sp), connect = "6-connect")
  expect_equal(num_clusters(cc2$index), 1L)
  expect_equal(sum(as.array(cc2$size) > 0), 2L)
})

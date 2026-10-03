# These six backends share voxel-series and voxel-matrix layouts. Independent
# binary fixtures catch decoding and square-matrix orientation regressions.
test_that("sparse downsampling keeps a square time-by-voxel result oriented", {
  mask <- array(FALSE, c(4, 4, 4))
  mask[c(1, 2, 3, 4)] <- TRUE
  # Two pairs of input voxels average to two output voxels. No symmetric or
  # constant matrix can expose an accidental transpose here.
  input <- rbind(c(1, 3, 5, 7), c(10, 30, 50, 70))
  x <- SparseNeuroVec(input, NeuroSpace(c(4, 4, 4, 2)), mask,
                     orientation = "time_x_voxels")
  expect_warning(y <- downsample(x, factor = 0.5), NA)
  expect_equal(series(y, indices(y)), rbind(c(2, 6), c(20, 60)))
})

test_that("sparse producers preserve non-symmetric square time-by-voxel values", {
  mask <- array(FALSE, c(2, 2, 2)); mask[c(1, 3, 6)] <- TRUE
  vox <- c(1L, 3L, 6L)
  input <- matrix(c(1, 3, 9, 10, 20, 50, 11, 30, 80), 3, 3)
  make <- function(m) SparseNeuroVec(m, NeuroSpace(c(2, 2, 2, nrow(m))),
                                    mask, orientation = "time_x_voxels")
  x <- make(input)
  expect_warning(scaled <- scale_series(x, TRUE, TRUE), NA)
  expect_equal(series(scaled, vox), unname(base::scale(input)))
  expect_warning(sum <- x + x, NA)
  expect_equal(series(sum, vox), input * 2)
  a <- make(input[1, , drop = FALSE])
  b <- make(input[2:3, , drop = FALSE])
  expect_warning(joined <- concat(a, b), NA)
  expect_equal(series(joined, vox), input)
  expect_warning(joined <- concat(a, make(input[2, , drop = FALSE]),
                                  make(input[3, , drop = FALSE])), NA)
  expect_equal(series(joined, vox), input)
  roi <- ROIVec(space(x), arrayInd(vox, c(2, 2, 2)), input)
  expect_warning(converted <- as(roi, "SparseNeuroVec"), NA)
  expect_equal(series(converted, vox), input)
})

for (encoding in c("FLOAT", "SHORT", "UBYTE")) {
  for (nt in c(3L, 6L)) {
    test_that(paste("backend accessor matrix", encoding, "times", nt), {
      f <- backend_binary_fixture(encoding, c(3L, 2L, 1L, nt))
      on.exit(unlink(f$path), add = TRUE)
      ref <- matrix(f$expected, 6, nt)
      dense <- read_vec(f$path, mode = "normal")
      mapped <- MappedNeuroVec(f$path)
      on.exit(mmap::munmap(mapped@filemap), add = TRUE, after = FALSE)
      fb <- FileBackedNeuroVec(f$path)
      fullmask <- array(TRUE, c(3, 2, 1))
      sparse <- SparseNeuroVec(f$expected, space(dense), fullmask)
      # bigstatsr is an Imports dependency; absence is a failure, not a skip.
      backing <- tempfile()
      on.exit(unlink(paste0(backing, c(".bk", ".rds"))), add = TRUE)
      big <- BigNeuroVec(f$expected, space(dense), fullmask, backingfile = backing)
      seq <- NeuroVecSeq(mapped, sparse, fb)
      backends <- list(dense = dense, sparse = sparse, big = big,
                       mapped = mapped, filebacked = fb, sequence = seq)
      expect_equal(as.matrix(dense), ref)
      for (name in names(backends)) {
        x <- backends[[name]]
        want <- if (name == "sequence") cbind(ref, ref, ref) else ref
        idx <- c(length(want), 1L, 7L, 6L, 1L)
        vox <- c(6L, 1L, 6L) # square for each three-timepoint component
        times <- c(ncol(want), 1L, nt, nt + (name == "sequence"))
        expect_equal(linear_access(x, idx), as.vector(want)[idx], info = name)
        expect_equal(x[idx], as.vector(want)[idx], info = name)
        expect_equal(as.matrix(x), want, info = name)
        expect_equal(series(x, vox, drop = FALSE), t(want[vox, , drop = FALSE]), info = name)
        expect_equal(series(x, 1:6, drop = FALSE), t(want), info = name)
        expect_equal(x[, , , times, drop = FALSE],
                     array(want[, times], c(3, 2, 1, length(times))), info = name)
        fac <- factor(rep(1:2, each = 3))
        reduced <- rbind(colMeans(want[1:3, , drop = FALSE]),
                         colMeans(want[4:6, , drop = FALSE]))
        rownames(reduced) <- levels(fac)
        expect_equal(split_reduce(x, fac), reduced, info = name)
      }
    })
  }
}

test_that("sparse, FBM and mixed sequence keep support and missing values distinct", {
  arr <- array(as.numeric(1:18), c(3, 2, 1, 3))
  arr[1] <- NA_real_
  mask <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE), c(3, 2, 1))
  sp <- NeuroSpace(dim(arr))
  sparse <- SparseNeuroVec(arr, sp, mask)
  backing <- tempfile()
  on.exit(unlink(paste0(backing, c(".bk", ".rds"))))
  big <- BigNeuroVec(arr, sp, mask, backingfile = backing)
  ref <- matrix(arr, 6, 3)
  ref[!as.vector(mask), ] <- 0
  for (x in list(sparse, big, NeuroVecSeq(sparse, big))) {
    want <- if (is(x, "NeuroVecSeq")) cbind(ref, ref) else ref
    expect_equal(as.matrix(x), want)
    expect_equal(series(x, c(2L, 1L, 2L), drop = FALSE), t(want[c(2, 1, 2), , drop = FALSE]))
    idx <- c(1, 2, length(want), 2)
    expect_equal(linear_access(x, idx), as.vector(want)[idx])
  }
})

# The applicability table is dev/backend-accessor-contracts.md. These six
# backends share voxel-series and voxel-matrix layouts; no shape inference.
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

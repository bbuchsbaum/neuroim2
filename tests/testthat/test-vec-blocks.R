block_test_sources <- function(arr, backing, labels = character()) {
  sp <- NeuroSpace(dim(arr), spacing = c(2, 3, 4))
  mask <- array(rep_len(c(TRUE, FALSE, TRUE), prod(dim(arr)[1:3])), dim(arr)[1:3])
  list(dense = DenseNeuroVec(arr, sp, volume_labels = labels),
       sparse = SparseNeuroVec(arr, sp, mask, volume_labels = labels),
       big = BigNeuroVec(arr, sp, mask, backingfile = backing, volume_labels = labels))
}

test_that("blocks reconstruct each supported store along either axis", {
  for (d in list(c(3, 2, 1, 5), c(1, 1, 1, 1), c(2, 2, 1, 4))) {
    arr <- array(as.numeric(seq_len(prod(d))), d)
    arr[1] <- NA_real_
    labels <- paste0("t", seq_len(d[4]))
    backing <- tempfile()
    sources <- block_test_sources(arr, backing, labels)
    for (name in names(sources)) for (axis in c("auto", "time", "voxel")) {
      x <- sources[[name]]
      chosen <- if (axis == "auto") if (name == "dense") "time" else "voxel" else axis
      width <- if (chosen == "time") prod(d[1:3]) else d[4]
      # About two units per block, deliberately leaving a final partial block.
      budget <- 8192 + 32 * width + 2 * (128 * width + 32) + 128
      blocks <- vec_blocks(x, axis, budget)
      plan <- attr(blocks, "plan")
      expect_identical(plan$along, chosen)
      expect_lte(plan$estimated_buffer_bytes, budget)
      expect_false(attr(blocks, "memoised"))
      want <- matrix(arr, prod(d[1:3]), d[4])
      if (name != "dense") want[!as.vector(x@mask), ] <- 0
      got <- matrix(NA_real_, nrow(want), ncol(want))
      visits <- matrix(0L, nrow(want), ncol(want))
      for (j in rev(seq_along(blocks))) {
        b <- blocks[[j]]
        expect_identical(dim(b$values), c(length(b$times), length(b$voxels)))
        expect_equal(b$values, t(want[b$voxels, b$times, drop = FALSE]))
        expect_equal(b$space, drop_dim(space(x)))
        expect_identical(b$volume_labels, labels[b$times])
        got[b$voxels, b$times] <- t(b$values)
        visits[b$voxels, b$times] <- visits[b$voxels, b$times] + 1L
      }
      expect_equal(got, want)
      expect_true(all(visits == 1L))
      expect_equal(blocks[[1]], blocks[[1]])
      expect_equal(as.matrix(x), want) # source remains unchanged
    }
    unlink(paste0(backing, c(".bk", ".rds")))
  }
})

test_that("construction is metadata-only, realization is lazy and not memoised", {
  arr <- array(as.numeric(1:120), c(3, 2, 1, 20))
  backing <- tempfile()
  on.exit(unlink(paste0(backing, c(".bk", ".rds"))))
  x <- block_test_sources(arr, backing)$sparse
  calls <- 0L
  original <- series
  testthat::local_mocked_bindings(series = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "neuroim2")
  blocks <- vec_blocks(x, budget = 12000)
  expect_identical(calls, 0L)
  before <- serialize(unclass(blocks)[seq_along(blocks)], NULL)
  b <- blocks[[1]]
  expect_identical(calls, 1L)
  expect_equal(blocks[[1]], b)
  expect_identical(calls, 2L)
  expect_identical(serialize(unclass(blocks)[seq_along(blocks)], NULL), before)
  # deflist's underlying cells remain NULL; no realized matrices are retained.
  expect_true(all(vapply(unclass(blocks), is.null, logical(1))))
})

test_that("bad budgets, unsupported classes and too-small blocks fail before reads", {
  x <- DenseNeuroVec(array(as.numeric(1:120), c(3, 2, 1, 20)),
                    NeuroSpace(c(3, 2, 1, 20)))
  for (budget in list(0, -1, NA_real_, Inf, numeric(), c(1, 2), "1000")) {
    expect_error(vec_blocks(x, budget = budget), "positive finite")
  }
  expect_error(vec_blocks(x, "other"))
  expect_error(vec_blocks(x, budget = 1), "minimum legal block")
  expect_error(vec_blocks(x, "voxel", budget = 9000), "minimum legal block")
  expect_error(vec_blocks(array(0, c(2, 2, 2, 2))), "does not yet support")
  expect_error(vec_blocks(NeuroVecSeq(x, x)), "NeuroVecSeq")
  blocks <- vec_blocks(x, budget = 12000)
  expect_error(blocks[[1.5]], "whole number")
  expect_error(blocks[[length(blocks) + 1]], "out of bounds")
  f <- backend_binary_fixture()
  on.exit(unlink(f$path), add = TRUE)
  mapped <- MappedNeuroVec(f$path)
  on.exit(mmap::munmap(mapped@filemap), add = TRUE, after = FALSE)
  expect_error(vec_blocks(mapped), "MappedNeuroVec")
  expect_error(vec_blocks(FileBackedNeuroVec(f$path)), "FileBackedNeuroVec")
})

test_that("fixed-budget realization does not allocate a full growing source", {
  # Rprofmem measures allocations, not RSS. Warm each path before measuring.
  budget <- 65536
  for (name in c("dense", "sparse", "big")) {
    maxima <- numeric()
    for (nv in c(128L, 4096L)) {
      arr <- array(as.numeric(seq_len(nv * 8L)), c(nv, 1L, 1L, 8L))
      backing <- tempfile()
      sources <- block_test_sources(arr, backing)
      x <- sources[[name]]
      # Voxel traversal fixes the minimum legal width as the source grows.
      blocks <- vec_blocks(x, along = "voxel", budget = budget)
      invisible(blocks[[1]])
      path <- tempfile()
      Rprofmem(path)
      b <- blocks[[1]]
      Rprofmem(NULL)
      sizes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
      maxima <- c(maxima, max(sizes, na.rm = TRUE))
      expect_lte(tail(maxima, 1), budget)
      expect_lte(as.numeric(object.size(b$values)), budget)
      unlink(c(path, paste0(backing, c(".bk", ".rds"))))
    }
    expect_lte(maxima[2], maxima[1]) # more list bookkeeping can shrink blocks
  }
})

test_that("time blocks stay bounded as the number of volumes grows", {
  for (name in c("dense", "sparse", "big")) {
    maxima <- numeric()
    for (nt in c(8L, 256L)) {
      arr <- array(as.numeric(seq_len(128L * nt)), c(128L, 1L, 1L, nt))
      backing <- tempfile()
      x <- block_test_sources(arr, backing)[[name]]
      blocks <- vec_blocks(x, "time", budget = 65536)
      invisible(blocks[[1]])
      path <- tempfile()
      Rprofmem(path)
      b <- blocks[[1]]
      Rprofmem(NULL)
      sizes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
      maxima <- c(maxima, max(sizes, na.rm = TRUE))
      expect_lte(tail(maxima, 1), 65536)
      expect_equal(b$values, t(as.matrix(x)[b$voxels, b$times, drop = FALSE]))
      unlink(c(path, paste0(backing, c(".bk", ".rds"))))
    }
    expect_lte(maxima[2], maxima[1])
  }
})

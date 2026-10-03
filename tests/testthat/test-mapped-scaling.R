# Write the header and SHORT payload independently of neuroim2's writer.
mapped_short_fixture <- function(slope = 2.5, intercept = -4) {
  path <- tempfile(fileext = ".nii")
  dims <- c(3L, 2L, 2L, 3L)
  stored <- as.integer(seq_len(prod(dims)) - 18L)
  con <- file(path, "w+b")
  on.exit(close(con))
  writeBin(raw(352), con)
  put <- function(offset, value, size) {
    seek(con, offset, origin = "start", rw = "write")
    writeBin(value, con, size = size, endian = .Platform$endian)
  }
  put(0, 348L, 4)
  put(40, c(4L, dims, 1L, 1L, 1L), 2)
  put(70, 4L, 2) # NIfTI INT16
  put(72, 16L, 2)
  put(76, rep(1, 8), 4)
  put(108, c(352, slope, intercept), 4)
  seek(con, 344, rw = "write")
  writeBin(as.raw(c(110, 43, 49, 0)), con) # n+1\0
  put(352, stored, 2)
  list(path = path, dims = dims, stored = stored)
}

test_that("mapped public accessors decode scaled SHORT values exactly once", {
  for (pars in list(c(2.5, -4), c(-2, 9))) {
    fixture <- mapped_short_fixture(pars[1], pars[2])
    mapped <- MappedNeuroVec(fixture$path)
    expected <- array(fixture$stored * pars[1] + pars[2], fixture$dims)
    dense <- read_vec(fixture$path, mode = "normal")
    expect_equal(as.matrix(dense), matrix(expected, 12, 3))

    idx <- c(36, 1, 13, 1, 25, 12)
    expect_equal(linear_access(mapped, idx), as.numeric(expected)[idx])
    expect_equal(mapped[idx], as.numeric(expected)[idx])
    expect_equal(mapped[, , , c(3, 1, 3)], expected[, , , c(3, 1, 3)])
    expect_equal(series(mapped, c(12L, 1L, 12L)),
                 t(matrix(expected, 12, 3)[c(12, 1, 12), , drop = FALSE]))
    expect_equal(as.matrix(mapped), matrix(expected, 12, 3))
    expect_equal(as.matrix(sub_vector(mapped, c(3, 1, 3))),
                 matrix(expected, 12, 3)[, c(3, 1, 3)])
    expect_equal(split_reduce(mapped, factor(rep(1:3, each = 4))),
                 split_reduce(dense, factor(rep(1:3, each = 4))))
    mmap::munmap(mapped@filemap)
    unlink(fixture$path)
  }
})

test_that("mapped reads share scale normalization and per-volume selection", {
  fixture <- mapped_short_fixture(1, 0)
  on.exit(unlink(fixture$path))
  meta <- read_header(fixture$path)
  cases <- list(
    list(s = 1, b = 0, want_s = 1, want_b = 0),
    list(s = 0, b = 99, want_s = 1, want_b = 0),
    list(s = NA_real_, b = 99, want_s = 1, want_b = 0),
    list(s = NaN, b = 99, want_s = 1, want_b = 0),
    list(s = Inf, b = 99, want_s = 1, want_b = 0),
    list(s = numeric(), b = numeric(), want_s = 1, want_b = 0),
    list(s = 2, b = NaN, want_s = 2, want_b = 0),
    list(s = c(1, -2, 3), b = c(0, 10, -4),
         want_s = c(1, -2, 3), want_b = c(0, 10, -4))
  )
  idx <- c(36, 1, 13, 25, 24, 1)
  time <- (idx - 1) %/% 12 + 1
  for (case in cases) {
    meta@slope <- case$s
    meta@intercept <- case$b
    mapped <- load_data(new("MappedNeuroVecSource", meta_info = meta))
    expected <- fixture$stored[idx] * rep_len(case$want_s, 3)[time] +
      rep_len(case$want_b, 3)[time]
    expect_equal(linear_access(mapped, idx), expected)
    decoded <- matrix(fixture$stored, 12, 3)
    decoded <- sweep(decoded, 2, rep_len(case$want_s, 3), `*`)
    decoded <- sweep(decoded, 2, rep_len(case$want_b, 3), `+`)
    times <- c(3L, 1L, 3L, 2L)
    expect_equal(mapped[, , , times, drop = FALSE],
                 array(decoded[, times], c(3, 2, 2, length(times))))
    expect_equal(as.matrix(sub_vector(mapped, times)), decoded[, times])
    expect_length(linear_access(mapped, numeric()), 0)
    expect_error(linear_access(mapped, 0), "out of bounds")
    expect_error(linear_access(mapped, 37), "out of bounds")
    expect_error(linear_access(mapped, 1.5), "whole-number")
    expect_error(linear_access(mapped, Inf), "whole-number")
    mmap::munmap(mapped@filemap)
  }
})

make_filebacked_seq_fixture <- function() {
  fixture_dir <- tempfile("filebacked-neurovecseq-")
  dir.create(fixture_dir)

  dims <- c(3L, 4L, 2L, 2L)
  space <- NeuroSpace(dims)
  arrays <- list(
    array(seq_len(prod(dims)), dims),
    array(1000 + seq_len(prod(dims)), dims)
  )
  files <- file.path(fixture_dir, c("first.nii", "second.nii"))

  write_vec(NeuroVec(arrays[[1]], space), files[[1]])
  write_vec(NeuroVec(arrays[[2]], space), files[[2]])

  list(dir = fixture_dir, files = files, arrays = arrays)
}

test_that("multiple file-backed images form a lazy NeuroVecSeq", {
  fixture <- make_filebacked_seq_fixture()
  on.exit(unlink(fixture$dir, recursive = TRUE), add = TRUE)

  seq <- read_vec(fixture$files, mode = "filebacked")

  expect_s4_class(seq, "NeuroVecSeq")
  expect_equal(dim(seq), c(3L, 4L, 2L, 4L))
  expect_true(validObject(seq))
  expect_length(seq@vecs, 2L)
  expect_true(all(vapply(
    seq@vecs,
    inherits,
    logical(1),
    what = "FileBackedNeuroVec"
  )))
})

test_that("NeuroVecSeq validation does not read file-backed image data", {
  fixture <- make_filebacked_seq_fixture()
  on.exit(unlink(fixture$dir, recursive = TRUE), add = TRUE)

  parts <- lapply(fixture$files, FileBackedNeuroVec)
  unlink(fixture$files)

  expect_no_error(seq <- do.call(NeuroVecSeq, parts))
  expect_s4_class(seq, "NeuroVecSeq")
  expect_equal(dim(seq), c(3L, 4L, 2L, 4L))
  expect_true(validObject(seq))
})

test_that("NeuroVecSeq construction never invokes dense materialization", {
  fixture <- make_filebacked_seq_fixture()
  on.exit(unlink(fixture$dir, recursive = TRUE), add = TRUE)

  local_mocked_bindings(
    as.dense = function(x) stop("unexpected dense materialization"),
    .package = "neuroim2"
  )

  expect_no_error(seq <- read_vec(fixture$files, mode = "filebacked"))
  expect_s4_class(seq, "NeuroVecSeq")
  expect_true(all(vapply(
    seq@vecs,
    inherits,
    logical(1),
    what = "FileBackedNeuroVec"
  )))
})

test_that("FileBackedNeuroVec matrix conversion reads the expected values", {
  fixture <- make_filebacked_seq_fixture()
  on.exit(unlink(fixture$dir, recursive = TRUE), add = TRUE)

  part <- FileBackedNeuroVec(fixture$files[[1]])
  observed <- as.matrix(part)

  expect_equal(dim(observed), c(24L, 2L))
  expect_equal(observed, as(part, "matrix"))
  expect_equal(observed, matrix(fixture$arrays[[1]], nrow = 24L, ncol = 2L))
})

test_that("file-backed NeuroVecSeq materializes explicitly across segments", {
  fixture <- make_filebacked_seq_fixture()
  on.exit(unlink(fixture$dir, recursive = TRUE), add = TRUE)

  seq <- read_vec(fixture$files, mode = "filebacked")
  expected <- cbind(
    matrix(fixture$arrays[[1]], nrow = 24L, ncol = 2L),
    matrix(fixture$arrays[[2]], nrow = 24L, ncol = 2L)
  )

  expect_equal(as.matrix(seq), expected)
  expect_equal(series(seq, 1L), expected[1L, ])
  expect_equal(as.matrix(as.dense(seq)), expected)
})

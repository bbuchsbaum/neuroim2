test_that("NeuroVec exposes per-volume labels and supports character extraction", {
  vec <- make_vec(dim = c(2L, 2L, 2L), ntime = 3L)
  vec@volume_labels <- c("baseline", "task", "rest")

  expect_equal(volume_labels(vec), c("baseline", "task", "rest"))
  expect_equal(vec[["task"]]@.Data, vec[[2]]@.Data)

  sub <- sub_vector(vec, c("rest", "baseline"))
  expect_s4_class(sub, "DenseNeuroVec")
  expect_equal(dim(sub)[4], 2L)
  expect_equal(volume_labels(sub), c("rest", "baseline"))
  expect_equal(sub[[1]]@.Data, vec[[3]]@.Data)
  expect_equal(sub[[2]]@.Data, vec[[1]]@.Data)
})

test_that("concat preserves volume labels without silently rewriting duplicates", {
  x <- make_vec(dim = c(2L, 2L, 2L), ntime = 2L)
  y <- make_vec(dim = c(2L, 2L, 2L), ntime = 2L)
  x@volume_labels <- c("run1_a", "shared")
  y@volume_labels <- c("shared", "run2_b")

  out <- concat(x, y)

  expect_s4_class(out, "DenseNeuroVec")
  expect_equal(volume_labels(out), c("run1_a", "shared", "shared", "run2_b"))
  expect_error(out[["shared"]], "not unique")
})

test_that("concat uses empty strings for unlabeled inputs when labels are partially present", {
  x <- make_vec(dim = c(2L, 2L, 2L), ntime = 2L)
  y <- make_vec(dim = c(2L, 2L, 2L), ntime = 1L)
  x@volume_labels <- c("a", "b")

  out <- concat(x, y)
  expect_equal(volume_labels(out), c("a", "b", ""))

  z <- concat(make_vec(dim = c(2L, 2L, 2L), ntime = 1L), make_vec(dim = c(2L, 2L, 2L), ntime = 1L))
  expect_equal(volume_labels(z), character())
})

test_that("NeuroVecSeq concatenates labels across segments", {
  x <- make_vec(dim = c(2L, 2L, 2L), ntime = 2L)
  y <- make_vec(dim = c(2L, 2L, 2L), ntime = 2L)
  x@volume_labels <- c("x1", "x2")
  y@volume_labels <- c("y1", "y2")

  seq_vec <- NeuroVecSeq(x, y)

  expect_equal(volume_labels(seq_vec), c("x1", "x2", "y1", "y2"))
  expect_equal(seq_vec[["y1"]]@.Data, y[[1]]@.Data)

  sub <- sub_vector(seq_vec, c("x2", "y2"))
  expect_equal(volume_labels(sub), c("x2", "y2"))
})

test_that("volume labels round-trip through NIfTI and are available in low-footprint readers", {
  vec <- make_vec(dim = c(2L, 2L, 2L), ntime = 3L)
  vec@volume_labels <- c("baseline", "task A", "rest/2")

  tmp <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp), add = TRUE)

  write_vec(vec, tmp)

  hdr <- read_header(tmp)
  expect_equal(
    neuroim2:::nifti_volume_labels(hdr@header, expected_length = 3L),
    c("baseline", "task A", "rest/2")
  )

  dense <- read_vec(tmp)
  fb <- read_vec(tmp, mode = "filebacked")

  expect_equal(volume_labels(dense), c("baseline", "task A", "rest/2"))
  expect_equal(volume_labels(fb), c("baseline", "task A", "rest/2"))
  expect_equal(dense[["task A"]]@.Data, vec[[2]]@.Data, tolerance = 1e-6)
  expect_equal(fb[["rest/2"]]@.Data, vec[[3]]@.Data, tolerance = 1e-6)
})

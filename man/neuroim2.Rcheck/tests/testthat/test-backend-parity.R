library(testthat)
library(neuroim2)

context("backend parity")

backend_fixture <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")

load_backend_vecs <- function() {
  list(
    normal = read_vec(backend_fixture, mode = "normal"),
    mmap = read_vec(backend_fixture, mode = "mmap"),
    filebacked = read_vec(backend_fixture, mode = "filebacked")
  )
}

test_that("read_vec backends expose consistent dimensions and classes", {
  vecs <- load_backend_vecs()

  expect_s4_class(vecs$normal, "DenseNeuroVec")
  expect_s4_class(vecs$mmap, "MappedNeuroVec")
  expect_s4_class(vecs$filebacked, "FileBackedNeuroVec")

  expect_equal(dim(vecs$normal), dim(vecs$mmap))
  expect_equal(dim(vecs$normal), dim(vecs$filebacked))
  expect_equal(spacing(vecs$normal), spacing(vecs$mmap))
  expect_equal(spacing(vecs$normal), spacing(vecs$filebacked))
  expect_equal(origin(vecs$normal), origin(vecs$mmap))
  expect_equal(origin(vecs$normal), origin(vecs$filebacked))
  expect_equal(trans(vecs$normal), trans(vecs$mmap))
  expect_equal(trans(vecs$normal), trans(vecs$filebacked))
})

test_that("read_vec backends agree on direct indexing", {
  vecs <- load_backend_vecs()
  ref <- vecs$normal
  ind <- c(1L, 25L, 100L, 250L)
  cds <- index_to_grid(space(ref), ind)

  expect_equal(vecs$mmap[1, 1, 1, ], ref[1, 1, 1, ])
  expect_equal(vecs$filebacked[1, 1, 1, ], ref[1, 1, 1, ])

  expect_equal(vecs$mmap[1:2, 2, 3, 3:4], ref[1:2, 2, 3, 3:4])
  expect_equal(vecs$filebacked[1:2, 2, 3, 3:4], ref[1:2, 2, 3, 3:4])

  expect_equal(vecs$mmap[ind], ref[ind])
  expect_equal(vecs$filebacked[ind], ref[ind])
  expect_equal(vecs$mmap[cds], ref[cds])
  expect_equal(vecs$filebacked[cds], ref[cds])
})

test_that("read_vec backends agree on series extraction", {
  vecs <- load_backend_vecs()
  ref <- vecs$normal
  voxel_coords <- rbind(
    c(1, 1, 1),
    c(2, 2, 2),
    c(3, 3, 3)
  )
  linear_idx <- c(1L, 10L, 100L)

  expect_equal(series(vecs$mmap, 1, 1, 1), series(ref, 1, 1, 1))
  expect_equal(series(vecs$filebacked, 1, 1, 1), series(ref, 1, 1, 1))

  expect_equal(series(vecs$mmap, voxel_coords), series(ref, voxel_coords))
  expect_equal(series(vecs$filebacked, voxel_coords), series(ref, voxel_coords))

  expect_equal(series(vecs$mmap, linear_idx), series(ref, linear_idx))
  expect_equal(series(vecs$filebacked, linear_idx), series(ref, linear_idx))
})

test_that("read_vec backends agree on sub_vector and matrix conversion", {
  vecs <- load_backend_vecs()
  ref <- vecs$normal

  expect_equal(as.matrix(sub_vector(vecs$mmap, 1:3)), as.matrix(sub_vector(ref, 1:3)))
  expect_equal(as.matrix(sub_vector(vecs$filebacked, 1:3)), as.matrix(sub_vector(ref, 1:3)))

  expect_equal(as(vecs$mmap, "matrix"), as(ref, "matrix"))
  expect_equal(as(vecs$filebacked, "matrix"), as(ref, "matrix"))
})

test_that("read_vec backends agree on ROI extraction paths", {
  vecs <- load_backend_vecs()
  ref <- vecs$normal
  roi_idx <- 1:10
  roi_ref <- series_roi(ref, roi_idx)
  roi_vol <- ROIVol(space(drop_dim(space(ref))), coords(roi_ref), data = rep(1, nrow(coords(roi_ref))))

  expect_equal(values(series_roi(vecs$mmap, roi_idx)), values(roi_ref))
  expect_equal(values(series_roi(vecs$filebacked, roi_idx)), values(roi_ref))

  expect_equal(values(series_roi(vecs$mmap, roi_vol)), values(series_roi(ref, roi_vol)))
  expect_equal(values(series_roi(vecs$filebacked, roi_vol)), values(series_roi(ref, roi_vol)))
})

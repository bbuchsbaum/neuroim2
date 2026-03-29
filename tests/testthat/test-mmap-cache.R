context("mmap cache")

library(neuroim2)

test_that("memory-mapped reads reuse a cached handle for the same file", {
  ns <- asNamespace("neuroim2")
  clear_cache <- get(".clear_mmap_cache", envir = ns)
  cache_size <- get(".mmap_cache_size", envir = ns)
  read_data <- get("read_mapped_data", envir = ns)
  read_series <- get("read_mapped_series", envir = ns)
  read_vols <- get("read_mapped_vols", envir = ns)

  clear_cache()
  on.exit(clear_cache(), add = TRUE)

  fname <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  meta1 <- FileBackedNeuroVec(fname)@meta
  meta2 <- FileBackedNeuroVec(fname)@meta

  expect_equal(cache_size(), 0L)

  expect_length(read_data(meta1, c(1L, 2L, 3L)), 3L)
  expect_equal(cache_size(), 1L)

  expect_equal(dim(read_series(meta1, c(1L, 2L))), c(meta1@dims[4], 2L))
  expect_equal(cache_size(), 1L)

  expect_equal(dim(read_vols(meta2, 1:2)), c(2L, prod(meta2@dims[1:3])))
  expect_equal(cache_size(), 1L)

  clear_cache()
  expect_equal(cache_size(), 0L)
})

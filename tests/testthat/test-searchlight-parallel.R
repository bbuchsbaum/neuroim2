context("searchlight parallel")

library(neuroim2)

make_mask_vol <- function() {
  sp <- NeuroSpace(c(5, 5, 5))
  arr <- array(1, dim = c(5, 5, 5))
  NeuroVol(arr, sp)
}

test_that("searchlight eager parallel matches sequential results", {
  skip_on_cran()

  mask <- make_mask_vol()
  seq_res <- searchlight(mask, radius = 1, eager = TRUE, nonzero = FALSE, cores = 0)
  par_res <- searchlight(mask, radius = 1, eager = TRUE, nonzero = FALSE, cores = 2)

  expect_equal(length(par_res), length(seq_res))
  expect_equal(par_res[1:10], seq_res[1:10])
  expect_equal(vapply(par_res[1:10], attr, integer(1), which = "mask_index"),
               seq_len(10))
})

test_that("searchlight_coords parallel matches sequential neighborhoods", {
  skip_on_cran()

  mask <- make_mask_vol()
  seq_res <- searchlight_coords(mask, radius = 1, nonzero = FALSE, cores = 0)
  par_res <- searchlight_coords(mask, radius = 1, nonzero = FALSE, cores = 2)

  expect_equal(length(par_res), length(seq_res))
  expect_equal(par_res[1:10], seq_res[1:10])
})

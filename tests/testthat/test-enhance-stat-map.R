library(testthat)
library(neuroim2)

# Build a deterministic noisy statistical map:
#   - a coherent signal blob in the centre,
#   - a Gaussian noise floor,
#   - sparse high-amplitude "salt-and-pepper" impulses.
make_noisy_stat <- function(dims = c(20L, 20L, 20L), seed = 42L,
                            blob_amp = 6, spike_amp = 10, n_spikes = 40L) {
  set.seed(seed)
  sp <- NeuroSpace(dims)
  arr <- array(stats::rnorm(prod(dims)), dims)

  # signal blob (centre 6x6x6)
  cn <- dims %/% 2L
  bx <- (cn[1] - 2L):(cn[1] + 3L)
  by <- (cn[2] - 2L):(cn[2] + 3L)
  bz <- (cn[3] - 2L):(cn[3] + 3L)
  arr[bx, by, bz] <- arr[bx, by, bz] + blob_amp

  # salt-and-pepper impulses, kept away from the blob core
  all_idx <- seq_len(prod(dims))
  blob_lin <- as.vector(outer(bx, (by - 1L) * dims[1], "+"))
  blob_lin <- as.vector(outer(blob_lin, (bz - 1L) * dims[1] * dims[2], "+"))
  cand <- setdiff(all_idx, blob_lin)
  spk <- cand[seq(1L, length(cand), length.out = n_spikes)]
  arr[spk] <- spike_amp * rep_len(c(1, -1), n_spikes)

  list(vol = NeuroVol(arr, sp), space = sp, dims = dims,
       blob = list(x = bx, y = by, z = bz), spikes = spk)
}

test_that("enhance_stat_map returns a same-grid NeuroVol", {
  d <- make_noisy_stat()
  out <- enhance_stat_map(d$vol)
  expect_s4_class(out, "NeuroVol")
  expect_equal(dim(out), d$dims)
  expect_equal(space(out), d$space)
  expect_false(any(is.na(as.array(out))))
})

test_that("salt-and-pepper impulses are suppressed", {
  d <- make_noisy_stat()
  out <- enhance_stat_map(d$vol)
  orig_spike <- mean(abs(d$vol[d$spikes]))
  new_spike <- mean(abs(out[d$spikes]))
  # impulses should be knocked well down toward their local (near-zero) median
  expect_lt(new_spike, 0.5 * orig_spike)
})

test_that("signal peaks are preserved, not depressed", {
  d <- make_noisy_stat()
  core <- list(d$blob$x, d$blob$y, d$blob$z)
  orig_core <- mean(d$vol[core[[1]], core[[2]], core[[3]]])
  out <- enhance_stat_map(d$vol)
  new_core <- mean(out[core[[1]], core[[2]], core[[3]]])
  # the blob carries real structure: its amplitude must not be depressed
  expect_gt(new_core, 0.8 * orig_core)
})

test_that("the noise floor is denoised", {
  d <- make_noisy_stat()
  # background voxels: everything except the blob and the spikes
  is_blob <- array(FALSE, d$dims)
  is_blob[d$blob$x, d$blob$y, d$blob$z] <- TRUE
  bg <- setdiff(which(!is_blob), d$spikes)
  out <- enhance_stat_map(d$vol)
  expect_lt(stats::sd(out[bg]), stats::sd(d$vol[bg]))
})

test_that("output is restricted to the mask", {
  d <- make_noisy_stat()
  arr_mask <- array(FALSE, d$dims)
  arr_mask[1:10, , ] <- TRUE
  mask <- LogicalNeuroVol(arr_mask, d$space)
  out <- enhance_stat_map(d$vol, mask)
  out_arr <- as.array(out)
  expect_true(all(out_arr[11:20, , ] == 0))
  expect_true(any(out_arr[1:10, , ] != 0))
})

test_that("all base methods run and return NeuroVol", {
  d <- make_noisy_stat()
  for (m in c("guided", "bilateral", "gaussian")) {
    out <- enhance_stat_map(d$vol, method = m)
    expect_s4_class(out, "NeuroVol")
    expect_equal(dim(out), d$dims)
    expect_false(any(is.na(as.array(out))), info = m)
  }
})

test_that("detail_gain = 0 gives a pure denoise (still smooths, no NaN)", {
  d <- make_noisy_stat()
  out <- enhance_stat_map(d$vol, detail_gain = 0)
  expect_s4_class(out, "NeuroVol")
  expect_false(any(is.na(as.array(out))))
  idx <- which(d$vol != 0)
  expect_lt(stats::sd(out[idx]), stats::sd(d$vol[idx]))
})

test_that("despike = FALSE leaves impulses largely intact", {
  d <- make_noisy_stat()
  with_despike <- enhance_stat_map(d$vol, despike = TRUE)
  no_despike <- enhance_stat_map(d$vol, despike = FALSE)
  expect_gt(mean(abs(no_despike[d$spikes])), mean(abs(with_despike[d$spikes])))
})

test_that("enhance_stat_map validates its arguments", {
  d <- make_noisy_stat()
  expect_error(enhance_stat_map(42), "NeuroVol")
  expect_error(enhance_stat_map(d$vol, despike_k = -1), "despike_k")
  expect_error(enhance_stat_map(d$vol, radius = 0), "radius")
  expect_error(enhance_stat_map(d$vol, detail_gain = -1), "detail_gain")
  expect_error(enhance_stat_map(d$vol, epsilon = -2), "epsilon")
})

test_that("plot_overlay and plot_ortho accept the enhance argument", {
  d <- make_noisy_stat()
  bg <- NeuroVol(array(1, d$dims), d$space)

  p_true <- plot_overlay(bg, d$vol, enhance = TRUE, draw = FALSE)
  expect_false(is.null(p_true))

  p_list <- plot_overlay(bg, d$vol,
                         enhance = list(detail_gain = 2, method = "bilateral"),
                         draw = FALSE)
  expect_false(is.null(p_list))

  p_off <- plot_overlay(bg, d$vol, enhance = FALSE, draw = FALSE)
  expect_false(is.null(p_off))

  p_ortho <- plot_ortho(d$vol, enhance = TRUE, draw = FALSE)
  expect_false(is.null(p_ortho))

  expect_error(plot_overlay(bg, d$vol, enhance = "yes", draw = FALSE),
               "enhance")
})

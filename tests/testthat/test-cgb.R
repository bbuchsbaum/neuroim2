test_that("cgb graph builds and smooths basic runs", {
  skip_on_cran()
  set.seed(123)
  nx <- 8; ny <- 6; nz <- 4; nt1 <- 60; nt2 <- 45
  sp1 <- NeuroSpace(c(nx, ny, nz, nt1), spacing = c(2, 2, 2))
  sp2 <- NeuroSpace(c(nx, ny, nz, nt2), spacing = c(2, 2, 2))

  arr1 <- array(rnorm(nx * ny * nz * nt1, sd = 0.5), c(nx, ny, nz, nt1))
  arr2 <- array(rnorm(nx * ny * nz * nt2, sd = 0.5), c(nx, ny, nz, nt2))

  pat1 <- base::scale(rnorm(nt1))
  pat2 <- base::scale(rnorm(nt2))
  arr1[3:4, 3:4, 2:3, ] <- arr1[3:4, 3:4, 2:3, ] + rep(pat1, each = 8)
  arr2[3:4, 3:4, 2:3, ] <- arr2[3:4, 3:4, 2:3, ] + rep(pat2, each = 8)

  run1 <- DenseNeuroVec(arr1, sp1)
  run2 <- DenseNeuroVec(arr2, sp2)
  mask <- LogicalNeuroVol(array(TRUE, dim(drop_dim(sp1))), drop_dim(sp1))

  g <- cgb_make_graph(list(run1, run2),
                      mask = mask,
                      window = 1L,
                      spatial_sigma = 2,
                      corr_map = "power",
                      corr_param = 2,
                      topk = 10)

  expect_true(is.list(g))
  expect_equal(length(g$row_ptr), length(g$mask_idx) + 1L)
  expect_true(length(g$val) >= length(g$mask_idx))

  sample_rows <- sample.int(length(g$mask_idx), min(25L, length(g$mask_idx)))
  row_sums <- vapply(sample_rows, function(i) {
    rng <- (g$row_ptr[i] + 1L):g$row_ptr[i + 1L]
    sum(g$val[rng])
  }, numeric(1))
  expect_true(all(abs(row_sums - 1) < 1e-8))

  smoothed <- cgb_smooth(run1, g, passes = 1L, lambda = 1)
  expect_s4_class(smoothed, "DenseNeuroVec")
  expect_equal(dim(smoothed), dim(run1))
})

test_that("cgb graph supports time weights, confounds, and LORO", {
  skip_on_cran()
  set.seed(456)
  nx <- 6; ny <- 5; nz <- 4; nt1 <- 40; nt2 <- 32
  sp1 <- NeuroSpace(c(nx, ny, nz, nt1), spacing = c(2, 2, 2))
  sp2 <- NeuroSpace(c(nx, ny, nz, nt2), spacing = c(2, 2, 2))

  arr1 <- array(rnorm(nx * ny * nz * nt1), c(nx, ny, nz, nt1))
  arr2 <- array(rnorm(nx * ny * nz * nt2), c(nx, ny, nz, nt2))
  run1 <- DenseNeuroVec(arr1, sp1)
  run2 <- DenseNeuroVec(arr2, sp2)
  mask <- LogicalNeuroVol(array(TRUE, dim(drop_dim(sp1))), drop_dim(sp1))

  w1 <- make_time_weights(fd = abs(rnorm(nt1)))
  w2 <- make_time_weights(fd = abs(rnorm(nt2)))
  conf1 <- cbind(seq_len(nt1), as.numeric(base::scale(seq_len(nt1))))
  conf2 <- cbind(seq_len(nt2), as.numeric(base::scale(seq_len(nt2))))

  graphs <- cgb_make_graph(list(run1, run2),
                           mask = mask,
                           window = 1L,
                           spatial_sigma = 1.5,
                           corr_map = "power",
                           corr_param = 2,
                           time_weights = list(w1, w2),
                           confounds = list(conf1, conf2),
                           leave_one_out = TRUE,
                           robust = "huber",
                           robust_c = 1.5)

  expect_length(graphs, 2L)
  lapply(graphs, function(g) {
    expect_equal(length(g$row_ptr), length(g$mask_idx) + 1L)
    sample_rows <- sample.int(length(g$mask_idx), min(10L, length(g$mask_idx)))
    sums <- vapply(sample_rows, function(i) {
      rng <- (g$row_ptr[i] + 1L):g$row_ptr[i + 1L]
      sum(g$val[rng])
    }, numeric(1))
    expect_true(all(abs(sums - 1) < 1e-8))
    invisible()
  })

loro <- cgb_smooth_loro(list(run1, run2), graphs, passes = 1L, lambda = 0.8)
  expect_length(loro, 2L)
  expect_s4_class(loro[[1]], "DenseNeuroVec")
  expect_equal(dim(loro[[2]]), dim(run2))
})

test_that("cgb smoothing preserves voxels outside the mask", {
  skip_on_cran()
  set.seed(99)
  nx <- 5; ny <- 4; nz <- 3; nt <- 20
  sp <- NeuroSpace(c(nx, ny, nz, nt), spacing = c(2, 2, 2))
  arr <- array(rnorm(nx * ny * nz * nt), c(nx, ny, nz, nt))
  # Create a mask excluding a corner block
  mask_arr <- array(TRUE, dim(drop_dim(sp)))
  mask_arr[1:2, 1:2, 1:2] <- FALSE
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp))
  run <- DenseNeuroVec(arr, sp)

  g <- cgb_make_graph(run,
                      mask = mask,
                      window = 1L,
                      spatial_sigma = 1.5,
                      corr_map = "power",
                      corr_param = 2,
                      topk = 8)
  smoothed <- cgb_smooth(run, g, passes = 1L, lambda = 1)
  arr_out <- as.array(smoothed)
  arr_in <- as.array(run)
  expect_equal(arr_out[!mask_arr], arr_in[!mask_arr])
})

test_that("self-edge is retained even with tiny top-k", {
  skip_on_cran()
  set.seed(42)
  nx <- 4; ny <- 4; nz <- 3; nt <- 15
  sp <- NeuroSpace(c(nx, ny, nz, nt), spacing = c(2, 2, 2))
  run <- DenseNeuroVec(array(rnorm(nx * ny * nz * nt), c(nx, ny, nz, nt)), sp)
  mask <- LogicalNeuroVol(array(TRUE, c(nx, ny, nz)), drop_dim(sp))

  g <- cgb_make_graph(run,
                      mask = mask,
                      window = 1L,
                      spatial_sigma = 2,
                      corr_map = "power",
                      corr_param = 2,
                      topk = 1,
                      add_self = TRUE)

  center_idx <- g$mask_idx
  rp <- g$row_ptr
  for (i in seq_along(center_idx)) {
    cols <- g$col_ind[rp[i]:(rp[i + 1] - 1L)]
    expect_true((center_idx[i] - 1L) %in% cols)
  }
})

test_that("intercept-only confounds reproduce baseline smoothing", {
  skip_on_cran()
  set.seed(321)
  nx <- 4; ny <- 4; nz <- 3; nt <- 25
  sp <- NeuroSpace(c(nx, ny, nz, nt), spacing = c(2, 2, 2))
  run <- DenseNeuroVec(array(rnorm(nx * ny * nz * nt), c(nx, ny, nz, nt)), sp)
  mask <- LogicalNeuroVol(array(TRUE, c(nx, ny, nz)), drop_dim(sp))

  g_plain <- cgb_make_graph(run,
                            mask = mask,
                            window = 1L,
                            spatial_sigma = 1.5,
                            corr_map = "power",
                            corr_param = 2,
                            topk = 6,
                            add_self = TRUE)

  intercepts <- list(matrix(0, nrow = nt, ncol = 0))
  g_proj <- cgb_make_graph(run,
                           mask = mask,
                           window = 1L,
                           spatial_sigma = 1.5,
                           corr_map = "power",
                           corr_param = 2,
                           topk = 6,
                           add_self = TRUE,
                           confounds = intercepts)
  # Compare smoothing outputs rather than CSR internals (ties/top-k ordering
  # may lead to tiny structural differences with equivalent effect)
  s_plain <- cgb_smooth(run, g_plain, passes = 1L, lambda = 1)
  s_proj  <- cgb_smooth(run, g_proj,  passes = 1L, lambda = 1)
  diff <- as.array(s_plain) - as.array(s_proj)
  expect_lt(max(abs(diff)), 1e-6)
})

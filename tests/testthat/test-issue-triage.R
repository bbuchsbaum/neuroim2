test_that("explicit sparse orientation preserves non-symmetric square series", {
  mask <- array(seq_len(8) %in% c(1, 3, 6), c(2, 2, 2))
  ix <- c(1L, 3L, 6L)
  for (nt in c(3L, 4L)) {
    sp <- NeuroSpace(c(2, 2, 2, nt))
    y <- matrix(seq_len(nt * 3), nt, 3)
    a <- SparseNeuroVec(y, sp, mask, orientation = "time_x_voxels")
    b <- SparseNeuroVec(t(y), sp, mask, orientation = "voxels_x_time")
    expect_equal(series(a, ix), y)
    expect_equal(series(b, ix), y)
    full <- matrix(0, 8, nt); full[ix, ] <- t(y)
    expect_equal(as.matrix(a), full)
    expect_equal(as.matrix(b), full)
    # Preserve the existing auto convention for square and rectangular input.
    expect_equal(series(SparseNeuroVec(t(y), sp, mask), ix), y)
  }
  expect_error(SparseNeuroVec(matrix(1, 4, 3), NeuroSpace(c(2,2,2,4)), mask,
                             orientation = "voxels_x_time"), "orientation")
  expect_error(SparseNeuroVec(matrix(1, 4, 3), NeuroSpace(c(2,2,2,4)), mask,
                             orientation = "bad"), "arg")
})

test_that("default background limits retain bright tissue; robust stays opt-in", {
  a <- array(0, c(20, 20, 1)); a[5:15,5:15,1] <- 100; a[9:10,9:10,1] <- 200
  s <- NeuroSpace(dim(a)); bg <- NeuroVol(a,s); ov <- NeuroVol(array(0,dim(a)),s)
  make <- function(...) plot_overlay(bg, ov, zlevels=1, draw=FALSE, assemble=FALSE, ...)
  expect_equal(make()[[1]]$scales$get_scales("fill")$limits, c(0,200))
  expect_equal(make(bg_range="robust")[[1]]$scales$get_scales("fill")$limits, c(0,100))
  expect_equal(make(bg_range=c(10,150))[[1]]$scales$get_scales("fill")$limits, c(10,150))
})

test_that("soft alpha controls affect rendered pixels independently of color", {
  dims <- c(3, 2, 1); s <- NeuroSpace(dims)
  bg <- NeuroVol(array(1, dims),s)
  vals <- c(0, 1, 2, 3, 8, NA_real_)
  ov <- NeuroVol(array(vals, dims),s)
  make <- function(...) plot_overlay(bg, ov, zlevels=1, draw=FALSE,
    ov_alpha_mode="soft", ov_thresh=1, ov_alpha=0.8,
    alpha_knee=1, alpha_cap=3, alpha_floor=0.25, alpha_gamma=1, ...)
  p <- make(ov_range=c(0,8), assemble=FALSE)
  pars <- attr(p,"soft_alpha")
  expect_equal(pars, list(lo=1,hi=3,gamma=1,alpha_floor=0.25))
  # Check actual raster bytes, independently of the curve helper.
  raster <- p[[1]]$layers[[2]]$geom_params$grob$raster
  alpha <- grDevices::col2rgb(as.vector(raster), alpha=TRUE)[4,]/255
  expect_equal(sort(unname(alpha)), sort(c(0,0.2,0.5,0.8,0.8,0)), tolerance=1/255)
  expect_identical(attr(make(ov_range=c(0,20),assemble=FALSE),"soft_alpha"),pars)
  expect_identical(attr(make(ov_range=c(0,8),assemble=TRUE),"soft_alpha"),pars)
  # Fixed anchors and exponent are independent of the supplied distribution.
  expect_identical(soft_alpha_params(1:100,knee=1,cap=3,gamma=1,alpha_floor=.25),pars)
})

test_that("soft alpha validates controls and supports transparent empty input", {
  for (args in list(list(gamma=0),list(knee=-1),list(cap=Inf),list(alpha_floor=2),
                    list(alpha_mid=1),list(gamma_min=2,gamma_max=1),list(knee=3,cap=2))) {
    expect_error(do.call(soft_alpha_params,c(list(mags=1:10),args)))
  }
  p <- soft_alpha_params(c(NA,Inf,0))
  expect_true(all(is.finite(unlist(p))))
  expect_equal(p$lo,0)
  expect_equal(soft_alpha_params(1:30,knee=0,cap=30,gamma_min=.5,gamma_max=.5)$gamma,.5)
})

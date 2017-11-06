

library(neuroim2)
library(testthat)

test_that("can compute connected components for random masks", { 
  
  D <- 64
  for (i in 1:10) {
    mat <- array(sample(c(0, 1), size = D^3, replace = TRUE), dim = rep(D, 3))
    nspace <- NeuroSpace(c(D,D,D))
    vol <- NeuroVol(mat, nspace)
    mask <- as.logical(vol)
    cc <- conn_comp(mask)
  }
  
  TRUE
  
})

test_that("conn_comp finds the correct number of clusters", {
  arr <- array(0,rep(50,3))
  spc <- NeuroSpace(rep(50,3), spacing=c(1,1,1))
  vol <- NeuroVol(arr, spc)
  
  sphere1 <- spherical_roi(vol, c(25,25,25), 2, fill=1)
  vol[coords(sphere1)] <- 1
  
  sphere2 <- spherical_roi(vol, c(8,8,8), 1, fill=1)
  vol[coords(sphere2)] <- 1
  
  cc <- conn_comp(as.logical(vol))
  
  expect_equal(sum(cc$index==1), nrow(coords(sphere1)))
  expect_equal(sum(cc$index == 2), nrow(coords(sphere2)))
  expect_equal(max(cc$index), 2)
  expect_equal(max(cc$size), max(c(nrow(coords(sphere1)), nrow(coords(sphere2)))))
  
})

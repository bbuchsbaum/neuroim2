library(assertthat)
library(purrr)

context("connected components")

test_that("can compute connected components for random masks", {
  # Test with a mask that guarantees multiple components
  D <- 20
  mat <- array(0, dim=rep(D,3))
  
  # Create two distinct components
  mat[5:7, 5:7, 5:7] <- 1  # Component 1
  mat[15:17, 15:17, 15:17] <- 1  # Component 2
  
  nspace <- NeuroSpace(c(D,D,D))
  vol <- NeuroVol(mat, nspace)
  mask <- as.logical(vol)
  cc <- conn_comp(mask)
  
  # Check we have exactly 2 components
  expect_equal(max(cc$index@data), 2)
  
  # Now test random masks with controlled density
  for (i in 1:5) {
    # Use lower density to encourage multiple components
    prob <- runif(1, 0.1, 0.3)
    mat <- array(sample(c(0, 1), size = D^3, replace = TRUE, prob=c(1-prob, prob)), 
                dim = rep(D, 3))
    
    # Ensure at least two components by forcing disconnected regions
    mat[1:3, 1:3, 1:3] <- 1
    mat[D-2:D, D-2:D, D-2:D] <- 1
    
    nspace <- NeuroSpace(c(D,D,D))
    vol <- NeuroVol(mat, nspace)
    mask <- as.logical(vol)
    cc <- conn_comp(mask)
    
    # Should have at least 2 components
    expect_gte(max(cc$index@data), 2)
  }
})

test_that("conn_comp finds the correct number of clusters", {
  arr <- array(0, rep(50,3))
  spc <- NeuroSpace(rep(50,3), spacing=c(1,1,1))
  vol <- NeuroVol(arr, spc)
  
  # Create two spheres far apart to ensure separation
  sphere1 <- spherical_roi(vol, c(25,25,25), 2, fill=1)
  vol[coords(sphere1)] <- 1
  
  sphere2 <- spherical_roi(vol, c(8,8,8), 1, fill=1)
  vol[coords(sphere2)] <- 1
  
  cc <- conn_comp(as.logical(vol))
  
  expect_equal(sum(cc$index==1), nrow(coords(sphere1)))
  expect_equal(sum(cc$index==2), nrow(coords(sphere2)))
  expect_equal(max(cc$index), 2)
  expect_equal(max(cc$size), max(c(nrow(coords(sphere1)), nrow(coords(sphere2)))))
})

test_that("conn_comp handles edge cases correctly", {
  # Test empty mask
  arr_empty <- array(FALSE, c(10,10,10))
  cc_empty <- conn_comp_3D(arr_empty)
  expect_equal(max(cc_empty$index), 0)
  
  # Test single voxel
  arr_single <- array(FALSE, c(10,10,10))
  arr_single[5,5,5] <- TRUE
  cc_single <- conn_comp_3D(arr_single)
  expect_equal(max(cc_single$index), 1)
  expect_equal(max(cc_single$size), 1)
  
  # Test fully connected mask
  arr_full <- array(FALSE, c(5,5,5))
  arr_full[2:4, 2:4, 2:4] <- TRUE
  cc_full <- conn_comp_3D(arr_full)
  expect_equal(max(cc_full$index), 1)
  expect_equal(max(cc_full$size), sum(arr_full))
})

test_that("conn_comp respects different connectivity patterns", {
  # Create a mask that will have different components under different connectivity
  arr <- array(FALSE, c(5,5,5))
  arr[2,2,2] <- TRUE
  arr[3,3,3] <- TRUE
  
  # With 6-connectivity these should be separate
  cc6 <- conn_comp_3D(arr, connect="6-connect")
  expect_equal(max(cc6$index), 2)

  
  # With 26-connectivity these should be connected
  cc26 <- conn_comp_3D(arr, connect="26-connect")
  expect_equal(max(cc26$index), 1)
})

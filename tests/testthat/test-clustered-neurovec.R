# tests/testthat/test-clustered-neurovec.R

library(neuroim2)

test_that("ClusteredNeuroVec: construction, dims, and broadcast are correct", {
  # 3D space (2x2x2), mask covers all 8 voxels
  sp3   <- NeuroSpace(c(2,2,2), spacing = c(1,1,1))
  mask  <- LogicalNeuroVol(array(TRUE, dim = c(2,2,2)), sp3)

  # Two clusters: first 4 voxels -> 1, next 4 voxels -> 2
  cids  <- rep(c(1L, 2L), each = 4L)
  cvol  <- ClusteredNeuroVol(mask, cids)

  # 4D source data with T = 3 time points
  sp4   <- NeuroSpace(c(2,2,2,3), spacing = c(1,1,1))
  arr   <- array(NA_real_, dim = c(2,2,2,3))

  # Explicit voxel values by time (vectorized in R column-major order)
  v_t1 <- c(10,12,  8, 6,   1,3, 5,7)   # means: cl1=9, cl2=4
  v_t2 <- c(20,18, 22,20,   2,4, 6,8)   # means: cl1=20, cl2=5
  v_t3 <- c( 0, 0,  0,10,  10,10,10,10) # means: cl1=2.5, cl2=10
  arr[,,,1] <- array(v_t1, dim = c(2,2,2))
  arr[,,,2] <- array(v_t2, dim = c(2,2,2))
  arr[,,,3] <- array(v_t3, dim = c(2,2,2))

  vec  <- NeuroVec(arr, sp4)
  cv   <- ClusteredNeuroVec(vec, cvol, FUN = mean)

  # 4D dims preserved (3D from mask + time)
  expect_equal(dim(cv), c(2,2,2,3))

  # Broadcast check at each timepoint: constant within clusters
  vol1 <- cv[,,,1]; m1 <- mean(v_t1[1:4]); m2 <- mean(v_t1[5:8])
  vol2 <- cv[,,,2]; n1 <- mean(v_t2[1:4]); n2 <- mean(v_t2[5:8])
  vol3 <- cv[,,,3]; p1 <- mean(v_t3[1:4]); p2 <- mean(v_t3[5:8])

  a1 <- as.vector(vol1); a2 <- as.vector(vol2); a3 <- as.vector(vol3)
  expect_true(all(a1[which(cids==1L)] == m1) && all(a1[which(cids==2L)] == m2))
  expect_true(all(a2[which(cids==1L)] == n1) && all(a2[which(cids==2L)] == n2))
  expect_true(all(a3[which(cids==1L)] == p1) && all(a3[which(cids==2L)] == p2))
})

test_that("ClusteredNeuroVec: cluster-level time-series matrix is correct (T x K)", {
  sp3   <- NeuroSpace(c(2,2,2), spacing = c(1,1,1))
  mask  <- LogicalNeuroVol(array(TRUE, dim = c(2,2,2)), sp3)
  cids  <- rep(c(1L, 2L), each = 4L)
  cvol  <- ClusteredNeuroVol(mask, cids)

  sp4   <- NeuroSpace(c(2,2,2,3), spacing = c(1,1,1))
  arr   <- array(NA_real_, dim = c(2,2,2,3))
  v_t1 <- c(10,12,  8, 6,   1,3, 5,7)   # cl1=9, cl2=4
  v_t2 <- c(20,18, 22,20,   2,4, 6,8)   # cl1=20, cl2=5
  v_t3 <- c( 0, 0,  0,10,  10,10,10,10) # cl1=2.5, cl2=10
  arr[,,,1] <- array(v_t1, dim = c(2,2,2))
  arr[,,,2] <- array(v_t2, dim = c(2,2,2))
  arr[,,,3] <- array(v_t3, dim = c(2,2,2))
  vec  <- NeuroVec(arr, sp4)

  cv   <- ClusteredNeuroVec(vec, cvol, FUN = mean)
  # Expected T x K matrix
  expected <- cbind(
    c(mean(v_t1[1:4]), mean(v_t2[1:4]), mean(v_t3[1:4])),  # cluster 1
    c(mean(v_t1[5:8]), mean(v_t2[5:8]), mean(v_t3[5:8]))   # cluster 2
  )
  expect_equal(cv@ts, expected, tolerance = 1e-10)
})

test_that("ClusteredNeuroVec: voxels in the same cluster match across time", {
  sp3   <- NeuroSpace(c(2,2,2), spacing = c(1,1,1))
  mask  <- LogicalNeuroVol(array(TRUE, dim = c(2,2,2)), sp3)
  cids  <- rep(c(1L, 2L), each = 4L)
  cvol  <- ClusteredNeuroVol(mask, cids)

  sp4   <- NeuroSpace(c(2,2,2,3), spacing = c(1,1,1))
  arr   <- array(seq_len(2*2*2*3), dim = c(2,2,2,3))  # arbitrary but deterministic
  vec   <- NeuroVec(arr, sp4)
  cv    <- ClusteredNeuroVec(vec, cvol, FUN = mean)

  # Two coordinates known to fall in the first 4 linear voxels -> cluster 1
  for (tt in 1:3) {
    vol <- cv[,,,tt]
    expect_equal(vol[1,1,1], vol[2,2,1])  # both in cluster 1 => identical
  }
})

test_that("cluster_searchlight_series (k-NN): returns ROIs with correct shapes and seed-first", {
  # Make a 2x2x1 lattice with 4 single-voxel clusters at the corners.
  sp3   <- NeuroSpace(c(2,2,1), spacing = c(1,1,1))
  mask  <- LogicalNeuroVol(array(TRUE, dim = c(2,2,1)), sp3)
  cids  <- as.integer(1:4)
  cvol  <- ClusteredNeuroVol(mask, cids)

  # Build a vec whose voxel series are unique per cluster.
  sp4   <- NeuroSpace(c(2,2,1,2), spacing = c(1,1,1))
  arr   <- array(NA_real_, dim = c(2,2,1,2))
  arr[,,,1] <- array(c(10,11,12,13), dim = c(2,2,1))  # t1
  arr[,,,2] <- array(c(20,21,22,23), dim = c(2,2,1))  # t2
  vec   <- NeuroVec(arr, sp4)
  cv    <- ClusteredNeuroVec(vec, cvol, FUN = mean)

  wins  <- cluster_searchlight_series(cv, k = 2L)  # self + nearest neighbor
  expect_length(wins, 4)

  # Check first seed window
  roi1  <- wins[[1]]
  expect_true(inherits(roi1, "ROIVec"))
  V     <- values(roi1)  # (time) x (#neighbors)
  expect_equal(dim(V), c(2, 2))         # 2 timepoints, 2 neighbors
  expect_equal(V[,1], c(10, 20))        # seed (self) first column equals cluster-1 series
})

test_that("cluster_searchlight_series (radius): selects by centroid distance", {
  # Same 2x2x1 four-cluster setup as above
  sp3   <- NeuroSpace(c(2,2,1), spacing = c(1,1,1))
  mask  <- LogicalNeuroVol(array(TRUE, dim = c(2,2,1)), sp3)
  cids  <- as.integer(1:4)
  cvol  <- ClusteredNeuroVol(mask, cids)

  sp4   <- NeuroSpace(c(2,2,1,2), spacing = c(1,1,1))
  arr   <- array(NA_real_, dim = c(2,2,1,2))
  arr[,,,1] <- array(c(10,11,12,13), dim = c(2,2,1))
  arr[,,,2] <- array(c(20,21,22,23), dim = c(2,2,1))
  vec   <- NeuroVec(arr, sp4)
  cv    <- ClusteredNeuroVec(vec, cvol, FUN = mean)

  # Centroids are at (1,1,1),(2,1,1),(1,2,1),(2,2,1)
  # With radius=1.1 from seed (1,1,1), we include self + two orthogonal neighbors (distance 1)
  winsr <- cluster_searchlight_series(cv, radius = 1.1)
  roi1  <- winsr[[1]]
  V1    <- values(roi1)
  expect_true(inherits(roi1, "ROIVec"))
  expect_equal(nrow(V1), 2)     # time = 2
  expect_equal(ncol(V1), 3)     # self + two neighbors (exclude diagonal at ~1.414)
})
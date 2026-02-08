library(testthat)

# Test ClusteredNeuroVol function
test_that("ClusteredNeuroVol works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  expect_s4_class(clusvol, "ClusteredNeuroVol")
  # Add more tests to check the correctness of the output, e.g., dimensions, etc.
})

test_that("ClusteredNeuroVol supports as.array and as.vector", {
  sp <- NeuroSpace(c(2, 2, 2))
  mask_data <- array(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE), dim = c(2, 2, 2))
  mask <- LogicalNeuroVol(mask_data, sp)

  clusters <- as.integer(seq_len(sum(mask_data)))
  cvol <- ClusteredNeuroVol(mask, clusters)

  arr <- as.array(cvol)
  expect_equal(dim(arr), dim(sp))
  expect_equal(arr[mask_data], clusters)
  expect_true(all(arr[!mask_data] == 0))

  vec <- as.vector(cvol)
  expect_equal(vec, as.vector(arr))
})

# Test as.DenseNeuroVol function
test_that("as.DenseNeuroVol works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)
  # Use the previously created clusvol object from the ClusteredNeuroVol test
  dense_vol <- as(clusvol, "DenseNeuroVol")

  expect_s4_class(dense_vol, "DenseNeuroVol")
  # Add more tests to check the correctness of the output, e.g., dimensions, etc.
})

# Test show method for ClusteredNeuroVol
test_that("show method for ClusteredNeuroVol works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)
  # Capture the output of the show method
  output <- capture.output(show(clusvol))

  # Check if the output contains the necessary information
  expect_true(any(grepl("NeuroVol", output)))
  expect_true(any(grepl("Type", output)))
  # Add more tests to check the correctness of the output
})

# Test centroids function
test_that("centroids works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  centroids_com <- centroids(clusvol, type = "center_of_mass")
  centroids_medoid <- centroids(clusvol, type = "medoid")

  expect_equal(ncol(centroids_com), 3)
  expect_equal(nrow(centroids_com), num_clusters(clusvol))

  expect_equal(ncol(centroids_medoid), 3)
  expect_equal(nrow(centroids_medoid), num_clusters(clusvol))
})

# Test split_clusters function
test_that("split_clusters works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)
  vol <- NeuroVol(array(runif(16^3), c(16, 16, 16)), bspace)

  clusters_split <- split_clusters(vol, clusvol)

  expect_equal(length(clusters_split), num_clusters(clusvol))
  # Add more tests to check the correctness of the output
})

# Test num_clusters function
test_that("num_clusters works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  num_clus <- num_clusters(clusvol)

  expect_equal(num_clus, 10)
})

# Test as.dense function
test_that("as.dense works correctly", {
  skip_on_cran()
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  dense_vol <- as.dense(clusvol)

  expect_s4_class(dense_vol, "NeuroVol")
  # Add more tests to check the correctness of the output, e.g., dimensions, etc.
})

test_that("split_clusters handles non-contiguous cluster ids for NeuroVol and NeuroVec", {
  sp <- NeuroSpace(c(2,2,2,2))
  vals <- array(1:16, dim = c(2,2,2,2))
  vec <- NeuroVec(vals, sp)
  vol <- vec[[1]]  # first volume

  # mask all voxels; assign cluster ids 2,4,6
  mask <- LogicalNeuroVol(array(TRUE, dim=dim(vol)), space(vol))
  clusters_mask <- c(2,4,6,2,4,6,2,4)  # length 8
  cvol <- ClusteredNeuroVol(mask, clusters_mask)

  # NeuroVol + ClusteredNeuroVol path
  split_vol <- split_clusters(vol, cvol)
  expect_equal(length(split_vol), 3)
  means_vol <- sapply(split_vol, function(r) mean(values(r)))
  expect_equal(means_vol, c(mean(vol[cvol@cluster_map[["2"]]]),
                            mean(vol[cvol@cluster_map[["4"]]]),
                            mean(vol[cvol@cluster_map[["6"]]])))

  # NeuroVec + ClusteredNeuroVol path
  split_vec <- split_clusters(vec, cvol)
  expect_equal(length(split_vec), 3)
  # compare first volume timepoint to NeuroVol means
  means_vec <- sapply(split_vec, function(r) mean(values(r)[1,]))
  expect_equal(means_vec, means_vol)
})

test_that("split_clusters handles non-contiguous cluster ids for integer clusters", {
  sp <- NeuroSpace(c(2,2,2))
  vol <- NeuroVol(array(1:8, dim=c(2,2,2)), sp)
  clusters_full <- rep(0, 8)
  clusters_full[c(1,4)] <- 2
  clusters_full[c(2,6)] <- 4
  clusters_full[c(3,8)] <- 6

  # NeuroVol + integer
  split_int_vol <- split_clusters(vol, as.integer(clusters_full))
  expect_equal(length(split_int_vol), 3)
  expect_true(all(sapply(split_int_vol, function(r) length(values(r)) > 0)))

  # NeuroVec + integer
  vec <- NeuroVec(array(1:16, dim=c(2,2,2,2)), NeuroSpace(c(2,2,2,2)))
  split_int_vec <- split_clusters(vec, as.integer(clusters_full))
  expect_equal(length(split_int_vec), 3)
  expect_true(all(sapply(split_int_vec, function(r) length(values(r)) > 0)))
})

test_that("split_clusters for ClusteredNeuroVol (missing clusters arg) returns cluster ROIs", {
  sp <- NeuroSpace(c(2,2,2))
  mask <- LogicalNeuroVol(array(TRUE, dim=c(2,2,2)), sp)
  clusters_mask <- c(2,4,6,2,4,6,2,4)
  cvol <- ClusteredNeuroVol(mask, clusters_mask)

  rois <- split_clusters(cvol)
  expect_equal(length(rois), 3)
  expect_true(all(sapply(rois, function(r) inherits(r, "ROIVol"))))
  # each ROI should carry its cluster id as data
  expect_setequal(unlist(lapply(rois, values)), c(2,4,6))
})

# ---------------------------------------------------------------------------
# sub_clusters
# ---------------------------------------------------------------------------

test_that("sub_clusters by integer ID returns correct ClusteredNeuroVol", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(c(rep(TRUE, 500), rep(FALSE, 500)), c(10,10,10)), sp)
  clusters <- rep(1:5, length.out = 500)
  lm <- list(A = 1, B = 2, C = 3, D = 4, E = 5)
  cvol <- ClusteredNeuroVol(mask, clusters, label_map = lm)

  sub <- sub_clusters(cvol, c(1L, 3L))
  expect_s4_class(sub, "ClusteredNeuroVol")
  expect_equal(sort(unique(sub@clusters)), c(1L, 3L))
  expect_equal(num_clusters(sub), 2)
  expect_equal(sum(sub@mask), sum(clusters == 1) + sum(clusters == 3))
  expect_equal(sort(names(sub@label_map)), c("A", "C"))
})

test_that("sub_clusters by numeric coerces to integer", {
  sp <- NeuroSpace(c(6L, 6L, 6L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(6,6,6)), sp)
  clusters <- rep(1:3, length.out = 216)
  cvol <- ClusteredNeuroVol(mask, clusters)

  sub <- sub_clusters(cvol, c(2, 3))
  expect_equal(sort(unique(sub@clusters)), c(2L, 3L))
})

test_that("sub_clusters by name looks up label_map", {
  sp <- NeuroSpace(c(6L, 6L, 6L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(6,6,6)), sp)
  clusters <- rep(1:3, length.out = 216)
  lm <- list(FFA = 1, PPA = 2, V1 = 3)
  cvol <- ClusteredNeuroVol(mask, clusters, label_map = lm)

  sub <- sub_clusters(cvol, c("FFA", "V1"))
  expect_equal(sort(unique(sub@clusters)), c(1L, 3L))
  expect_equal(sort(names(sub@label_map)), c("FFA", "V1"))
})

test_that("sub_clusters errors on invalid IDs", {
  sp <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(4,4,4)), sp)
  cvol <- ClusteredNeuroVol(mask, rep(1:2, length.out = 64))

  expect_error(sub_clusters(cvol, 99L), "not found")
})

test_that("sub_clusters errors on invalid names", {
  sp <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(4,4,4)), sp)
  lm <- list(X = 1, Y = 2)
  cvol <- ClusteredNeuroVol(mask, rep(1:2, length.out = 64), label_map = lm)

  expect_error(sub_clusters(cvol, "ZZZ"), "not found")
})

test_that("sub_clusters works on ClusteredNeuroVec by integer", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(c(rep(TRUE, 500), rep(FALSE, 500)), c(10,10,10)), sp)
  clusters <- rep(1:5, length.out = 500)
  lm <- list(A = 1, B = 2, C = 3, D = 4, E = 5)
  cvol <- ClusteredNeuroVol(mask, clusters, label_map = lm)

  sp4 <- NeuroSpace(c(10L, 10L, 10L, 20L), c(1, 1, 1))
  vec <- NeuroVec(array(rnorm(10^3 * 20), c(10,10,10,20)), sp4)
  cvec <- ClusteredNeuroVec(vec, cvol)

  sub <- sub_clusters(cvec, c(2L, 4L))
  expect_s4_class(sub, "ClusteredNeuroVec")
  expect_equal(ncol(sub@ts), 2)
  expect_equal(nrow(sub@ts), 20)
  # cl_map should have contiguous column indices
  expect_equal(sort(unique(sub@cl_map[sub@cl_map > 0])), c(1L, 2L))
  # series extraction should work
  ts_val <- series(sub, 1, 1, 1)
  expect_equal(length(ts_val), 20)
})

test_that("sub_clusters works on ClusteredNeuroVec by name", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(c(rep(TRUE, 500), rep(FALSE, 500)), c(10,10,10)), sp)
  clusters <- rep(1:5, length.out = 500)
  lm <- list(A = 1, B = 2, C = 3, D = 4, E = 5)
  cvol <- ClusteredNeuroVol(mask, clusters, label_map = lm)

  sp4 <- NeuroSpace(c(10L, 10L, 10L, 20L), c(1, 1, 1))
  vec <- NeuroVec(array(rnorm(10^3 * 20), c(10,10,10,20)), sp4)
  cvec <- ClusteredNeuroVec(vec, cvol)

  sub <- sub_clusters(cvec, c("A", "C", "E"))
  expect_s4_class(sub, "ClusteredNeuroVec")
  expect_equal(ncol(sub@ts), 3)
  expect_equal(sort(names(sub@cvol@label_map)), c("A", "C", "E"))
})

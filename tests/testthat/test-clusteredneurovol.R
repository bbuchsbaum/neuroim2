library(testthat)

# Test ClusteredNeuroVol function
test_that("ClusteredNeuroVol works correctly", {
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  expect_s4_class(clusvol, "ClusteredNeuroVol")
  # Add more tests to check the correctness of the output, e.g., dimensions, etc.
})

# Test as.DenseNeuroVol function
test_that("as.DenseNeuroVol works correctly", {
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
  bspace <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
  grid <- index_to_grid(bspace, 1:(16 * 16 * 16))
  kres <- kmeans(grid, centers = 10)
  mask <- NeuroVol(rep(1, 16^3), bspace)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)

  dense_vol <- as.dense(clusvol)

  expect_s4_class(dense_vol, "NeuroVol")
  # Add more tests to check the correctness of the output, e.g., dimensions, etc.
})



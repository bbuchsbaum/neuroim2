# Independent values encode both voxel and time; singleton axes must survive
# conversion just as they do in the public voxel-by-time matrix contract.
for (nt in c(1L, 3L)) for (nv in c(1L, 3L)) {
  for (mask_type in c("numeric", "logical volume")) {
    test_that(paste("dense-to-sparse preserves", nt, "times and", nv,
                    "voxels with a", mask_type, "mask"), {
      dims <- c(3L, 3L, 1L, nt)
      sp <- NeuroSpace(dims, c(2, 3, 4), c(10, 20, 30))
      values <- outer(seq_len(9L), seq_len(nt), function(v, t) 100 * t + v)
      dense <- DenseNeuroVec(array(values, dims), sp)
      voxels <- c(1L, 4L, 8L)[seq_len(nv)]
      keep <- array(FALSE, dims[1:3])
      keep[voxels] <- TRUE
      selection <- if (mask_type == "numeric") as.numeric(voxels) else
        LogicalNeuroVol(keep, drop_dim(sp))

      expect_warning(sparse <- as.sparse(dense, selection), NA)
      expect_s4_class(sparse, "SparseNeuroVec")
      expect_equal(dim(sparse), dims)
      expect_equal(space(sparse), sp)
      expect_equal(indices(sparse), voxels)

      expected_matrix <- values
      expected_matrix[-voxels, ] <- 0
      expect_equal(as.matrix(sparse), expected_matrix)
      expected_series <- outer(seq_len(nt), voxels, function(t, v) 100 * t + v)
      expect_equal(series(sparse, voxels, drop = FALSE), expected_series)
    })
  }
}

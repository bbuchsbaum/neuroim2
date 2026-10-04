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

# Expected support is specified independently of mask construction. In
# particular, values must stay at their original voxels for unsorted masks.
mask_cases <- list(
  reversed_pair = list(mask = c(8, 1), support = c(1L, 8L)),
  reversed = list(mask = c(8, 4, 1), support = c(1L, 4L, 8L)),
  permuted = list(mask = c(4, 1, 8), support = c(1L, 4L, 8L)),
  duplicates = list(mask = c(8, 1, 8), support = c(1L, 8L)),
  repeated_singleton = list(mask = c(4, 4), support = 4L),
  with_zeros = list(mask = c(0, 8, 1, 0), support = c(1L, 8L)),
  negative = list(mask = -c(1, 3, 5), support = c(2L, 4L, 6L, 7L, 8L, 9L)),
  negative_duplicates = list(mask = c(-8, 0, -1, -8), support = c(2:7, 9L)),
  negative_all = list(mask = -(1:9), support = integer()),
  negative_outside = list(mask = -10, support = 1:9),
  empty = list(mask = numeric(), support = integer()),
  zero = list(mask = 0, support = integer()),
  all_reversed = list(mask = 9:1, support = 1:9),
  fractional_R_indices = list(mask = c(8.9, 1.2), support = c(1L, 8L))
)
for (nt in c(1L, 2L, 3L, 5L)) for (case in names(mask_cases)) {
  test_that(paste("numeric sparse mask keeps voxel identity:", case, "times", nt), {
    dims <- c(3L, 3L, 1L, nt)
    sp <- NeuroSpace(dims, c(2, 3, 4), c(10, 20, 30))
    values <- outer(seq_len(9L), seq_len(nt), function(v, t) 100 * t + v)
    dense <- DenseNeuroVec(array(values, dims), sp)
    selection <- mask_cases[[case]]
    expected <- matrix(0, 9L, nt)
    expected[selection$support, ] <- values[selection$support, , drop = FALSE]
    keep <- array(FALSE, dims[1:3])
    keep[selection$support] <- TRUE

    expect_warning(sparse <- as.sparse(dense, selection$mask), NA)
    expect_true(validObject(sparse))
    expect_equal(dim(sparse), dims)
    expect_equal(space(sparse), sp)
    expect_equal(indices(sparse), selection$support)
    expect_equal(as.matrix(sparse), expected)
    logical_sparse <- as.sparse(dense, LogicalNeuroVol(keep, drop_dim(sp)))
    expect_equal(as.matrix(sparse), as.matrix(logical_sparse))
    query <- c(8L, 1L, 8L, 4L)
    expect_equal(series(sparse, query, drop = FALSE), t(expected)[, query, drop = FALSE])
  })
}

test_that("numeric sparse masks reject unknown and out-of-bounds voxel indices", {
  sp <- NeuroSpace(c(3L, 3L, 1L, 1L))
  dense <- DenseNeuroVec(array(101:109, dim(sp)), sp)
  for (mask in list(NA_real_, NaN, Inf, -Inf, c(1, NA_real_), 10, c(1, 10))) {
    expect_error(as.sparse(dense, mask))
  }
  expect_error(as.sparse(dense, c(1, -2)))
})

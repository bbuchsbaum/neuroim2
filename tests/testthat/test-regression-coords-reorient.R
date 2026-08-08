library(testthat)
library(neuroim2)

# Regression tests for two coordinate defects fixed in 0.17.0.
#
# BUG 1 -- half-voxel offset.
#   index_to_coord() applied `index_to_grid(x, idx) - 0.5` and coord_to_index()
#   applied `grid + 0.5`, where grid_to_coord()/coord_to_grid() correctly use 1
#   (1-based grid to 0-based affine input). The shortcut conversions therefore
#   disagreed with the equivalent two-step path by exactly spacing/2, and
#   coord_to_index() could return a neighbouring voxel. On a 2mm grid,
#   coord_to_index(grid_to_coord(index_to_grid(12345))) returned 10264.
#
# BUG 2 -- reorient() ignored the source orientation and inverted the codes.
#   It built its transform from the target axis codes alone, so the operation
#   was absolute rather than relative and asking for an image's existing
#   orientation was not a no-op. It also passed the codes straight to
#   findAnatomy3D(), which names an axis by where it STARTS ("R" means
#   Right-to-Left) rather than where it increases towards, inverting every
#   request: reorient(ras_space, c("R","A","S")) returned LPI.
#
# The assertions below are written to fail loudly under either old behaviour.

# ---- fixtures ---------------------------------------------------------------

iso_space <- function() {
  NeuroSpace(c(8L, 9L, 7L), spacing = c(2, 2, 2), origin = c(-90, -126, -72))
}

aniso_space <- function() {
  NeuroSpace(c(7L, 6L, 5L), spacing = c(3, 1.5, 4), origin = c(10, -20, 30))
}

oblique_space <- function() {
  aff <- matrix(c(
    2.0, 0.3, 0.0, -90,
    0.0, 2.0, 0.2, -126,
    0.0, 0.0, 2.5, -72,
    0.0, 0.0, 0.0, 1
  ), nrow = 4, byrow = TRUE)
  NeuroSpace(dim = c(6L, 7L, 5L), trans = aff)
}

real_space <- function() {
  space(read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")))
}

# All 48 valid orientations: each image axis mapped to a distinct anatomical
# axis (3! orderings) with either sign (2^3).
all_orientations <- function() {
  pairs <- list(c("R", "L"), c("A", "P"), c("S", "I"))
  out <- list()
  for (perm in list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))) {
    for (s1 in 1:2) {
      for (s2 in 1:2) {
        for (s3 in 1:2) {
          out[[length(out) + 1]] <- c(
            pairs[[perm[1]]][s1], pairs[[perm[2]]][s2], pairs[[perm[3]]][s3]
          )
        }
      }
    }
  }
  out
}

# =============================================================================
# BUG 1: index_to_coord / coord_to_index half-voxel offset
# =============================================================================

test_that("index_to_coord matches the grid path for every voxel", {
  # Exhaustive rather than sampled: the old bug was a uniform offset, but a
  # partial fix could leave particular axes or edges wrong.
  for (sp in list(iso_space(), aniso_space(), oblique_space())) {
    idx <- seq_len(prod(dim(sp)))
    via_shortcut <- index_to_coord(sp, as.integer(idx))
    via_grid <- grid_to_coord(sp, index_to_grid(sp, as.integer(idx)))
    expect_equal(via_shortcut, via_grid)
  }
})

test_that("coord_to_index inverts index_to_coord for every voxel", {
  for (sp in list(iso_space(), aniso_space(), oblique_space())) {
    idx <- seq_len(prod(dim(sp)))
    expect_equal(
      as.vector(coord_to_index(sp, index_to_coord(sp, as.integer(idx)))),
      idx
    )
  }
})

test_that("coord_to_index matches grid_to_index on the same coordinates", {
  for (sp in list(iso_space(), aniso_space(), real_space())) {
    idx <- as.integer(seq(1, prod(dim(sp)), length.out = 200))
    world <- index_to_coord(sp, idx)
    expect_equal(
      as.vector(coord_to_index(sp, world)),
      as.vector(grid_to_index(sp, coord_to_grid(sp, world)))
    )
  }
})

test_that("voxel (1,1,1) lands exactly on the origin by either route", {
  for (sp in list(iso_space(), aniso_space(), real_space())) {
    expect_equal(as.vector(index_to_coord(sp, 1L)), origin(sp))
    expect_equal(as.vector(grid_to_coord(sp, matrix(c(1, 1, 1), nrow = 1))), origin(sp))
  }
})

test_that("the offset is a whole voxel, not half of one", {
  # Directly pins the defect: under the old code each coordinate was displaced
  # by spacing/2, so this difference was non-zero on every axis.
  sp <- aniso_space()
  idx <- as.integer(c(1, 5, 40, prod(dim(sp))))

  delta <- index_to_coord(sp, idx) - grid_to_coord(sp, index_to_grid(sp, idx))
  expect_true(all(abs(delta) < 1e-9))

  # And the value the old code would have produced is now demonstrably wrong.
  old_style <- grid_to_coord(sp, index_to_grid(sp, idx) + 0.5)
  expect_false(isTRUE(all.equal(index_to_coord(sp, idx), old_style)))
})

test_that("all index_to_coord dispatch paths agree", {
  sp <- aniso_space()
  arr <- array(0, dim(sp))
  vol <- NeuroVol(arr, sp)
  vec <- NeuroVec(array(0, c(dim(sp), 3L)), add_dim(sp, 3L))

  idx <- as.integer(c(1, 17, 63))
  ref <- grid_to_coord(sp, index_to_grid(sp, idx))

  expect_equal(index_to_coord(sp, idx), ref) # NeuroSpace, integer
  expect_equal(index_to_coord(sp, as.numeric(idx)), ref) # NeuroSpace, numeric
  expect_equal(index_to_coord(vol, idx), ref) # NeuroVol
  expect_equal(index_to_coord(vec, idx), ref) # NeuroVec
})

test_that("all coord_to_index dispatch paths agree", {
  sp <- aniso_space()
  vol <- NeuroVol(array(0, dim(sp)), sp)

  idx <- as.integer(c(1, 17, 63))
  world <- index_to_coord(sp, idx)

  expect_equal(as.vector(coord_to_index(sp, world)), idx) # NeuroSpace, matrix
  expect_equal(as.vector(coord_to_index(vol, world)), idx) # NeuroVol, matrix
  expect_equal(coord_to_index(sp, as.vector(world[1, ])), idx[1]) # NeuroSpace, numeric
})

test_that("vectorised conversion equals the scalar one", {
  sp <- real_space()
  idx <- as.integer(c(1, 999, 40000))

  batch <- index_to_coord(sp, idx)
  one_at_a_time <- t(vapply(idx, function(i) as.vector(index_to_coord(sp, i)), numeric(3)))

  expect_equal(unname(batch), unname(one_at_a_time))
})

test_that("conversions survive a negative-determinant (LAS) affine", {
  # The shipped mask runs LAS, so its first axis has a negative direction --
  # the case where a sign error in the offset would be easiest to miss.
  sp <- real_space()
  expect_lt(det(trans(sp)[1:3, 1:3]), 0)

  idx <- as.integer(c(1, 500, 12345, prod(dim(sp))))
  expect_equal(index_to_coord(sp, idx), grid_to_coord(sp, index_to_grid(sp, idx)))
  expect_equal(as.vector(coord_to_index(sp, index_to_coord(sp, idx))), idx)
})

test_that("grid and linear index round-trip independently of the affine", {
  for (sp in list(iso_space(), aniso_space(), oblique_space(), real_space())) {
    idx <- as.integer(seq(1, prod(dim(sp)), length.out = 50))
    expect_equal(as.vector(grid_to_index(sp, index_to_grid(sp, idx))), idx)
  }
})

# =============================================================================
# BUG 2: reorient()
# =============================================================================

test_that("reorient reaches all 48 orientations from any starting point", {
  starts <- list(iso_space(), aniso_space(), real_space())
  targets <- all_orientations()
  expect_length(targets, 48)

  for (sp in starts) {
    for (tgt in targets) {
      got <- affine_to_axcodes(trans(reorient(sp, tgt)))
      expect_identical(got, tgt,
        info = paste0(
          "from ", paste(affine_to_axcodes(trans(sp)), collapse = ""),
          " to ", paste(tgt, collapse = "")
        )
      )
    }
  }
})

test_that("asking for the orientation an image already has is a no-op", {
  # The sharpest form of the old bug: this returned LPI for an RAS input.
  for (sp in list(iso_space(), aniso_space(), real_space(), oblique_space())) {
    current <- affine_to_axcodes(trans(sp))
    expect_equal(trans(reorient(sp, current)), trans(sp))
  }
})

test_that("reorient does not return the inverse of the request", {
  # Pins the specific inversion: every code came back flipped.
  flip <- c(R = "L", L = "R", A = "P", P = "A", S = "I", I = "S")
  sp <- iso_space()

  for (tgt in all_orientations()) {
    got <- affine_to_axcodes(trans(reorient(sp, tgt)))
    expect_false(identical(got, unname(flip[tgt])),
      info = paste("target", paste(tgt, collapse = ""))
    )
  }
})

test_that("reorient preserves dimensions, voxel sizes and handedness magnitude", {
  sp <- aniso_space()
  for (tgt in all_orientations()) {
    out <- reorient(sp, tgt)
    expect_identical(dim(out), dim(sp))
    expect_equal(sort(spacing(out)), sort(spacing(sp)))
    expect_equal(abs(det(trans(out)[1:3, 1:3])), abs(det(trans(sp)[1:3, 1:3])))
  }
})

test_that("reorient is idempotent and composes", {
  sp <- real_space()
  for (tgt in list(c("R", "A", "S"), c("L", "P", "I"), c("A", "R", "S"), c("I", "P", "L"))) {
    once <- reorient(sp, tgt)
    twice <- reorient(once, tgt)
    expect_equal(trans(twice), trans(once))

    # Reaching a target via an intermediate orientation gives the same codes.
    via <- reorient(reorient(sp, c("L", "P", "I")), tgt)
    expect_identical(affine_to_axcodes(trans(via)), tgt)
  }
})

test_that("reorient accepts lower-case codes", {
  sp <- real_space()
  expect_identical(
    affine_to_axcodes(trans(reorient(sp, c("r", "a", "s")))),
    c("R", "A", "S")
  )
})

test_that("reorient works on a space carrying a time axis", {
  sp4 <- add_dim(iso_space(), 12L)
  out <- reorient(sp4, c("L", "P", "I"))

  expect_identical(affine_to_axcodes(trans(out)), c("L", "P", "I"))
  expect_equal(dim(out)[1:3], dim(sp4)[1:3])
})

test_that("reorient rejects malformed orientation codes", {
  sp <- iso_space()

  expect_error(reorient(sp, c("R", "A", "Q")))
  expect_error(reorient(sp, c("R", "A")))
  expect_error(reorient(sp, c("R", "A", "S", "T")))
  expect_error(reorient(sp, c("", "A", "S")))
})

test_that("reorient relabels the grid rather than moving the array", {
  # It changes where a voxel sits in world space -- that is the property that
  # makes it the wrong tool for rearranging data, and it should stay true.
  sp <- real_space()
  ras <- reorient(sp, c("R", "A", "S"))

  before <- as.vector(grid_to_coord(sp, matrix(c(1, 1, 1), nrow = 1)))
  after <- as.vector(grid_to_coord(ras, matrix(c(1, 1, 1), nrow = 1)))

  expect_false(isTRUE(all.equal(before, after)))
  expect_identical(dim(ras), dim(sp))
})

test_that("coordinate conversions stay consistent after reorienting", {
  # The two fixes interact: conversions must hold on a reoriented space too.
  sp <- reorient(real_space(), c("R", "A", "S"))
  idx <- as.integer(c(1, 777, 20000))

  expect_equal(index_to_coord(sp, idx), grid_to_coord(sp, index_to_grid(sp, idx)))
  expect_equal(as.vector(coord_to_index(sp, index_to_coord(sp, idx))), idx)
})

# =============================================================================
# Precision: grid_to_index() truncation vs affine rounding
# =============================================================================

test_that("grid_to_index still truncates genuinely fractional input", {
  # The near-integer snap added for the coord_to_index fix must not turn
  # grid_to_index into a rounding function. R's own array indexing truncates,
  # and callers passing 2.6 mean voxel 2.
  sp <- iso_space()

  expect_equal(grid_to_index(sp, matrix(c(2.0, 1, 1), nrow = 1)), 2)
  expect_equal(grid_to_index(sp, matrix(c(2.4, 1, 1), nrow = 1)), 2)
  expect_equal(grid_to_index(sp, matrix(c(2.6, 1, 1), nrow = 1)), 2)
  expect_equal(grid_to_index(sp, matrix(c(2.9, 1, 1), nrow = 1)), 2)
})

test_that("grid_to_index absorbs the affine's own rounding error", {
  # signif(trans, 7) means a coordinate that is arithmetically an exact voxel
  # centre can arrive just under the integer. Truncating that is an off-by-one.
  sp <- iso_space()

  expect_equal(grid_to_index(sp, matrix(c(3 - 1e-7, 1, 1), nrow = 1)), 3)
  expect_equal(grid_to_index(sp, matrix(c(3 - 1e-5, 1, 1), nrow = 1)), 3)
  expect_equal(grid_to_index(sp, matrix(c(3 + 1e-5, 1, 1), nrow = 1)), 3)
})

test_that("the two-step and shortcut world-to-index paths agree on real data", {
  # This is the pairing that broke: coord_to_index() rounds, while
  # grid_to_index(coord_to_grid(...)) truncated an affine-rounded value.
  sp <- real_space()
  idx <- as.integer(seq(1, prod(dim(sp)), length.out = 500))
  world <- index_to_coord(sp, idx)

  expect_equal(as.vector(coord_to_index(sp, world)), idx)
  expect_equal(as.vector(grid_to_index(sp, coord_to_grid(sp, world))), idx)
})

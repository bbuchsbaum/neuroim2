context("connected components: compiled labeller vs the R reference")

library(neuroim2)

# The interpreted two-pass union-find that conn_comp_3D() used before the
# compiled labeller landed, kept verbatim as the equivalence reference. It is
# ~1000x slower, which is why it is not the implementation any more, but it is
# the definition of the output the compiled version has to reproduce -- labels,
# sizes and the size-descending numbering with its tie-break.
ref_conn_comp_3D <- function(mask, connect = c("26-connect", "18-connect", "6-connect")) {
  # Input validation with more informative messages
  if (!is.array(mask) || length(dim(mask)) != 3) {
    cli::cli_abort("{.arg mask} must be a 3D array.")
  }
  if (!is.logical(mask[1])) {
    cli::cli_abort("{.arg mask} must be a logical array.")
  }

  connect <- match.arg(connect)

  # Use integer arrays for better memory efficiency
  nodes <- seq_len(length(mask) %/% 9)  # Initialize nodes with identity mapping
  labels <- array(0L, dim(mask))        # Use integer array
  
  DIM <- dim(mask)
  
  # Pre-compute neighborhood patterns more efficiently
  local.mask <- switch(connect,
    "6-connect" = {
      matrix(c(
        -1,0,0,  1,0,0,  0,-1,0,  0,1,0,  0,0,-1,  0,0,1
      ), ncol=3, byrow=TRUE)
    },
    "18-connect" = {
      matrix(c(
        -1,0,0,  1,0,0,   0,-1,0,  0,1,0,   0,0,-1,  0,0,1,  # 6-neighbors
        -1,-1,0, -1,1,0,  1,-1,0,  1,1,0,                     # edge neighbors xy
        -1,0,-1, -1,0,1,  1,0,-1,  1,0,1,                     # edge neighbors xz
        0,-1,-1, 0,-1,1,  0,1,-1,  0,1,1                      # edge neighbors yz
      ), ncol=3, byrow=TRUE)
    },
    # 26-connect: more efficient than expand.grid
    matrix(c(
      -1,-1,-1, -1,-1,0, -1,-1,1,  -1,0,-1, -1,0,0, -1,0,1,  -1,1,-1, -1,1,0, -1,1,1,
       0,-1,-1,  0,-1,0,  0,-1,1,   0,0,-1,         0,0,1,    0,1,-1,  0,1,0,  0,1,1,
       1,-1,-1,  1,-1,0,  1,-1,1,   1,0,-1,  1,0,0, 1,0,1,    1,1,-1,  1,1,0,  1,1,1
    ), ncol=3, byrow=TRUE)
  )
  
  # Transpose once for efficiency
  tlocal.mask <- t(local.mask)
  
  # Optimized neighbor function with fewer allocations
  neighbors <- function(vox) {
    # Add current voxel to each neighbor offset
    vox.hood <- t(tlocal.mask + vox)
    
    # Handle boundary cases more efficiently
    if (any(vox == 1) || any(vox == DIM)) {
      valid <- vox.hood[,1] >= 1 & vox.hood[,1] <= DIM[1] &
               vox.hood[,2] >= 1 & vox.hood[,2] <= DIM[2] &
               vox.hood[,3] >= 1 & vox.hood[,3] <= DIM[3]
      vox.hood <- vox.hood[valid, , drop=FALSE]
    }
    
    # Get labeled neighbors
    has_label <- labels[cbind(vox.hood[,1], vox.hood[,2], vox.hood[,3])] != 0
    vox.hood[has_label, , drop=FALSE]
  }
  
  # Path compression in find() for better performance
  find <- function(i) {
    root <- i
    # Find root
    while (nodes[root] != root) {
      root <- nodes[root]
    }
    # Path compression
    while (nodes[i] != root) {
      parent <- nodes[i]
      nodes[i] <- root
      i <- parent
    }
    root
  }
  
  # First pass: Initial labeling
  nextlabel <- 1L
  grid <- which(mask > 0, arr.ind=TRUE)  # More efficient than .indexToGrid
  
  for (i in seq_len(nrow(grid))) {
    vox <- grid[i,]
    nabes <- neighbors(vox)
    
    if (nrow(nabes) == 0) {
      nodes[nextlabel] <- nextlabel
      labels[vox[1], vox[2], vox[3]] <- nextlabel
    } else {
      L <- labels[cbind(nabes[,1], nabes[,2], nabes[,3])]
      ML <- min(L)
      labels[vox[1], vox[2], vox[3]] <- ML
      nodes[nextlabel] <- ML
      
      # Union operation with unique labels for efficiency
      for (lab in unique(L)) {
        rootx <- find(lab)
        nodes[rootx] <- find(ML)
      }
    }
    
    nextlabel <- nextlabel + 1L
  }
  
  # Second pass: Resolve labels more efficiently
  label_indices <- which(labels > 0)
  labels[label_indices] <- vapply(labels[label_indices], find, integer(1))
  
  # Create output volumes efficiently
  labs <- labels[label_indices]
  clusters <- sort(table(labs), decreasing=TRUE)
  
  # Handle case of no components
  if (length(clusters) == 0) {
    return(list(
      index = array(0L, DIM),
      size = array(0L, DIM)
    ))
  }
  
  # Size volume
  SVol <- array(0L, DIM)
  SVol[label_indices] <- clusters[as.character(labs)]
  
  # Index volume with consecutive indices
  indices <- seq_along(clusters)
  names(indices) <- names(clusters)
  IVol <- array(0L, DIM)
  IVol[label_indices] <- indices[as.character(labs)]
  
  list(index=IVol, size=SVol)
}

random_mask <- function(dim, frac, seed) {
  set.seed(seed)
  array(runif(prod(dim)) < frac, dim)
}

test_that("the compiled labeller reproduces the R reference exactly", {
  cases <- list(
    list(dim = c(8L, 9L, 7L),  frac = 0.30),
    list(dim = c(8L, 9L, 7L),  frac = 0.05),
    list(dim = c(8L, 9L, 7L),  frac = 0.75),
    list(dim = c(12L, 12L, 12L), frac = 0.20),
    list(dim = c(1L, 9L, 7L),  frac = 0.40),   # degenerate axis
    list(dim = c(9L, 1L, 1L),  frac = 0.50),
    list(dim = c(5L, 4L, 20L), frac = 0.15),
    list(dim = c(3L, 3L, 3L),  frac = 1.00),   # everything in mask
    list(dim = c(4L, 4L, 4L),  frac = 0.00)    # nothing in mask
  )
  for (ci in seq_along(cases)) {
    cs <- cases[[ci]]
    for (seed in 1:3) {
      m <- random_mask(cs$dim, cs$frac, seed * 100 + ci)
      for (cn in c("6-connect", "18-connect", "26-connect")) {
        got <- conn_comp_3D(m, connect = cn)
        want <- ref_conn_comp_3D(m, connect = cn)
        info <- paste(paste(cs$dim, collapse = "x"), cs$frac, seed, cn)
        expect_identical(got$index, want$index, info = info)
        expect_identical(got$size, want$size, info = info)
      }
    }
  }
})

test_that("connectivity means what it says", {
  # Two voxels touching only at a corner: one component under 26-connect,
  # two under 18- and 6-connect.
  m <- array(FALSE, c(3L, 3L, 3L))
  m[1, 1, 1] <- TRUE; m[2, 2, 2] <- TRUE
  expect_equal(max(conn_comp_3D(m, "26-connect")$index), 1L)
  expect_equal(max(conn_comp_3D(m, "18-connect")$index), 2L)
  expect_equal(max(conn_comp_3D(m, "6-connect")$index), 2L)

  # Touching along an edge: one under 26 and 18, two under 6.
  m2 <- array(FALSE, c(3L, 3L, 3L))
  m2[1, 1, 1] <- TRUE; m2[2, 2, 1] <- TRUE
  expect_equal(max(conn_comp_3D(m2, "26-connect")$index), 1L)
  expect_equal(max(conn_comp_3D(m2, "18-connect")$index), 1L)
  expect_equal(max(conn_comp_3D(m2, "6-connect")$index), 2L)
})

test_that("components are numbered by decreasing size", {
  m <- array(FALSE, c(10L, 10L, 10L))
  m[1:4, 1:4, 1:4] <- TRUE     # 64 voxels
  m[9:10, 9:10, 9:10] <- TRUE  # 8 voxels
  m[1, 10, 1] <- TRUE          # 1 voxel
  cc <- conn_comp_3D(m, "6-connect")
  expect_equal(max(cc$index), 3L)
  expect_equal(sum(cc$index == 1L), 64L)
  expect_equal(sum(cc$index == 2L), 8L)
  expect_equal(sum(cc$index == 3L), 1L)
  expect_equal(cc$size[1, 1, 1], 64L)
  expect_equal(cc$size[1, 10, 1], 1L)
  expect_true(all(cc$index[!m] == 0L))
  expect_true(all(cc$size[!m] == 0L))
})

test_that("an empty mask yields empty output", {
  m <- array(FALSE, c(4L, 5L, 3L))
  cc <- conn_comp_3D(m)
  expect_equal(dim(cc$index), c(4L, 5L, 3L))
  expect_true(all(cc$index == 0L))
  expect_true(all(cc$size == 0L))
})

test_that("the labeller rejects hostile arguments", {
  expect_error(conn_comp_3D(array(1, c(2, 2, 2))), "logical")
  expect_error(conn_comp_3D(array(TRUE, c(2, 2))), "3D array")
  expect_error(conn_comp_3D(array(TRUE, c(2, 2, 2)), connect = "nope"))
  ns <- asNamespace("neuroim2")
  expect_error(ns$conn_comp_labels_cpp(rep(TRUE, 8), c(2L, 2L, 2L), 7L),
               "must be 6, 18 or 26")
  expect_error(ns$conn_comp_labels_cpp(rep(TRUE, 8), c(2L, 2L), 26L),
               "length 3")
  expect_error(ns$conn_comp_labels_cpp(rep(TRUE, 8), c(3L, 2L, 2L), 26L),
               "elements but")
})

# ---------------------------------------------------------------------------
# conn_comp(): the wrapper around the labeller was rewritten at the same time,
# to stop re-deriving coordinate matrices per component. Its output must be
# unchanged, so the previous body is kept here as the reference.
# ---------------------------------------------------------------------------

ref_conn_comp <- function(x, threshold = 0, cluster_table = TRUE,
                          local_maxima = TRUE, local_maxima_dist = 15, ...) {
  prune <- get(".pruneCoords", envir = asNamespace("neuroim2"))
  mask <- (x > threshold)
  stopifnot(any(mask))
  comps <- conn_comp_3D(mask, ...)
  grid <- as.data.frame(index_to_grid(mask, which(mask > 0)))
  colnames(grid) <- c("x", "y", "z")
  locations <- split(grid, comps$index[comps$index > 0])
  ret <- list(index = ClusteredNeuroVol(mask, clusters = comps$index[mask > 0]),
              size = NeuroVol(comps$size, space(x)), voxels = locations)
  if (cluster_table) {
    maxima <- do.call(rbind, lapply(locations, function(loc) {
      if (nrow(loc) == 1) loc else loc[which.max(x[as.matrix(loc)]), ]
    }))
    N <- comps$size[as.matrix(maxima)]
    ret$cluster_table <- data.frame(index = 1:NROW(maxima), x = maxima[, 1],
                                    y = maxima[, 2], z = maxima[, 3], N = N,
                                    Area = N * prod(spacing(x)),
                                    value = x[as.matrix(maxima)])
  }
  if (local_maxima) {
    coord.sets <- lapply(locations, function(loc) sweep(as.matrix(loc), 2, spacing(x), "*"))
    loc.max <- do.call(rbind, mapply(function(cset, i) {
      idx <- prune(as.matrix(cset), x[as.matrix(locations[[i]])], mindist = local_maxima_dist)
      cbind(i, as.matrix(locations[[i]])[idx, , drop = FALSE])
    }, coord.sets, seq_along(coord.sets), SIMPLIFY = FALSE))
    loc.max <- cbind(loc.max, x[loc.max[, 2:4, drop = FALSE]])
    row.names(loc.max) <- 1:NROW(loc.max)
    colnames(loc.max) <- c("index", "x", "y", "z", "value")
    ret$local_maxima <- loc.max
  }
  ret
}

same_conn_comp <- function(a, b) {
  identical(a$voxels, b$voxels) &&
    identical(as.array(a$size), as.array(b$size)) &&
    identical(a$index@clusters, b$index@clusters) &&
    isTRUE(all.equal(a$cluster_table, b$cluster_table)) &&
    isTRUE(all.equal(a$local_maxima, b$local_maxima)) &&
    identical(sort(ls(a$index@cluster_map)), sort(ls(b$index@cluster_map))) &&
    all(vapply(ls(a$index@cluster_map),
               function(k) identical(a$index@cluster_map[[k]], b$index@cluster_map[[k]]),
               logical(1)))
}

test_that("conn_comp() output is unchanged by the rewrite", {
  for (dm in list(c(12L, 13L, 11L), c(8L, 9L, 7L))) {
    for (seed in 1:2) {
      set.seed(seed * 31 + dm[1])
      sp <- NeuroSpace(dm, c(2, 2, 2.5))
      v <- NeuroVol(array(rnorm(prod(dm)), dm), sp)
      for (thr in c(-0.5, 0, 0.8, 1.5)) {
        for (cn in c("6-connect", "26-connect")) {
          got <- suppressWarnings(conn_comp(v, threshold = thr, connect = cn))
          want <- suppressWarnings(ref_conn_comp(v, threshold = thr, connect = cn))
          expect_true(same_conn_comp(got, want),
                      info = paste(paste(dm, collapse = "x"), seed, thr, cn))
        }
      }
    }
  }
})

test_that("a single cluster still gets integer row names in the cluster table", {
  # A 1-row matrix subset picks up the column name, which used to turn the
  # cluster table's row names from 1 into "x".
  sp <- NeuroSpace(c(6L, 6L, 6L), c(2, 2, 2))
  d <- array(0, c(6, 6, 6)); d[2:4, 2:4, 2:4] <- 1
  cc <- suppressWarnings(conn_comp(NeuroVol(d, sp), threshold = 0.5))
  expect_equal(nrow(cc$cluster_table), 1L)
  expect_identical(attr(cc$cluster_table, "row.names"), 1L)
})

test_that("ClusteredNeuroVol builds the same cluster map as the name-by-name loop", {
  mask <- LogicalNeuroVol(array(TRUE, c(8L, 8L, 8L)), NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1)))
  set.seed(3)
  clusters <- sample(1:40, 512, replace = TRUE)
  cv <- ClusteredNeuroVol(mask, clusters)

  want <- split(which(as.array(mask)), clusters)
  expect_setequal(ls(cv@cluster_map), names(want))
  for (k in names(want)) expect_identical(cv@cluster_map[[k]], want[[k]])
})

#' Iterate over neuroimaging data in storage-aware blocks
#'
#' Construct a lazy, repeatable list of decoded data blocks. The caller owns
#' computation and any accumulated results.
#'
#' @param x A concrete DenseNeuroVec, SparseNeuroVec or BigNeuroVec.
#' @param along Traversal direction: `"auto"`, `"time"` or `"voxel"`.
#' @param budget Positive finite working-buffer budget in bytes. Default: 256 MiB.
#'
#' @return A non-memoising [deflist::deflist]. Each element is a list with
#'   `values` (time-by-voxel matrix), `voxels` and `times` (original one-based
#'   indices), `along`, `space` (the original 3D spatial geometry), and
#'   `volume_labels` for its times. The `plan` attribute reports traversal,
#'   block length, block count, budget and estimated maximum buffer bytes.
#'
#' @details
#' Construction reads no intensity data. Auto traverses dense arrays by time
#' and sparse/FBM stores by voxel, matching their storage layout. An explicit
#' direction is honored. Sparse positions outside the mask remain zero.
#'
#' The planner conservatively allows 128 bytes per selected value for decoded
#' output, indexing, gather and transpose buffers, 32 bytes per axis index, and
#' the deferred-list pointer array. This is a buffer estimate, not a process-RSS
#' limit: existing source storage, R object/allocator/GC overhead, mapped pages,
#' and blocks retained by the caller are excluded. Source opening and eager FBM
#' construction precede iteration. At least one full volume (time traversal) or
#' one complete time course (voxel traversal) must fit. Too-small budgets fail
#' during construction, before reading. Blocks are recomputed on repeated access;
#' the iterator does not cache realized blocks. Keep the source unchanged while
#' iterating, including an FBM's externally mutable backing file.
#'
#' MappedNeuroVec, FileBackedNeuroVec, NeuroVecSeq, 5D, clustered, bucket and
#' custom subclasses are currently rejected explicitly. Their access or scratch
#' costs require separate adapters; no implicit dense conversion is performed.
#'
#' @examples
#' x <- DenseNeuroVec(array(as.numeric(1:48), c(2, 3, 2, 4)),
#'                   NeuroSpace(c(2, 3, 2, 4)))
#' blocks <- vec_blocks(x, budget = 16384)
#' totals <- numeric(12)
#' for (i in seq_along(blocks)) {
#'   b <- blocks[[i]]
#'   totals[b$voxels] <- totals[b$voxels] + colSums(b$values)
#' }
#' totals / dim(x)[4]
#' @md
#' @export
vec_blocks <- function(x, along = c("auto", "time", "voxel"),
                       budget = 256 * 1024^2) {
  along <- match.arg(along)
  if (!is.numeric(budget) || length(budget) != 1L ||
      !is.finite(budget) || budget <= 0) {
    stop("'budget' must be one positive finite number of bytes.", call. = FALSE)
  }
  storage <- .vec_block_storage(x)
  d <- dim(x)
  if (length(d) != 4L || any(!is.finite(d)) || any(d <= 0)) {
    stop("vec_blocks requires a non-empty 4D source.", call. = FALSE)
  }
  axis <- if (along == "auto") storage$along else along
  nv <- prod(d[1:3]); nt <- d[4]
  count <- if (axis == "time") nt else nv
  width <- if (axis == "time") nv else nt
  # The fixed opposite-axis indices and a small fixed workspace allowance.
  fixed <- 8192 + 32 * width
  per_unit <- 128 * width + 32
  size <- min(count, floor((budget - fixed) / per_unit))
  if (size < 1) {
    stop("'budget' cannot hold the minimum legal block (one full ",
         if (axis == "time") "volume" else "time course", ").", call. = FALSE)
  }
  # Account for deflist's O(number of blocks) pointer array, without creating it.
  # Shrinking monotonically reaches a feasible block or rejects the budget.
  repeat {
    nblocks <- ceiling(count / size)
    estimated <- fixed + per_unit * size + 8 * nblocks
    if (estimated <= budget) break
    next_size <- floor((budget - fixed - 8 * nblocks) / per_unit)
    if (next_size < 1) {
      stop("'budget' cannot hold a block and the deferred-list index.", call. = FALSE)
    }
    size <- next_size
  }
  plan <- list(along = axis, block_length = size, nblocks = nblocks,
               budget = budget, estimated_buffer_bytes = estimated,
               nvoxels = nv, ntimes = nt, backend = storage$name)
  f <- .vec_block_reader(x, plan, drop_dim(space(x)), volume_labels(x))
  blocks <- deflist::deflist(f, nblocks, memoise = FALSE)
  attr(blocks, "plan") <- plan
  blocks
}

# Deliberately exact-class checks: unknown lazy subclasses may override access
# with unbounded work, so their memory costs cannot inherit these guarantees.
.vec_block_storage <- function(x) {
  for (name in c("DenseNeuroVec", "SparseNeuroVec", "BigNeuroVec")) {
    if (identical(class(x), methods::getClass(name)@className)) {
      return(list(name = name, along = if (name == "DenseNeuroVec") "time" else "voxel"))
    }
  }
  stop("vec_blocks does not yet support class '", class(x)[1L], "'.", call. = FALSE)
}

.vec_block_reader <- function(x, plan, geometry, labels) {
  force(x); force(plan); force(geometry); force(labels)
  function(i) {
    if (length(i) != 1L || !is.finite(i) || i != floor(i) ||
        i < 1 || i > plan$nblocks) {
      stop("Block index must be one whole number within the iterator.", call. = FALSE)
    }
    first <- (i - 1) * plan$block_length + 1
    count <- if (plan$along == "time") plan$ntimes else plan$nvoxels
    selected <- seq.int(first, min(count, first + plan$block_length - 1))
    voxels <- if (plan$along == "voxel") selected else seq_len(plan$nvoxels)
    times <- if (plan$along == "time") selected else seq_len(plan$ntimes)
    if (plan$along == "voxel" && plan$backend != "DenseNeuroVec") {
      values <- matrix(series(x, as.integer(voxels), drop = FALSE),
                       nrow = length(times), ncol = length(voxels))
    } else {
      # Paired linear elements of this Cartesian block, volume-major. Dense
      # primitive subsetting avoids extracting/coercing its entire data part.
      idx <- rep(voxels, times = length(times)) +
        rep((times - 1) * plan$nvoxels, each = length(voxels))
      vals <- if (plan$backend == "DenseNeuroVec") .subset(x, idx) else linear_access(x, idx)
      values <- t(matrix(vals, nrow = length(voxels), ncol = length(times)))
    }
    list(values = values, voxels = voxels, times = times, along = plan$along,
         space = geometry, volume_labels = .subset_volume_labels(labels, times))
  }
}

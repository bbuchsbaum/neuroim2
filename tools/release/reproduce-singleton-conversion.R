args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L, args[1] %in% c("before", "after"))
if (length(args) > 1L) .libPaths(c(args[2], .libPaths()))
library(neuroim2)
cat("neuroim2:", as.character(packageVersion("neuroim2")), "\n")

results <- list()
for (nt in c(1L, 2L, 3L, 5L)) for (nv in c(1L, 2L, 3L, 5L)) {
  dims <- c(3L, 3L, 1L, nt)
  sp <- NeuroSpace(dims)
  values <- outer(seq_len(9L), seq_len(nt), function(v, t) 100 * t + v)
  x <- DenseNeuroVec(array(values, dims), sp)
  voxels <- seq_len(nv) * 2L - 1L
  keep <- array(FALSE, dims[1:3]); keep[voxels] <- TRUE
  for (type in c("numeric", "logical volume")) {
    selection <- if (type == "numeric") as.numeric(voxels) else
      LogicalNeuroVol(keep, drop_dim(sp))
    error <- ""
    passed <- tryCatch({
      y <- as.sparse(x, selection)
      full <- values; full[-voxels, ] <- 0
      expected <- outer(seq_len(nt), voxels, function(t, v) 100 * t + v)
      isTRUE(all.equal(as.matrix(y), full)) &&
        isTRUE(all.equal(series(y, voxels, drop = FALSE), expected)) &&
        identical(as.integer(dim(y)), dims)
    }, error = function(e) { error <<- conditionMessage(e); FALSE })
    results[[length(results) + 1L]] <- data.frame(
      times = nt, voxels = nv, mask = type, passed = passed, error = error)
  }
}
results <- do.call(rbind, results)
print(results, row.names = FALSE)
cat("Passed:", sum(results$passed), "Failed:", sum(!results$passed), "\n")
if (args[1] == "before") {
  # Both singleton axes fail on master; rectangular/square controls work.
  stopifnot(all(results$passed == (results$times > 1L & results$voxels > 1L)))
}
if (args[1] == "after") stopifnot(all(results$passed))

# Run from package root: Rscript dev/bench/iterator/read-costs.R
# Synthetic data only; preparation is timed separately from traversal.
devtools::load_all(quiet = TRUE)
source("tests/testthat/helper-backend-fixture.R")
run_costs <- function() {
  f <- backend_binary_fixture("SHORT", c(16L, 16L, 8L, 24L))
  on.exit(unlink(f$path), add = TRUE)
  backing <- tempfile()
  on.exit(unlink(paste0(backing, c(".bk", ".rds"))), add = TRUE)
  prepare <- function(fun) {
    elapsed <- system.time(x <- fun())[["elapsed"]]
    list(x = x, elapsed = elapsed)
  }
  dense <- prepare(function() read_vec(f$path))
  stores <- list(dense = dense,
    mapped = prepare(function() MappedNeuroVec(f$path)),
    filebacked = prepare(function() FileBackedNeuroVec(f$path)),
    sparse = prepare(function() SparseNeuroVec(f$expected, space(dense$x),
                                               array(TRUE, dim(f$expected)[1:3]))),
    big = prepare(function() BigNeuroVec(f$expected, space(dense$x),
                                        array(TRUE, dim(f$expected)[1:3]),
                                        backingfile = backing)))
  on.exit(mmap::munmap(stores$mapped$x@filemap), add = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$volumes <- 0
  # Trace counts requested raw volumes at the real compiled-reader boundary.
  options(neuroim2.cost.counter = counter)
  on.exit(options(neuroim2.cost.counter = NULL), add = TRUE)
  trace("nifti_read_volumes_cpp", where = asNamespace("neuroim2"),
        tracer = quote({
          cc <- getOption("neuroim2.cost.counter")
          cc$volumes <- cc$volumes + length(vols)
        }), print = FALSE)
  on.exit(untrace("nifti_read_volumes_cpp", where = asNamespace("neuroim2")), add = TRUE)
  rows <- list()
  for (name in names(stores)) for (along in c("time", "voxel")) {
    x <- stores[[name]]$x
    nv <- prod(dim(x)[1:3]); nt <- dim(x)[4]
    groups <- if (along == "time") split(seq_len(nt), ceiling(seq_len(nt) / 3)) else
      split(seq_len(nv), ceiling(seq_len(nv) / 256))
    counter$volumes <- 0
    mem <- tempfile()
    Rprofmem(mem)
    timing <- system.time(for (g in groups) {
      # Same public-accessor routes contemplated by the block planner.
      value <- if (along == "time") {
        x[, , , g, drop = FALSE]
      } else series(x, as.integer(g), drop = FALSE)
      stopifnot(length(value) == if (along == "time") nv * length(g) else nt * length(g))
    })[["elapsed"]]
    Rprofmem(NULL)
    alloc <- suppressWarnings(as.numeric(sub(" .*", "", readLines(mem))))
    unlink(mem)
    rows[[length(rows) + 1L]] <- data.frame(
      backend = name, along = along, preparation_seconds = stores[[name]]$elapsed,
      elapsed_seconds = timing, allocated_bytes = sum(alloc, na.rm = TRUE),
      largest_allocation = max(alloc, na.rm = TRUE),
      raw_volumes = counter$volumes, raw_payload_bytes = counter$volumes * nv * 2)
  }
  print(do.call(rbind, rows), row.names = FALSE)
  # Audit primitive dense subsetting: no full-source data-part extraction.
  mem <- tempfile(); Rprofmem(mem)
  z <- .subset(dense$x, 1:2048)
  Rprofmem(NULL)
  stopifnot(identical(as.numeric(z), as.numeric(f$expected)[1:2048]))
  a <- suppressWarnings(as.numeric(sub(" .*", "", readLines(mem))))
  cat("dense .subset largest allocation:", max(a, na.rm = TRUE), "\n")
  unlink(mem)
}
run_costs()

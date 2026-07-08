## Micro-benchmarks for hot data-access paths in
## array_like.R / neurovec.R / sparse_neurovec.R
## Skip attaching if a dev build was already loaded via pkgload::load_all().
if (!isTRUE(getOption("neuroim2.bench.loaded"))) {
  suppressPackageStartupMessages(library(neuroim2))
}

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  bench <- function(expr, times = 50) {
    expr <- substitute(expr)
    t <- replicate(times, system.time(eval(expr, parent.frame()))[["elapsed"]])
    median(t) * 1000
  }
} else {
  bench <- function(expr, times = 50) {
    expr <- substitute(expr)
    res <- microbenchmark::microbenchmark(list = list(eval = expr),
                                          times = times, unit = "ms")
    stats::median(res$time) / 1e6
  }
}

set.seed(42)

## ---- Dense 4D ----
dims <- c(64, 64, 40, 60)
arr <- array(rnorm(prod(dims)), dims)
dvec <- DenseNeuroVec(arr, NeuroSpace(dims))

## Dense 3D volume
vol <- dvec[[1]]

## Sparse
mask_arr <- array(runif(prod(dims[1:3])) > 0.85, dims[1:3])
svec <- as.sparse(dvec, LogicalNeuroVol(mask_arr, drop_dim(space(dvec))))

## Coordinate sets
ncoord <- 2000
coords <- cbind(sample.int(dims[1], ncoord, replace = TRUE),
                sample.int(dims[2], ncoord, replace = TRUE),
                sample.int(dims[3], ncoord, replace = TRUE))
lin_idx <- sample.int(prod(dims), 50000)
spatial_lin <- sample.int(prod(dims[1:3]), 20000)

res <- list()
res[["dense_vol_slice[i,]"]]        <- bench(vol[10:50, ])
res[["dense_vol_full[,,]"]]         <- bench(vol[, , ], times = 30)
res[["dense_vec_4d[i,j,k,m]"]]      <- bench(dvec[10:40, 10:40, 5:30, 1:20], times = 20)
res[["dense_vec_linear[idx]"]]      <- bench(dvec[lin_idx])
res[["dense_series_matrix"]]        <- bench(series(dvec, coords), times = 20)
res[["sparse_linear_access"]]       <- tryCatch(bench(linear_access(svec, lin_idx)),
                                                error = function(e) NA_real_)
res[["sparse_[i,j,k,m]"]]           <- bench(svec[10:40, 10:40, 5:30, 1:20], times = 20)
res[["sparse_series_matrix"]]       <- bench(series(svec, coords), times = 20)
res[["ilv_lookup"]]                 <- bench(neuroim2:::lookup(svec@map, spatial_lin))

cat(sprintf("%-32s %10s\n", "benchmark", "median_ms"))
for (nm in names(res)) cat(sprintf("%-32s %10.3f\n", nm, res[[nm]]))

## Correctness cross-checks against the base R array, so refactors stay honest.
dense_arr <- as.array(dvec)
sparse_arr <- as.array(svec)
stopifnot(identical(dim(dense_arr), as.integer(dims)))
chk <- function(label, a, b) {
  d <- max(abs(as.vector(a) - as.vector(b)))
  cat(sprintf("CHK %-28s maxdiff=%.3g\n", label, d))
  invisible(d)
}
chk("dense_vol_slice",  vol[10:50, ], dense_arr[10:50, , , 1])
chk("dense_vol_full",   vol[, , ],    dense_arr[, , , 1])
chk("dense_vec_4d",     dvec[10:40, 10:40, 5:30, 1:20], dense_arr[10:40, 10:40, 5:30, 1:20])
chk("dense_vec_linear", dvec[lin_idx], dense_arr[lin_idx])
chk("sparse_4d",        svec[10:40, 10:40, 5:30, 1:20], sparse_arr[10:40, 10:40, 5:30, 1:20])
chk("sparse_linear",    linear_access(svec, lin_idx), sparse_arr[lin_idx])
sc <- series(svec, coords)
sc_ref <- vapply(seq_len(nrow(coords)),
                 function(r) sparse_arr[coords[r,1], coords[r,2], coords[r,3], ],
                 numeric(dims[4]))
chk("sparse_series",    sc, sc_ref)
dc <- series(dvec, coords)
dc_ref <- vapply(seq_len(nrow(coords)),
                 function(r) dense_arr[coords[r,1], coords[r,2], coords[r,3], ],
                 numeric(dims[4]))
chk("dense_series",     dc, dc_ref)

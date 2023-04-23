

bench_chunk <- function(cdim, compression=0, nbit=TRUE) {
  bv <- array(rnorm(96*96), c(64,64,10,200))
  bv <- NeuroVec(bv, NeuroSpace(dim=dim(bv)))
  fname <- paste0(tempfile(), ".h5")
  hfile <- to_nih5_vec(bv,fname, chunk_dim=cdim, compression = compression, nbit=nbit)
  print(hfile[["data/elements"]]$get_storage_size())
  hfile$flush()
  hfile$close_all()


  hvec <- H5NeuroVec(fname)

  mask <- NeuroVol(array(rep(1, prod(dim(bv)[1:3])), dim(bv)[1:3]), NeuroSpace(dim=dim(bv)[1:3]) )
  grid <- index_to_grid(mask, 1:prod(dim(bv)[1:3]))
  kres <- kmeans(grid, centers=200, iter.max=500)
  kvol <- ClusteredNeuroVol(mask, kres$cluster)
  cds <- split_clusters(mask, kres$cluster)
  ret <- microbenchmark::microbenchmark(lapply(cds, function(x) {
    tmp <- series(hvec, coords(x))
    1
  }), times=1)

  ret$time/1000000000

}

g <- expand.grid(1:5, 1:5, 1:5)
ret <- lapply(1:nrow(g), function(i) {
  print(i)
  ret <- bench_chunk(c(g[i,1], g[i,2], g[i,3], 200), compression=0, nbit=FALSE)
})


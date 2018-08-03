

CachedSparseNeuroVec <- function(file_name, bucket_vol) {
  mask <- bucket_vol@mask
  meta <- read_header(file_name)
  assert_that(length(meta@dims) == 4)
  assert_that(all(meta@dims[1:3] == dim(bucket_vol)))

  indices <- which(mask> 0)
  map <- IndexLookupVol(space(mask), indices=indices)

  len <- prod(dim(mask))
  sp <- add_dim(space(mask), meta@dims[4])

  new("CachedSparseNeuroVec",
      space=sp,
      meta=meta,
      mask=mask,
      map=map,
      bucket_vol=bucket_vol,
      cache_size=as.integer(5),
      cache=sparseMatrix(i=numeric(), j=numeric(), dims=c(meta@dims[4], len)),
      cache_list=vector(length(vec@bucket_vol@cluster_map),mode="list"))

}


#' @export
#' @rdname series-methods
setMethod(f="series", signature=signature(x="CachedSparseNeuroVec", i="matrix"),
          def=function(x,i) {
            idx <- grid_to_index(x@mask, i)

            out <- matrix(0, dim(x)[4], length(idx))


            first <- min(idx)
            last <- max(idx)

            f <- series_reader(x@meta@data_file)
            m <- f(first, last)
            keep <- seq(first, last) %in% idx
            m[,keep]
          })



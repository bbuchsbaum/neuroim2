#' MappedNeuroVecSource
#'
#' constructs a MappedNeuroVecSource object
#'
#' @param file_name the name of the image file to be memory mapped
#' @export
#' @rdname MappedNeuroVecSource-class
MappedNeuroVecSource <- function(file_name) {
  meta <- read_header(file_name)
  stopifnot(length(dim(meta)) >= 3)
  new("MappedNeuroVecSource", meta_info=meta)
}


#' MappedNeuroVec
#'
#' constructs a MappedNeuroVec object
#'
#' @param file_name the name of the 4D image file that containing the memory-mapped data source.
#' @export
#' @examples
#' @rdname MappedNeuroVec-class
MappedNeuroVec <- function(file_name) {
  src <- MappedNeuroVecSource(file_name)
  load_data(src)
}


#' @export
#' @rdname load_data-methods
setMethod(f="load_data", signature=c("MappedNeuroVecSource"),
          def=function(x) {
            meta <- x@meta_info
            fmap <- mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), prot=mmap::mmapFlags("PROT_READ"))
            offset <- meta@data_offset/.getDataSize(meta@data_type)

            bspace <- NeuroSpace(dim(meta), meta@spacing,
                                 meta@origin, meta@spatial_axes)

            new("MappedNeuroVec", space=bspace, filemap=fmap, offset=as.integer(offset))

          })



setMethod(f="linear_access", signature=signature(x = "MappedNeuroVec", i = "numeric"),
          def=function (x, i) {
            idx <- i + x@offset
            x@filemap[idx]
          })



#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "MappedNeuroVec", i = "numeric", j = "numeric"),
          def = function (x, i, j, k, m, ..., drop = TRUE) {
            if (missing(k))
              k = 1:(dim(x)[3])
            if (missing(m)) {
              m <- 1:(dim(x)[4])
            }

            vmat <- expand.grid(i=i, j=j, k=k, m=m)
            idx <- .gridToIndex(dim(x), vmat)
            vals <- linear_access(x,idx)

            ret <- array(vals, c(length(i), length(j), length(k), length(m)))

            if (drop) {
              drop(ret)
            } else {
              ret
            }


          })

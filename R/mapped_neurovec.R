#' MappedNeuroVecSource
#'
#' constructs a MappedNeuroVecSource object
#'
#' @param file_name the name of the image file to be memory mapped
#' @export
#' @rdname MappedNeuroVecSource-class
#' @return a new instance of type \code{MappedNeuroVecSource}
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
#' @rdname MappedNeuroVec-class
#' @return a new instance of type \code{MappedNeuroVec}
MappedNeuroVec <- function(file_name) {
  src <- MappedNeuroVecSource(file_name)
  load_data(src)
}


#' @keywords internal
#' @noRd
setMethod(f="load_data", signature=c(x="MappedNeuroVecSource"),
          def=function(x) {
            meta <- x@meta_info
            fmap <- mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), prot=mmap::mmapFlags("PROT_READ"))
            offset <- meta@data_offset/.getDataSize(meta@data_type)

            bspace <- NeuroSpace(dim(meta), meta@spacing,
                                 meta@origin, meta@spatial_axes, trans=trans(meta))

            new("MappedNeuroVec", space=bspace, filemap=fmap, offset=as.integer(offset))

          })


#' @noRd
setMethod(f="linear_access", signature=signature(x = "MappedNeuroVec", i = "numeric"),
          def=function (x, i) {
            idx <- i + x@offset
            x@filemap[idx]
          })



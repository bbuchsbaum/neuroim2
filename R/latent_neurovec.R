
#' @keywords internal
LatentNeuroVecSource <- function(file_name) {
  ### could put logic in here inspect h5 object.
  new("LatentNeuroVecSource", file_name=file_name)
}



#' Construct a LatentNeuroVec object
#'
#' This function constructs a LatentNeuroVec object from a basis, associated loadings, a NeuroSpace instance, a mask, and an optional offset.
#'
#' @param basis A numeric n-by-k matrix containing the latent vectors forming the reduced space.
#' @param loadings A numeric p-by-k matrix of p loadings.
#' @param space A NeuroSpace instance defining the dimensions, spacing, origin, axes, and transformation of the neuroimaging space.
#' @param mask A 3D logical array, 1D logical vector, or an instance of LogicalNeuroVol class representing the brain mask.
#' @param offset An optional numeric 1-by-p offset vector. If not provided, it defaults to a zero vector.
#'
#' @return A new \code{\linkS4class{LatentNeuroVec}} instance representing the latent neuroimaging vectors.
#'
#' @section Usage:
#' \preformatted{
#' LatentNeuroVec(basis, loadings, space, mask, offset = NULL)
#' }
#'
#' @section Examples:
#' \preformatted{
#' bspace <- NeuroSpace(c(2,2,2,10), c(1,1,1))
#' mask <- array(rnorm(2*2*2) > -100, c(2,2,2))
#' mat <- matrix(rnorm(sum(mask)), 10, sum(mask))
#' pres <- prcomp(mat)
#' svec <- LatentNeuroVec(pres$x, pres$rotation, bspace, mask, offset=colMeans(mat))
#' svec2 <- SparseNeuroVec(mat, bspace, mask)
#' length(indices(svec)) == sum(mask)
#'
#' all.equal(svec2[1:prod(dim(mask))], svec[1:prod(dim(mask))])
#' }
#' @export
#' @rdname LatentNeuroVec-class
LatentNeuroVec <- function(basis, loadings, space, mask, offset=NULL) {
  stopifnot(inherits(space, "NeuroSpace"))

  if (!inherits(mask, "LogicalNeuroVol")) {
    mspace <- NeuroSpace(dim(space)[1:3],
                         spacing(space),
                         origin(space),
                         axes(space),
                         trans(space))
    mask <- LogicalNeuroVol(as.logical(mask), mspace)
  }

  cardinality <- sum(mask)

  if (is.null(offset)) {
    offset <- rep(0, nrow(loadings))
  }

  assert_that(nrow(loadings) == cardinality, msg="`loadings`` must have same number of rows as nonzero entries in `mask`")
  assert_that(ncol(loadings) == ncol(basis), msg="`basis` and `loadings` must have the same number of columns, `k`")
  assert_that(nrow(basis) == dim(space)[4], msg="`basis` must have same number of rows as 4th dimension of `space` argument")

  if (is.matrix(basis)) {
    basis <- Matrix(basis)
  }
  if (is.matrix(loadings)) {
    loadings <- Matrix(loadings)
  }

  new("LatentNeuroVec", basis=basis, loadings=loadings, space=space, mask=mask,
      map=IndexLookupVol(space(mask), as.integer(which(mask))), offset=offset)

}


#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "matrix"),
          def=function (x, i) {
            b1 <- x@basis[as.numeric(i[,1]),,drop=FALSE]
            b2 <- x@loadings[as.numeric(i[,2]),,drop=FALSE]
            rowSums(b1*b2) + x@offset[i[,2]]
          })

#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "integer"),
          def=function (x, i) {
            b1 <- x@basis
            b2 <- x@loadings[as.numeric(i),,drop=FALSE]
            out <- tcrossprod(b1,b2)
            as.matrix(sweep(out, 2, x@offset[i], "+"))
          })


#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "numeric"),
          def=function (x, i) {
            callGeneric(x, as.integer(i))
          })


#'
#' @rdname linear_access-methods
setMethod(f="linear_access", signature=signature(x = "LatentNeuroVec", i = "numeric"),
          def=function (x, i) {
            nels <- prod(dim(x)[1:3])
            n <- ceiling(i/nels)
            offset <- i %% nels
            offset[offset == 0] <- nels

            ll <- lookup(x@map, offset)
            nz <- which(ll > 0)

            if (length(nz) == 0) {
              return(numeric(length(i)))
            }

            idx2d <- cbind(n[nz], ll[nz])

            b1 <- x@basis[idx2d[,1],,drop=FALSE]
            b2 <- x@loadings[idx2d[,2],,drop=FALSE]
            vals <- rowSums(b1*b2) + x@offset[idx2d[,2]]


            ovals <- numeric(length(i))
            ovals[nz] <- vals
            ovals
          })


#' [[
#'
#' @rdname SparseNeuroVec-methods
#' @param x the object
#' @param i the volume index
#' @export
setMethod(f="[[", signature=signature(x="LatentNeuroVec", i="numeric"),
          def = function(x, i) {
            stopifnot(length(i) == 1)
            xs <- space(x)
            dat <- (tcrossprod(x@basis[i,,drop=FALSE], x@loadings))[1,]
            dat <- dat + x@offset
            #dat <- x@data[i,]
            newdim <- dim(xs)[1:3]
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            SparseNeuroVol(dat, bspace, indices=indices(x))
          })




#' @export write_vec
#' @rdname write_vec-methods
#' @param nbit set nbit compression
#' @param compression compression level 1 to 9
#' @param chunk_dim the dimensions of each chunk
setMethod(f="write_vec",signature=signature(x="LatentNeuroVec", file_name="character", format="missing", data_type="missing"),
          def=function(x, file_name, nbit=FALSE, compression=9, chunk_dim=NULL) {
            obj <- to_h5_latentvec(x, file_name, chunk_dim=chunk_dim, nbit=nbit, compression=compression)
            obj$close()
          })


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="LatentNeuroVec", y="LatentNeuroVec"),
          def=function(x,y,...) {
            do.call(NeuroVecSeq, list(x,y,...))
          })


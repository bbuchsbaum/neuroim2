#' @keywords internal
#' @noRd
LatentNeuroVecSource <- function(file_name) {
  ### could put logic in here inspect h5 object.
  new("LatentNeuroVecSource", file_name=file_name)
}



#' Construct a LatentNeuroVec Object
#'
#' @description
#' This function constructs a LatentNeuroVec object, which represents a latent space 
#' representation of neuroimaging data. It combines a basis, associated loadings, 
#' a NeuroSpace instance, a mask, and an optional offset.
#'
#' @param basis A numeric n-by-k matrix containing the latent vectors forming the reduced space.
#' @param loadings A numeric p-by-k matrix of p loadings.
#' @param space A \code{\linkS4class{NeuroSpace}} instance defining the dimensions, spacing, 
#'   origin, axes, and transformation of the neuroimaging space.
#' @param mask A 3D logical array, 1D logical vector, or an instance of 
#'   \code{\linkS4class{LogicalNeuroVol}} class representing the brain mask.
#' @param offset An optional numeric 1-by-p offset vector. If not provided, it defaults to a zero vector.
#'
#' @return A new \code{\linkS4class{LatentNeuroVec}} instance representing the latent neuroimaging vectors.
#'
#' @details
#' The function performs several checks to ensure the consistency of the input data:
#' \itemize{
#'   \item The 'space' argument must be a NeuroSpace object.
#'   \item If 'mask' is not a LogicalNeuroVol, it's converted to one.
#'   \item The number of rows in 'loadings' must match the number of non-zero entries in 'mask'.
#'   \item The number of columns in 'basis' and 'loadings' must be the same.
#'   \item The number of rows in 'basis' must match the 4th dimension of 'space'.
#' }
#' The function also converts 'basis' and 'loadings' to Matrix objects for efficient computation.
#'
#' @examples
#' # Create a simple NeuroSpace
#' bspace <- NeuroSpace(c(2,2,2,10), c(1,1,1))
#' 
#' # Create a mask
#' mask <- array(rnorm(2*2*2) > -100, c(2,2,2))
#' 
#' # Create some random data
#' mat <- matrix(rnorm(sum(mask) * 10), 10, sum(mask))
#' 
#' # Perform PCA
#' pres <- prcomp(mat)
#' 
#' # Create a LatentNeuroVec
#' svec <- LatentNeuroVec(pres$x, pres$rotation, bspace, 
#'                        mask, offset = colMeans(mat))
#' 
#' # Compare with SparseNeuroVec
#' svec2 <- SparseNeuroVec(mat, bspace, mask)
#' 
#' # Check consistency
#' stopifnot(length(indices(svec)) == sum(mask))
#' stopifnot(all.equal(svec2[1:prod(dim(mask))], svec[1:prod(dim(mask))]))
#'
#' @seealso 
#' \code{\link{NeuroSpace-class}}, \code{\link{LogicalNeuroVol-class}}, \code{\link{SparseNeuroVec-class}}
#'
#' @export
#' @rdname LatentNeuroVec-class
LatentNeuroVec <- function(basis, loadings, space, mask, offset = NULL) {
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


#' @noRd
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "matrix"),
          def=function (x, i) {
            b1 <- x@basis[as.numeric(i[,1]),,drop=FALSE]
            b2 <- x@loadings[as.numeric(i[,2]),,drop=FALSE]
            rowSums(b1*b2) + x@offset[i[,2]]
          })

#' @noRd
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "integer"),
          def=function (x, i) {
            b1 <- x@basis
            b2 <- x@loadings[as.numeric(i),,drop=FALSE]
            out <- tcrossprod(b1,b2)
            as.matrix(sweep(out, 2, x@offset[i], "+"))
          })


#' @noRd
setMethod(f="matricized_access", signature=signature(x = "LatentNeuroVec", i = "numeric"),
          def=function (x, i) {
            callGeneric(x, as.integer(i))
          })



#' @noRd
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


#' Extract a Single Volume from LatentNeuroVec
#'
#' @description
#' This method extracts a single volume from a LatentNeuroVec object and returns it as a SparseNeuroVol.
#'
#' @param x A LatentNeuroVec object.
#' @param i A numeric index specifying which volume to extract.
#'
#' @return A \code{\linkS4class{SparseNeuroVol}} object representing the extracted volume.
#'
#' @details
#' The method performs the following steps:
#' 1. Checks that the input index is a single number.
#' 2. Computes the data for the specified volume using the basis and loadings.
#' 3. Adds the offset to the computed data.
#' 4. Creates a new NeuroSpace object for the extracted volume.
#' 5. Returns a SparseNeuroVol object with the computed data and new space.
#'
#' @examples
#' # Assuming 'svec' is a previously created LatentNeuroVec object
#' \dontrun{
#' single_vol <- svec[[1]]  # Extract the first volume
#' }
#'
#' @seealso 
#' \code{\link{LatentNeuroVec-class}}, \code{\link{SparseNeuroVol-class}}
#'
#' @rdname SparseNeuroVec-methods
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




#' Write LatentNeuroVec to File
#'
#' @description
#' This method writes a LatentNeuroVec object to an HDF5 file.
#'
#' @param x A LatentNeuroVec object to be written.
#' @param file_name A character string specifying the output file name.
#' @param nbit Logical; whether to use N-bit compression. Default is FALSE.
#' @param compression An integer from 1 to 9 specifying the compression level. Default is 9.
#' @param chunk_dim An optional vector specifying the chunk dimensions for HDF5 storage.
#'
#' @return None (called for side effects).
#'
#' @details
#' This method uses the internal \code{to_h5_latentvec} function to write the LatentNeuroVec 
#' object to an HDF5 file. It allows for customization of compression settings and chunking.
#'
#' @examples
#' \dontrun{
#' # Assuming 'svec' is a previously created LatentNeuroVec object
#' write_vec(svec, "output.h5", compression = 6)
#' }
#'
#' @seealso 
#' \code{\link{LatentNeuroVec-class}}, \code{\link{read_vec}}
#'
#' @rdname write_vec-methods
#' @export
setMethod(f="write_vec",
          signature=signature(x="LatentNeuroVec", file_name="character",
                              format="missing", data_type="missing"),
          def=function(x, file_name, nbit=FALSE, compression=9, chunk_dim=NULL) {
            obj <- to_h5_latentvec(x, file_name, chunk_dim=chunk_dim,
                                   nbit=nbit, compression=compression)
            obj$close()
          })


#' Concatenate LatentNeuroVec Objects
#'
#' @description
#' This method concatenates two or more LatentNeuroVec objects.
#'
#' @param x A LatentNeuroVec object.
#' @param y Another LatentNeuroVec object.
#' @param ... Additional LatentNeuroVec objects to concatenate.
#'
#' @return A \code{\linkS4class{NeuroVecSeq}} object representing the concatenated vectors.
#'
#' @details
#' This method creates a NeuroVecSeq object from the input LatentNeuroVec objects, 
#' effectively concatenating them along the 4th dimension.
#'
#' @examples
#' \dontrun{
#' # Assuming 'svec1' and 'svec2' are previously created LatentNeuroVec objects
#' concatenated <- concat(svec1, svec2)
#' }
#'
#' @seealso 
#' \code{\link{LatentNeuroVec-class}}, \code{\link{NeuroVecSeq-class}}
#'
#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="LatentNeuroVec", y="LatentNeuroVec"),
          def=function(x,y,...) {
            do.call(NeuroVecSeq, list(x,y,...))
          })

#' resample an image to match the space of another image
#'
#' @export
#' @param the source image
#' @param the target image
#' @param ... extra args
setGeneric("resample", function(source, target, ...) standardGeneric("resample"))


#' @export
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' @export
setGeneric("drop", function(x) standardGeneric("drop"))

#' @export
setGeneric("as.matrix", function(x) standardGeneric("as.matrix"))


setGeneric("scale")



#' print an object
#'
#' @param x the object to print
#' @param ... additional arguments
#' @export
#' @rdname print-methods
setGeneric(name="print", def=function(x, ...) standardGeneric("print"))

#' extract data values of object
#'
#' @param x the object to get values from
#' @param ... additional arguments
#' @export
#' @rdname values-methods
setGeneric(name="values", def=function(x, ...) standardGeneric("values"))

#' extract values from an array-like object using linear indexing
#'
#' @param x a data source
#' @param i a vector of indices
#' @param ... additional arguments
#' @export
#' @rdname linear_access-methods
setGeneric(name="linear_access", def=function(x, i, ...) standardGeneric("linear_access"))


#' extract values from a matricized (x,y,z) of a 4D tensor using a space-time coordinate matrix
#'
#' @param x a data source
#' @param i an index matrix
#' @param ... additional arguments
#' @rdname matricized_access-methods
setGeneric(name="matricized_access", def=function(x, i, ...) standardGeneric("matricized_access"))



#' load data from a data source
#'
#' @param x a data source
#' @param ... additional arguments
#' @export
#' @rdname load_data-methods
setGeneric(name="load_data", def=function(x, ...) standardGeneric("load_data"))

#' apply a function to an object
#'
#' @param x the object that is mapped
#' @param m the mapping object
#' @param ... additional arguments
#' @export
#' @rdname map-methods
setGeneric(name="mapf", def=function(x, m, ...) standardGeneric("mapf"))

#' extract an ordered series of 3d volumes
#'
#' @param x the object that supplies the volume data
#' @param indices the subset of volumes to extract
#' @param ... additional arguments
#' @export
#' @rdname vols-methods
#' @examples
#'
#' vec <- read_vec(system.file("extdata", "global_mask_v5.nii", package="neuroim2"))
#' vs <- vols(vec)
#' length(vs) == dim(vec)[4]
#'
#' vs <- vols(vec, indices=1:3)
#' length(vs) == 3
setGeneric(name="vols", def=function(x, indices, ...) standardGeneric("vols"))

#' extract an ordered list of 1d vectors
#'
#' @param x the object that supplies the vector data
#' @param subset the subset of vectors to extract
#' @param ... additional arguments
#' @export
#' @rdname vectors-methods
setGeneric(name="vectors", def=function(x, subset, ...) standardGeneric("vectors"))

#' cut a vector-valued object into a list of sub-blocks
#'
#' @param x the object to split
#' @param indices a vector of indices for the blocks. Must match the length of the inut vector.
#' @param ... additional arguments
#' @export
#' @rdname split_blocks-methods
setGeneric(name="split_blocks", def=function(x, indices, ...) standardGeneric("split_blocks"))



#' partition an image into a set of disjoint clusters
#'
#' @param x the object to partition
#' @param k the number of clusters
#' @param ... additional arguments
#' @export
#' @rdname partition-methods
setGeneric(name="partition", def=function(x, k, ...) standardGeneric("partition"))



#' cut an object into a list of spatial or spatiotemporal clusters
#'
#' @param x the object to split
#' @param clusters a vector of cluster indices to split by
#' @param ... additional arguments
#' @export
#' @rdname split_clusters-methods
setGeneric(name="split_clusters", def=function(x, clusters, ...) standardGeneric("split_clusters"))


#' extract an ordered series of 2d slices
#'
#' @param x the object that supplies the slices
#' @param ... additional arguments
#' @export
#' @rdname slices-methods
setGeneric(name="slices", def=function(x, ...) standardGeneric("slices"))




#' extract the number of dimensions of an object
#'
#' @param x n-dimensional object
#' @param ... additional arguments
#' @export
#' @examples
#'
#' x = NeuroSpace(c(10,10,10), spacing=c(1,1,1))
#' ndim(x) == 3
#' x = NeuroSpace(c(10,10,10,3), spacing=c(1,1,1))
#' ndim(x) == 4
#'
#' @rdname ndim-methods
setGeneric(name="ndim", def=function(x, ...) standardGeneric("ndim"))


#' get dimensions of an axis
#'
#' @param x the object
#' @param axis the axis to return the dimension of
#' @rdname dim_of-methods
setGeneric(name="dim_of", def=function(x, axis) standardGeneric("dim_of"))

#' find dimensions of a given axis
#'
#' @param x the object
#' @param axis the axis to return dimension of
#' @rdname which_dim-methods
setGeneric(name="which_dim", def=function(x, axis) standardGeneric("which_dim"))


#' add a dimension to an object
#'
#' @param x a dimensioned object
#' @param n the size of the dimension to add
#' @export
#' @rdname add_dim-methods
#' @examples
#' x = NeuroSpace(c(10,10,10), c(1,1,1))
#' x1 <- add_dim(x, 10)
#' ndim(x1) == 4
#' dim(x1)[4] == 10
setGeneric(name="add_dim", def=function(x, n) standardGeneric("add_dim"))

#' drop a dimension from an object
#'
#' @param x a dimensioned object
#' @param dimnum the index of the dimension to drop
#' @export
#' @rdname drop_dim-methods
#' @examples
#' x = NeuroSpace(c(10,10,10), c(1,1,1))
#' x1 <- drop_dim(x)
#' ndim(x1) == 2
#' dim(x1)[2] == 10
setGeneric(name="drop_dim", def=function(x, dimnum) standardGeneric("drop_dim"))

#' extract geometric properties of an image.
#'
#' @param x the object to query, e.g. an instance of \code{\linkS4class{NeuroVol}} or \code{\linkS4class{NeuroVec}}
#' @param ... additional arguments
#' @return an object representing the geometric space of the image of type \code{\linkS4class{NeuroSpace}}
#' @export space
#' @examples
#' x = NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(rnorm(10*10*10), x)
#' identical(x,space(vol))
#'
#' @rdname space-methods
setGeneric(name="space", def=function(x, ...) standardGeneric("space"))

#' fill disjoint sets of values with the output of a function
#'
#' @param x the object to split
#' @param fac the \code{factor} to split by
#' @param FUN the function to summarize the the sets
#' @return a new object where the original values have been replaced by the function output
#' @export
#' @rdname split_fill-methods
#' @details \code{FUN} can either return a scalar for each input vector or a vector equal to the length of the input vector.
#' If it returns a scalar then every voxel in the set will be filled with that value in the output vector.
#' @examples
#'
#' ## summarize with mean -- FUN returns a scalar
#' x = NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(rnorm(10*10*10), x)
#' fac <- factor(rep(1:10, length.out=1000))
#' ovol.mean <- split_fill(vol, fac, mean)
#' identical(dim(ovol.mean), dim(vol))
#' length(unique(as.vector(ovol.mean))) == 10
#' ## transform by reversing vector -- FUN returns a vector.
#' ovol2 <- split_fill(vol, fac, rev)
setGeneric(name="split_fill", def=function(x, fac, FUN) standardGeneric("split_fill"))

#' map values from one set of values to another set using a user-supplied lookup table
#'
#' @param x the object to map values from
#' @param lookup the lookup table. The first column is the "key" the second column is the "value".
#' @return a new object where the original values have been filled in with the values in the lookup table
#' @export
#' @examples
#' x <- NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(sample(1:10, 10*10*10, replace=TRUE), x)
#'
#' ## lookup table is list
#' lookup <- lapply(1:10, function(i) i*10)
#' names(lookup) <- 1:10
#' ovol <- map_values(vol, lookup)
#'
#' ## lookup table is matrix. First column is key, second column is value
#' names(lookup) <- 1:length(lookup)
#' lookup.mat <- cbind(as.numeric(names(lookup)), unlist(lookup))
#' ovol2 <- map_values(vol, lookup.mat)
#' all.equal(as.vector(ovol2), as.vector(ovol))
#'
#' @rdname map_values-methods
setGeneric(name="map_values", def=function(x, lookup) standardGeneric("map_values"))



#' center/scale row-subsets of a matrix or matrix-like object
#'
#' @param x a numeric matrix or matrix-like object
#' @param f the splitting object, typically a \code{factor} or set of \code{integer} indices. must be equal to number of rows of matrix.
#' @param center should values within each submatrix be centered? (mean removed from each column of submatrix)
#' @param scale should values be scaled? (divide vector by standard deviation from each column of submatrix)
#' @return a new matrix or matrix-like object where the original rows have been grouped by \code{f} and then centered and/or scaled for each grouping
#' @docType methods
#' @export
#' @examples
#'
#' M <- matrix(rnorm(1000), 10, 100)
#' fac <- factor(rep(1:2, each=5))
#' Ms <- split_scale(M, fac)
#'
#' ## correctly centered
#' all(abs(apply(Ms[fac == 1,], 2, mean)) < .000001)
#' all(abs(apply(Ms[fac == 2,], 2, mean)) < .000001)
#'
#' # correctly scaled
#' all.equal(apply(Ms[fac == 1,], 2, sd), rep(1, ncol(Ms)))
#' all.equal(apply(Ms[fac == 2,], 2, sd), rep(1, ncol(Ms)))
#' @rdname split_scale-methods
setGeneric(name="split_scale", def=function(x, f, center, scale) standardGeneric("split_scale"))

#' summarize subsets of an object by first splitting by row and then "reducing" by a summary \code{function}
#'
#' @param x a numeric matrix(like) object
#' @param fac the factor to define subsets of the object
#' @param FUN the function to apply to each subset. if \code{FUN} is missing, than the mean of each sub-matrix column is computed.
#' @return a new \code{matrix} where the original values have been "reduced" by the supplied function.
#' @docType methods
#' @details if \code{FUN} is supplied it must take a vector and return a single scalar value. If it returns more than one value, an error will occur.
#'
#' if \code{x} is a \code{\linkS4class{NeuroVec}} instance then voxels (dims 1:3) are treated as columns and time-series (dim 4) as rows.
#' The summary function then is applied to groups of voxels. However, if the goal is to apply a function to groups of time-points,
#' then this can be achieved as follows:
#'
#' \code{ split_reduce(t(as.matrix(bvec)), fac) }
#'
#'
#' @export
#' @examples
#' mat = matrix(rnorm(100*100), 100, 100)
#' fac = factor(sample(1:3, nrow(mat), replace=TRUE))
#' ## compute column means of each sub-matrix
#' ms <- split_reduce(mat, fac)
#' all.equal(row.names(ms), levels(fac))
#'
#' ## compute column medians of each sub-matrix
#' ms <- split_reduce(mat, fac, median)
#'
#' ## compute time-series means grouped over voxels.
#' ## Here, \code{length(fac)} must equal the number of voxels: \code{prod(dim(bvec)[1:3]}
#' bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
#' fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
#' ms <- split_reduce(bvec, fac)
#' ms2 <- split_reduce(bvec, fac, mean)
#' all.equal(row.names(ms), levels(fac))
#' all.equal(ms,ms2)
#'
#' @rdname split_reduce-methods
setGeneric(name="split_reduce", def=function(x, fac, FUN) standardGeneric("split_reduce"))


#' extract the voxel dimensions of an image
#'
#' @param x the object
#' @return a numeric vector
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' all.equal(spacing(bspace), c(2,2,2))
#' @rdname spacing-methods
setGeneric(name="spacing", def=function(x) standardGeneric("spacing"))

#' extract the spatial bounds (origin + dim * spacing) of an image
#' param x the object

#' @param x the object with \code{bounds} property
#' @return a \code{matrix} where each row contains the min (column 1) and max (column 2) bounds of the image dimension from 1 to \code{ndim(image)}.
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' b <- bounds(bspace)
#' nrow(b) == ndim(bspace)
#' ncol(b) == 2
#' @rdname bounds-methods
setGeneric(name="bounds",     def=function(x) standardGeneric("bounds"))


#' extract image axes
#'
#' @param x an object with a set of axes
#' @export
#' @rdname axes-methods
setGeneric(name="axes",  def=function(x) standardGeneric("axes"))

#' extract image origin
#'
#' @param x an object with an origin
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' origin(bspace)
#'
#' @rdname origin-methods
setGeneric(name="origin", def=function(x) standardGeneric("origin"))


#' return the centroid of an object
#'
#' @param x an object with a centroid
#' @param ... extra args
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' centroid(bspace)
#'
#' @rdname centroid-methods
setGeneric(name="centroid", def=function(x, ...) standardGeneric("centroid"))

#' return a matrix of centroids of an object
#'
#' @param x an object with multiple centroids (e.g. a \code{ClusteredNeuroVol})
#' @param ... extra args
#' @export
#' @rdname centroids-methods
setGeneric(name="centroids", def=function(x, ...) standardGeneric("centroids"))



#' extract image coordinate transformation
#'
#' @param x an object with a transformation
#' @export
#' @details
#' This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
#' @details This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
#' see \code{\linkS4class{NeuroSpace}}
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' trans(bspace)
#' all.equal(dim(trans(bspace)), c(4,4))
#' @rdname trans-methods
setGeneric(name="trans",  def=function(x) standardGeneric("trans"))

#' extract inverse image coordinate transformation
#'
#' @param x an object
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' itrans <- inverse_trans(bspace)
#' identical(trans(bspace) %*% inverse_trans(bspace), diag(4))
#' @rdname inverse_trans-methods
setGeneric(name="inverse_trans", def=function(x) standardGeneric("inverse_trans"))

#' read a sequence of elements from an input source
#'
#' @param x the input channel
#' @param num_elements the number of elements to read
#' @return the elements as a vector
#' @export
#' @rdname read_elements-methods
setGeneric(name="read_elements", def=function(x, num_elements) standardGeneric("read_elements"))

#' read a set of column vector from an input source (e.g. \code{ColumnReader})
#' @param x the input channel
#' @param column_indices the column indices
#' @return a \code{matrix} consisting of the requested column vectors
#' @export
#' @rdname read_columns-methods
setGeneric(name="read_columns", def=function(x, column_indices) standardGeneric("read_columns"))


#' write a sequence of elements from an input source
#'
#' @param x the output channel
#' @param els the elements to write
#' @export
#' @rdname write_elements-methods
setGeneric(name="write_elements", def=function(x, els) standardGeneric("write_elements"))


#' write a 3d image volume to disk
#'
#' @param x an image object, typically a \code{\linkS4class{NeuroVol}} instance.
#' @param file_name output file name
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param data_type output data type, If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
#' Otherwise output format will be inferred from R the datatype of the image.
#' @export
#' @details
#'
#'  The output format will be inferred from file extension.
#' @details The output format will be inferred from file extension.
#'  \code{write_vol(x, "out.nii")} outputs a NIFTI file.
#'  \code{write_vol(x, "out.nii.gz")} outputs a gzipped NIFTI file.
#'
#' No other file output formats are currently supported.
#'
#'
#' @examples
#'
#' bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#' \dontrun{
#' write_vol(bvol, "out.nii")
#' write_vol(bvol, "out.nii.gz")
#' }
#' @rdname write_vol-methods
setGeneric(name="write_vol",  def=function(x, file_name, format, data_type) standardGeneric("write_vol"))


#' write a 4d image vector to disk
#'
#' @param x an image object, typically a \code{NeuroVec} instance.
#' @param file_name output file name.
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param data_type the numeric data type. If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
#' Otherwise output format will be inferred from R the datatype of the image.
#' @param ... extra args
#' @export
#' @examples
#'
#' bvec <- NeuroVec(array(0, c(10,10,10,10)), NeuroSpace(c(10,10,10,10), c(1,1,1)))
#' \dontrun{
#' writeVector(bvol, "out.nii")
#' writeVector(bvol, "out.nii.gz")
#' writeVector(bvec, "out.nii")
#' writeVector(bvec, "out.nii.gz")
#' }
#' @rdname write_vec-methods
setGeneric(name="write_vec",  def=function(x, file_name, format, data_type, ...) standardGeneric("write_vec"))



#' remap the grid-to-world coordinates mapping of an image.
#'
#' @param x the object
#' @param orient the orientation code indcating the "remapped" axes.
#' @return a reoriented space
#' @rdname reorient-methods
#' @export
setGeneric(name="reorient", def=function(x, orient) standardGeneric("reorient"))


#' convert 1d indices to n-dimensional grid coordinates
#' @param x the object
#' @param idx the 1d \code{vector} of indices
#' @return a matrix of grid coordinates
#' @export
#' @examples
#'
#'  bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#'  idx <- 1:10
#'  g <- index_to_grid(bvol, idx)
#'  bvol[g]
#'
#' @rdname index_to_grid-methods
setGeneric(name="index_to_grid",   def=function(x, idx) standardGeneric("index_to_grid"))

#' convert 1d indices to n-dimensional real world coordinates
#'
#' @param x the object
#' @param idx the 1D indices
#' @return a matrix of real coordinates
#' @export
#' @examples
#' bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#' idx <- 1:10
#' g <- index_to_coord(bvol, idx)
#' idx2 <- coord_to_index(bvol, g)
#' all.equal(idx, idx2)
#' @rdname index_to_coord-methods
setGeneric(name="index_to_coord",   def=function(x, idx) standardGeneric("index_to_coord"))

#' convert n-dimensional real world coordinates to 1D indices
#'
#' @param x the object
#' @param coords a matrix of real world coordinates
#' @return a vector of indices
#' @export
#' @rdname coord_to_index-methods
setGeneric(name="coord_to_index",   def=function(x, coords) standardGeneric("coord_to_index"))

#' convert n-dimensional real world coordinates to grid coordinates
#'
#' @param x the object
#' @param coords a matrix of real world coordinates
#' @return a matrix of grid coordinates
#' @export
#' @rdname coord_to_grid-methods
setGeneric(name="coord_to_grid",   def=function(x, coords) standardGeneric("coord_to_grid"))

#' Generic function to convert N-dimensional grid coordinate coordinates to real world coordinates
#'
#' @param x the object
#' @param coords a matrix of grid coordinates
#' @return a matrix of real coordinates
#' @export
#' @rdname grid_to_coord-methods
setGeneric(name="grid_to_coord",   def=function(x, coords) standardGeneric("grid_to_coord"))

#' Generic function to convert voxel coordinates in the reference space (LPI) to native array space.
#'
#' @param x the object
#' @param vox a matrix of LPI voxel coordinates
#' @return a matrix of native voxel coordinates
#' @export
#' @rdname grid_to_grid-methods
setGeneric(name="grid_to_grid",   def=function(x, vox) standardGeneric("grid_to_grid"))


#' Generic function to convert N-dimensional grid coordinate to 1D indices
#' @param x the object, typically a \code{NeuroVol} or \code{NeuroSpace} instance.
#' @param coords a matrix where each row is a coordinate or a vector of length equal to \code{ndim(x)}
#' @return a vector of indices
#' @export
#' @rdname grid_to_index-methods
setGeneric(name="grid_to_index",   def=function(x, coords) standardGeneric("grid_to_index"))



#' Generic function to extract a sub-vector from a \code{NeuroVec} object.
#' @param x four-dimensional image
#' @param i the indices of the volume(s) to extract
#' @param ... additional arguments
#' @return a  \code{NeuroVec} object that is a sub-sequence of the supplied object.
#' @export
#' @examples
#' bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
#' vec <- sub_vector(bvec,1:2)
#' all.equal(2, dim(vec)[4])
#'
#' vec <- sub_vector(bvec, c(1,3,5,7))
#' all.equal(4, dim(vec)[4])
#'
#' mask <- LogicalNeuroVol(rep(TRUE, 24*24*24), NeuroSpace(c(24,24,24), c(1,1,1)))
#' svec <- SparseNeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)),
#' NeuroSpace(c(24,24,24,24), c(1,1,1)), mask)
#' vec <- sub_vector(svec, c(1,3,5))
#' all.equal(3, dim(vec)[4])
#' @rdname sub_vector-methods
setGeneric(name="sub_vector", def=function(x, i, ...) standardGeneric("sub_vector"))


#' Generic functions to scale (center and/or normalize by standard deviation) each series of a 4D image
#' That is, if the 4th dimension is 'time' each series is a 1D time series.
#' @param x a four dimensional image
#' @param center a \code{logical} value indicating whether series should be centered
#' @param scale a \code{logical} value indicating whether series should be divided by standard deviation
#' @export
#' @examples
#' bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
#' res <- scale_series(bvec, TRUE, TRUE)
#' @rdname scale_series-methods
setGeneric(name="scale_series", def=function(x, center, scale) standardGeneric("scale_series"))


#' Convert to from dense to sparse representation
#'
#' @param x the object to make sparse, e.g. \code{DenseNeuroVol} or \code{DenseNeuroVec}
#' @param mask the elements to retain
#' @param ... additional arguments
#'
#' @details \code{mask} can be an integer vector of 1D indices or a mask volume of class \code{LogicalNeuroVol}
#' @export
#' @examples
#' bvol <- NeuroVol(array(runif(24*24*24), c(24,24,24)), NeuroSpace(c(24,24,24), c(1,1,1)))
#' indmask <- sort(sample(1:(24*24*24), 100))
#' svol <- as.sparse(bvol, indmask)
#'
#'
#' mask <- LogicalNeuroVol(runif(length(indmask)), space=space(bvol), indices=indmask)
#' sum(mask) == 100
#' @rdname as.sparse-methods
setGeneric(name="as.sparse", def=function(x, mask, ...) standardGeneric("as.sparse"))

#' Convert to dense representation
#'
#' @param x the object to densify
#' @rdname as.dense-methods
setGeneric(name="as.dense", def=function(x) standardGeneric("as.dense"))

#' Convert to a LogicalNeuroVol
#' @param x the object to binarize
#' @param indices the indices to set to TRUE
#' @export
#' @rdname as.mask-methods
setGeneric(name="as.mask", def=function(x, indices) standardGeneric("as.mask"))


#' Extract set of patches
#'
#' generate a set of coordinate "patches" of fixed size from an image object.
#'
#' @param x the object to extract patches from
#' @param dims a vector indicating the dimensions of the patches
#' @param mask mask indicating the valid patch area
#' @param ... additional args
#' @rdname patch_set-methods
setGeneric(name="patch_set", def=function(x, dims, mask, ...) standardGeneric("patch_set"))

#' num_clusters
#'
#' @param x the object to extract number of clusters
#' @export
#' @rdname num_clusters-methods
setGeneric(name="num_clusters", def=function(x) standardGeneric("num_clusters"))


#' Extract coordinates
#'
#' @param x the object to extract coordinates from
#' @param ... additional arguments
#' @export
#' @rdname coords-methods
setGeneric(name="coords", def=function(x, ...) standardGeneric("coords"))

#' Extract indices
#' @param x the object to extract indices
#' @export
#' @rdname indices-methods
setGeneric(name="indices", def=function(x) standardGeneric("indices"))

#' Index Lookup operation
#' @param x the object to query
#' @param i the index to lookup
#' @param ... additional arguments
#' @export
#' @rdname lookup-methods
setGeneric(name="lookup", def=function(x, i, ...) standardGeneric("lookup"))


#' Extract one or more series from object and return as ROI object
#'
#' @inheritParams series
#' @rdname series-methods
setGeneric(name="series_roi", def=function(x, i, ...) standardGeneric("series_roi"))

#' Extract one or more series from object
#' @param x the object
#' @param i the series indices
#' @param ... additional arguments
#' @export
#' @rdname series-methods
setGeneric(name="series", def=function(x, i, ...) standardGeneric("series"))

#' slice
#'
#' Extract a 2D slice from an image volume
#'
#' @param x the object
#' @param zlevel coordinate (in voxel units) along the sliced axis
#' @param along the axis along which to slice
#' @param orientation the target orientation of the 2D slice
#' @param ... additional arguments
#' @export
#' @rdname slice-methods
setGeneric(name="slice", def=function(x, zlevel, along, orientation, ...) standardGeneric("slice"))

#' Render an image to create a drawable image.
#' @param x the object, e.g. an instance of type \code{NeuroSlice}
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @rdname render-methods
setGeneric(name="render", def=function(x, width, height, colmap,...) standardGeneric("render"))


#' Render a slice at z coordinate
#' @param x the object, e.g. an instance of type \code{Layer} or \code{Overlay}
#' @param zpos the z coordinate to slice through.
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @rdname render_slice-methods
setGeneric(name="render_slice", def=function(x, zpos, width, height, colmap,...) standardGeneric("render_slice"))




#' Extract permutation matrix
#' @param x the object
#' @param ... additional arguments
#' @export
#' @rdname perm_mat-methods
setGeneric(name="perm_mat", def=function(x, ...) standardGeneric("perm_mat"))

#' Concatenate two objects in the time dimension
#'
#' @param x the first object, typically \code{NeuroVol} or \code{NeuroVec}
#' @param y the second object, typically \code{NeuroVol} or \code{NeuroVec}
#' @details The \code{x} and \code{y} images must have compatible dimensions. a \code{NeuroVol} can be concatenated to \code{NeuroVec}, and vice versa. See examples.
#' @param ... additional objects
#'
#' @examples
#' bv1 <- NeuroVol(rep(1,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
#' bv2 <- NeuroVol(rep(2,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
#' bv3 <- concat(bv1,bv2)
#' inherits(bv3, "NeuroVec")
#'
#' bv4 <- concat(bv3, bv1)
#' dim(bv4)[4] == 3
#' bv5 <- concat(bv1, bv3)
#' dim(bv4)[4] == 3
#'
#' bv6 <- concat(bv4,bv5)
#' dim(bv6)[4] == 6
#'
#' @export
#' @rdname concat-methods
setGeneric(name="concat", def=function(x,y, ...) standardGeneric("concat"))


#' Find connected components
#' @name conn_comp
#' @param x the image object
#' @param ... additonal arguments
#' @export
#' @rdname conn_comp-methods
setGeneric(name="conn_comp", def=function(x, ...) standardGeneric("conn_comp"))

#' extract voxel coordinates
#'
#' @param x the object to extract voxels from
#' @param ... additional arguments to function
#' @export
#' @rdname voxels-methods
setGeneric(name="voxels", def=function(x, ...) standardGeneric("voxels"))


if (!isGeneric("image"))
  setGeneric("image", function(x, ...) standardGeneric("image"))

if (!isGeneric("as.raster"))
  setGeneric("as.raster", function(x, ...) standardGeneric("as.raster"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{FileFormat}} instance.
#' Will test for match to either header file or data file
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @rdname file_matches-methods
#' @export
setGeneric(name="file_matches", def=function(x, file_name) standardGeneric("file_matches"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{FileFormat}} instance.
#' Will test for match to header file only
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @export
#' @rdname header_file_matches-methods
setGeneric(name="header_file_matches", def=function(x, file_name) standardGeneric("header_file_matches"))

#' Generic function to test whether a file name conforms to the given a \code{\linkS4class{FileFormat}} instance.
#' Will test for match to data file only
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @export
#' @rdname data_file_matches-methods
setGeneric(name="data_file_matches", def=function(x, file_name) standardGeneric("data_file_matches"))

#' Generic function to get the name of the header file, given a file name and a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return the correct header name
#' @rdname header_file-methods
#' @export
setGeneric(name="header_file", def=function(x, file_name) standardGeneric("header_file"))

#' Generic function to get the name of the data file, given a file name and a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return the correct header name
#' @export
#' @rdname data_file-methods
setGeneric(name="data_file", def=function(x, file_name) standardGeneric("data_file"))

#' Generic function to strip extension from file name, given a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return file_name without extension
#' @export
#' @rdname strip_extension-methods
setGeneric(name="strip_extension", def=function(x, file_name) standardGeneric("strip_extension"))

#' Generic function to read image meta info given a file and a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name contianing meta information
#' @export
#' @rdname read_meta_info-methods
setGeneric(name="read_meta_info", def=function(x, file_name) standardGeneric("read_meta_info"))


#' Generic function to position kernel in a postion in image space
#' @param x the kernel object
#' @param sp the space to embed the kernel
#' @param center_voxel the voxel marking the center of the kernel in the embedded space
#' @export
#' @rdname embed_kernel-methods
setGeneric("embed_kernel", def=function(x, sp, center_voxel, ...) standardGeneric("embed_kernel"))



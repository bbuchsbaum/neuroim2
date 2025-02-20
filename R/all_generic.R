#' @export
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' Generic Drop Method
#'
#' Provides a mechanism to remove dimensions or elements from an object.
#'
#' @param x An object.
#' @return A modified version of the input object with reduced dimensions or elements.
#'
#' @export
setGeneric("drop", function(x) standardGeneric("drop"))

#' Generic as.matrix Method
#'
#' Coerces an object to a matrix.
#'
#' @param x An object to be coerced to a matrix.
#' @param ... Additional arguments passed to methods.
#' @return A matrix representation of the input object.
#'
#' @export
setGeneric("as.matrix", function(x, ...) standardGeneric("as.matrix"))

#' Generic Scale Method
#'
#' Scales an object by (typically) subtracting the mean and dividing by the standard deviation.
#'
#' @param x The object to be scaled.
#' @param ... Additional arguments for scaling methods.
#' @return A scaled version of the input object.
#'
#' @export
setGeneric("scale", function(x, ...) standardGeneric("scale"))

#' Resample an Image to Match the Space of Another Image
#'
#' This function resamples a source image to match the spatial properties (dimensions, resolution, and orientation) of a target image.
#'
#' @param source An object representing the source image to be resampled. This could be a 3D or 4D image object, depending on the use case.
#' @param target An object representing the target image, whose spatial properties will be used as the reference for resampling the source image.
#' @param ... Additional arguments passed to the resampling function, such as interpolation method, boundary handling, or other resampling options.
#'
#' @return An object representing the resampled source image, with the same spatial properties as the target image.
#'
#' @examples
#'
#' img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' rspace <- space(img)
#'
#' ### normally, one would resample from two existing soource and target spaces.
#' ### But here we manually create the target space, which is a bit ugly.
#'
#' newtrans4X3 <- trans(img)[1:4, 1:3]
#' newtrans4X3 <- newtrans4X3 * c(.5,.5,.5,1)
#' newtrans <- cbind(newtrans4X3, c(space(img)@origin,1))
#'
#' rspace <- NeuroSpace(rspace@dim*2, rspace@spacing/2, origin=rspace@origin, trans=trans(img))
#' \donttest{
#' rvol <- resample(img, rspace)
#' }
#'
#'
#' @export
#' @rdname resample-methods
setGeneric("resample", function(source, target, ...) standardGeneric("resample"))


#' print an object
#'
#' @param x the object to print
#' @param ... additional arguments
#' @keywords internal
#' @noRd
setGeneric(name="print_", def=function(x, ...) standardGeneric("print_"))

#' Extract Data Values of an Object
#'
#' @param x the object to get values from
#' @param ... additional arguments
#' @return A vector or array containing the values extracted from the object.
#' @export
#' @rdname values-methods
setGeneric(name="values", def=function(x, ...) standardGeneric("values"))

#' Extract values from an array-like object using linear indexing.
#'
#' This function extracts the values of the elements in an array-like object using
#' linear indexing. Linear indexing is a way of indexing an array by a single index
#' that is computed from multiple indices using a formula.
#'
#' @param x a data source.
#' @param i a vector of indices.
#' @param ... additional arguments to be passed to methods.
#' @return A vector containing the values at the specified linear indices.
#' @export
setGeneric(name="linear_access", def=function(x, i, ...) standardGeneric("linear_access"))


#' Extract values from a matricized (x,y,z) of a 4D tensor using a space-time coordinate matrix.
#'
#' This function extracts the values of the elements in a 4D tensor, where the first three dimensions
#' (x,y,z) have been matricized into a single dimension, using a space-time coordinate matrix.
#' The input \code{i} is an index matrix, where each row specifies a (x,y,z,t) coordinate.
#' The output is a vector containing the values of the elements in \code{x} at the specified
#' space-time coordinates.
#'
#' @param x a data source.
#' @param i an index matrix specifying the space-time coordinates.
#' @param ... additional arguments to be passed to methods.
#' @return A vector containing the values at the specified space-time coordinates.
#' @rdname matricized_access-methods
#' @export
setGeneric(name="matricized_access", def=function(x, i, ...) standardGeneric("matricized_access"))



#' Load data from a data source.
#'
#' This function loads data from a data source and returns it in a format that is compatible with
#' other functions in the neuroim2 package. The format of the returned data depends on the type
#' of data source used.
#'
#' @param x a data source.
#' @param ... additional arguments to be passed to methods.
#' @return A data object in a format compatible with neuroim2 package functions.
#' @keywords internal
#' @export
setGeneric(name="load_data", def=function(x, ...) standardGeneric("load_data"))

#' Apply a function to an object.
#'
#' This function applies a function to an object, with additional arguments passed to the function
#' using the \code{...} argument. The mapping object specifies how the function is to be applied,
#' and can take many different forms, depending on the object and function used. The return value
#' depends on the function used.
#'
#' @param x the object that is mapped.
#' @param m the mapping object.
#' @param ... additional arguments to be passed to the function.
#' @return The result of applying the mapping function to the input object.
#' @export
#' @rdname map-methods
setGeneric(name="mapf", def=function(x, m, ...) standardGeneric("mapf"))

#' Extract an ordered series of 3D volumes.
#'
#' This function extracts an ordered series of 3D volumes from an object that supplies volume data.
#' The \code{indices} argument specifies the subset of volumes to extract, and can be a vector
#' of indices or a logical vector. The return value is a list containing the extracted volumes
#' in the same order as the specified indices.
#'
#' @param x the object that supplies the volume data.
#' @param indices the subset of volumes to extract.
#' @param ... additional arguments to be passed to methods.
#' @export
#' @rdname vols-methods
#' @examples
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' vs <- vols(vec)
#' length(vs) == dim(vec)[4]
#'
#' vs <- vols(vec, indices=1:3)
#' length(vs) == 3
setGeneric(name="vols", def=function(x, indices, ...) standardGeneric("vols"))


#' Extract an ordered list of 1D vectors.
#'
#' This function extracts an ordered list of 1D vectors from an object that supplies vector data.
#' The \code{subset} argument specifies the subset of vectors to extract, and can be a vector
#' of indices or a logical vector. The return value is a list containing the extracted vectors
#' in the same order as the specified indices.
#'
#' @param x the object that supplies the vector data.
#' @param subset the subset of vectors to extract.
#' @param ... additional arguments to be passed to methods.
#' @export
#' @rdname vectors-methods
setGeneric(name="vectors", def=function(x, subset, ...) standardGeneric("vectors"))

#' Cut a vector-valued object into a list of sub-blocks
#'
#' Splits a vector-valued object into a list of sub-blocks defined by a vector of indices.
#'
#' @param x a vector-valued object
#' @param indices a vector of indices defining the sub-blocks. Must match the length of the input vector.
#' @param ... additional arguments
#' @return A list of sub-blocks, where each sub-block contains the elements from the input object corresponding to the matching indices.
#' @export
#' @rdname split_blocks-methods
split_blocks <- function(x, indices, ...) standardGeneric("split_blocks")


#' Partition an image into a set of disjoint clusters
#'
#' This function partitions an image into a set of disjoint clusters using k-means clustering.
#'
#' @param x the image to partition, represented as a 3D array.
#' @param k the number of clusters to form.
#' @param ... additional arguments passed to the kmeans function.
#' @export
#' @rdname partition-methods
#' @examples
#' # Load an example 3D image
#' library(neuroim2)
#' img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Partition the image into 5 clusters using default options
#' clusters <- partition(img, 5)
#'
#'
#' @return a 3D array where each voxel is assigned to a cluster.
#' @seealso \code{\link{kmeans}}
setGeneric(name="partition", def=function(x, k, ...) standardGeneric("partition"))

#' Cut an object into a list of spatial or spatiotemporal clusters
#'
#' This function cuts an object into a list of sub-objects based on a vector of cluster indices.
#' The resulting list contains each of the clusters as separate objects.
#'
#' @param x The object to split. The input object must be a 4D tensor, where the first three dimensions correspond to the spatial dimensions of the data and the fourth dimension corresponds to time.
#' @param clusters A vector of cluster indices to split by.
#' @param ... Additional arguments to be passed to methods.
#' @return A list of sub-objects, where each sub-object corresponds to a unique cluster index.
#' @export
#' @rdname split_clusters-methods
setGeneric(name="split_clusters", def=function(x, clusters, ...) standardGeneric("split_clusters"))



#' Extract an ordered series of 2D slices from a 3D or 4D object
#'
#' This function extracts an ordered series of 2D slices from a 3D or 4D object. The returned slices are in the order they appear in the original object.
#'
#' @param x The 3D or 4D object to extract slices from
#' @param ... Additional arguments to be passed to the underlying methods
#'
#' @return A list of 2D matrices, each containing a slice from the input object
#'
#' @export
#' @rdname slices-methods
setGeneric(name="slices", def=function(x, ...) standardGeneric("slices"))




#' Extract the number of dimensions of an object
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
#' @return The number of dimensions of the input object `x`
#' @rdname ndim-methods
setGeneric(name="ndim", def=function(x, ...) standardGeneric("ndim"))


#' Get the length of a given dimension of an object
#'
#' This function returns the length of a given axis (dimension) of an object. The
#' axis can be specified using its position or name.
#'
#' @param x the object whose axis to query
#' @param axis an integer or character string indicating which axis to query
#' @export
#' @rdname dim_of-methods
#' @return An integer representing the length of the specified axis.
#' @examples
#'
#' x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
#' stopifnot(dim_of(x, x@axes@i) == 10)
setGeneric(name="dim_of", def=function(x, axis) standardGeneric("dim_of"))


#' Find Dimensions of a Given Axis
#'
#' This function returns the dimension of the specified axis for a given object, such as a matrix or an array.
#'
#' @param x An object representing the input data, such as a matrix or an array.
#' @param axis An integer representing the axis for which the dimension is requested. For example, 1 for rows, 2 for columns, or 3 for the depth in a 3D array.
#'
#' @return An integer representing the dimension of the specified axis for the given object.
#'
#' @export
#' @rdname which_dim-methods
#'
#' @examples
#'
#' x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
#' which_dim(x, x@axes@j) == 2
#'
setGeneric(name="which_dim", def=function(x, axis) standardGeneric("which_dim"))


#' Add a Dimension to an Object
#'
#' This function adds a new dimension to a given object, such as a matrix or an array.
#'
#' @param x A dimensioned object, such as a matrix, an array, or a NeuroSpace object.
#' @param n An integer representing the size of the dimension to add.
#'
#' @return An updated dimensioned object with the new dimension added.
#'
#' @examples
#' # Create a NeuroSpace object
#' x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#'
#' # Add a new dimension with size 10
#' x1 <- add_dim(x, 10)
#'
#' # Check the new dimension
#' ndim(x1) == 4
#' dim(x1)[4] == 10
#'
#' @export
#' @rdname add_dim-methods
setGeneric(name="add_dim", def=function(x, n) standardGeneric("add_dim"))

#' Drop a Dimension from an Object
#'
#' This function removes a specified dimension from a given object, such as a matrix or an array.
#'
#' @param x A dimensioned object, such as a NeuroSpace object.
#' @param dimnum An integer representing the index of the dimension to drop.
#'
#' @return An updated dimensioned object with the specified dimension removed.
#'
#' @examples
#' # Create a NeuroSpace object with dimensions (10, 10, 10)
#' x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#'
#' # Drop the first dimension
#' x1 <- drop_dim(x, 1)
#'
#' # Check the new dimensions
#' ndim(x1) == 2
#' dim(x1)[1] == 10
#'
#' @export
#' @rdname drop_dim-methods
setGeneric(name="drop_dim", def=function(x, dimnum) standardGeneric("drop_dim"))

#' Extract Geometric Properties of an Image
#'
#' This function retrieves the geometric properties of a given image, such as dimensions and voxel size.
#'
#' @param x The object to query, which can be an instance of \code{\linkS4class{NeuroVol}} or \code{\linkS4class{NeuroVec}}.
#' @param ... Additional arguments, if needed.
#'
#' @return An object representing the geometric space of the image, of type \code{\linkS4class{NeuroSpace}}.
#'
#' @examples
#' # Create a NeuroSpace object with dimensions (10, 10, 10) and voxel size (1, 1, 1)
#' x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#'
#' # Create a NeuroVol object with random data and the specified NeuroSpace
#' vol <- NeuroVol(rnorm(10 * 10 * 10), x)
#'
#' # Retrieve the geometric properties of the NeuroVol object
#' identical(x, space(vol))
#'
#' @export
#' @rdname space-methods
setGeneric(name="space", def=function(x, ...) standardGeneric("space"))


#' Fill Disjoint Sets of Values with the Output of a Function
#'
#' This function splits an object into disjoint sets of values based on a factor, applies a specified function to each set,
#' and returns a new object with the original values replaced by the function's output.
#'
#' @param x The object to split.
#' @param fac The \code{factor} to split by.
#' @param FUN The function used to summarize the sets.
#'
#' @return A new object in which the original values have been replaced by the output of the function.
#'
#' @export
#'
#' @details The \code{FUN} function can either return a scalar for each input vector or a vector equal to the length of the input vector.
#' If it returns a scalar, every voxel in the set will be filled with that value in the output vector.
#'
#' @examples
#' ## Summarize with mean -- FUN returns a scalar
#' x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#' vol <- NeuroVol(rnorm(10 * 10 * 10), x)
#' fac <- factor(rep(1:10, length.out=1000))
#' ovol.mean <- split_fill(vol, fac, mean)
#' identical(dim(ovol.mean), dim(vol))
#' length(unique(as.vector(ovol.mean))) == 10
#'
#' ## Transform by reversing vector -- FUN returns a vector
#' ovol2 <- split_fill(vol, fac, rev)
#'
#' @rdname split_fill-methods
setGeneric(name="split_fill", def=function(x, fac, FUN) standardGeneric("split_fill"))

#' Map Values from One Set to Another Using a User-supplied Lookup Table
#'
#' This function maps values from one set to another using a lookup table provided by the user.
#'
#' @param x The object from which values will be mapped.
#' @param lookup The lookup table. The first column is the "key" and the second column is the "value".
#'
#' @return A new object in which the original values have been replaced with the values from the lookup table.
#'
#' @export
#' @examples
#' x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#' vol <- NeuroVol(sample(1:10, 10 * 10 * 10, replace = TRUE), x)
#'
#' ## Lookup table is a list
#' lookup <- lapply(1:10, function(i) i * 10)
#' names(lookup) <- 1:10
#' ovol <- map_values(vol, lookup)
#'
#' ## Lookup table is a matrix. The first column is the key, and the second column is the value
#' names(lookup) <- 1:length(lookup)
#' lookup.mat <- cbind(as.numeric(names(lookup)), unlist(lookup))
#' ovol2 <- map_values(vol, lookup.mat)
#' all.equal(as.vector(ovol2), as.vector(ovol))
#'
#' @rdname map_values-methods
setGeneric(name="map_values", def=function(x, lookup) standardGeneric("map_values"))



#' Center and/or Scale Row-subsets of a Matrix or Matrix-like Object
#'
#' This function centers and/or scales the row-subsets of a numeric matrix or matrix-like object.
#'
#' @param x A numeric matrix or matrix-like object.
#' @param f The splitting object, typically a factor or a set of integer indices. Must be equal to the number of rows in the matrix.
#' @param center Should values within each submatrix be centered? If TRUE, the mean is removed from each column of the submatrix.
#' @param scale Should values be scaled? If TRUE, the vector is divided by the standard deviation for each column of the submatrix.
#'
#' @return A new matrix or matrix-like object in which the original rows have been grouped by 'f' and then centered and/or scaled for each grouping.
#'
#' @export
#' @examples
#'
#' M <- matrix(rnorm(1000), 10, 100)
#' fac <- factor(rep(1:2, each=5))
#' Ms <- split_scale(M, fac)
#'
#' ## Correctly centered
#' all(abs(apply(Ms[fac == 1,], 2, mean)) < .000001)
#' all(abs(apply(Ms[fac == 2,], 2, mean)) < .000001)
#'
#' ## Correctly scaled
#' all.equal(apply(Ms[fac == 1,], 2, sd), rep(1, ncol(Ms)))
#' all.equal(apply(Ms[fac == 2,], 2, sd), rep(1, ncol(Ms)))
#' @rdname split_scale-methods
setGeneric(name="split_scale", def=function(x, f, center, scale) standardGeneric("split_scale"))

#' Summarize Subsets of an Object by Splitting by Row and Applying a Summary Function
#'
#' This function summarizes subsets of a numeric matrix or matrix-like object by first splitting the object by row and then applying a summary function.
#'
#' @param x A numeric matrix or matrix-like object.
#' @param fac A factor to define subsets of the object.
#' @param FUN The summary function to apply to each subset. If not provided, the mean of each sub-matrix column is computed.
#'
#' @return A new matrix where the original values have been "reduced" by the supplied function.
#'
#' @details If 'FUN' is supplied, it must take a vector and return a single scalar value. If it returns more than one value, an error will occur.
#'
#' If 'x' is a NeuroVec instance, voxels (dimensions 1:3) are treated as columns and time-series (dimension 4) as rows.
#' The summary function is then applied to groups of voxels. However, if the goal is to apply a function to groups of time-points.
#'
#' @export
#' @examples
#' mat = matrix(rnorm(100*100), 100, 100)
#' fac = factor(sample(1:3, nrow(mat), replace=TRUE))
#' ## Compute column means of each sub-matrix
#' ms <- split_reduce(mat, fac)
#' all.equal(row.names(ms), levels(fac))
#'
#' ## Compute column medians of each sub-matrix
#' ms <- split_reduce(mat, fac, median)
#'
#' ## Compute time-series means grouped over voxels.
#' ## Here, 'length(fac)' must equal the number of voxels: 'prod(dim(bvec)[1:3])'
#' bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
#' fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
#' ms <- split_reduce(bvec, fac)
#' ms2 <- split_reduce(bvec, fac, mean)
#' all.equal(row.names(ms), levels(fac))
#' all.equal(ms, ms2)
#'
#' @rdname split_reduce-methods
setGeneric(name="split_reduce", def=function(x, fac, FUN) standardGeneric("split_reduce"))


#' Extract Voxel Dimensions of an Image
#'
#' This function extracts the voxel dimensions of an image represented by the input object.
#'
#' @param x The object representing the image.
#'
#' @return A numeric vector containing the voxel dimensions of the image.
#'
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
#' all.equal(spacing(bspace), c(2, 2, 2))
#'
#' @rdname spacing-methods
setGeneric(name="spacing", def=function(x) standardGeneric("spacing"))

#' Extract Spatial Bounds of an Image
#'
#' This function extracts the spatial bounds (origin + dim * spacing) of an image represented by the input object.
#'
#' @param x The object with the `bounds` property, typically an image.
#'
#' @return A matrix where each row contains the min (column 1) and max (column 2) bounds of the image dimension from 1 to `ndim(image)`.
#'
#' @examples
#' bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
#' b <- bounds(bspace)
#' nrow(b) == ndim(bspace)
#' ncol(b) == 2
#'
#' @rdname bounds-methods
setGeneric(name="bounds",     def=function(x) standardGeneric("bounds"))


#' Extract Image Axes
#'
#' @param x an object with a set of axes
#'
#' @return the `axes` associated with the object.
#'
#' @export
#' @examples
#' x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
#' class(axes(x)) == "AxisSet3D"
#'
#' @rdname axes-methods
setGeneric(name="axes",  def=function(x) standardGeneric("axes"))

#' Extract Image Origin
#'
#' @param x an object with an origin
#' @export
#'
#' @return the origin of the image
#'
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' stopifnot(origin(bspace) == c(0,0,0))
#'
#' @rdname origin-methods
setGeneric(name="origin", def=function(x) standardGeneric("origin"))


#' return the centroid of an object
#'
#' @param x an object with a centroid
#' @param ... extra args
#' @export
#' @return the centroid of the object
#' @examples
#'
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' centroid(bspace)
#'
#' @rdname centroid-methods
setGeneric(name="centroid", def=function(x, ...) standardGeneric("centroid"))

#' Return a matrix of centroids of an object
#'
#' @param x an object with multiple centroids (e.g. a \code{ClusteredNeuroVol})
#' @param ... extra args
#' @return A matrix where each row represents the coordinates of a centroid.
#' @export
#' @rdname centroids-methods
setGeneric(name="centroids", def=function(x, ...) standardGeneric("centroids"))



#' Extract image coordinate transformation
#'
#' @param x an object with a transformation
#' @return A 4x4 transformation matrix that maps from grid coordinates to real world coordinates.
#' @export
#' @details This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
#' see \code{\linkS4class{NeuroSpace}}
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' trans(bspace)
#' all.equal(dim(trans(bspace)), c(4,4))
#' @rdname trans-methods
setGeneric(name="trans",  def=function(x) standardGeneric("trans"))

#' Extract inverse image coordinate transformation
#'
#' @param x an object
#' @return A 4x4 transformation matrix that maps from real world coordinates back to grid coordinates.
#' @export
#' @examples
#' bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
#' itrans <- inverse_trans(bspace)
#' identical(trans(bspace) %*% inverse_trans(bspace), diag(4))
#' @rdname inverse_trans-methods
setGeneric(name="inverse_trans", def=function(x) standardGeneric("inverse_trans"))

#' Read a sequence of elements from an input source
#'
#' @param x the input channel
#' @param num_elements the number of elements to read
#' @return the elements as a vector
#' @export
#' @rdname read_elements-methods
setGeneric(name="read_elements", def=function(x, num_elements) standardGeneric("read_elements"))

#' Read a set of column vector from an input source (e.g. \code{ColumnReader})
#' @param x the input channel
#' @param column_indices the column indices
#' @return a \code{matrix} consisting of the requested column vectors
#' @export
#' @rdname read_columns-methods
setGeneric(name="read_columns", def=function(x, column_indices) standardGeneric("read_columns"))


#' Write a sequence of elements from an input source
#'
#' @param x the output channel
#' @param els the elements to write
#' @return Invisibly returns NULL after writing the elements.
#' @export
#' @rdname write_elements-methods
setGeneric(name="write_elements", def=function(x, els) standardGeneric("write_elements"))


#' Write a 3d image volume to disk
#'
#' @param x an image object, typically a \code{\linkS4class{NeuroVol}} instance.
#' @param file_name output file name
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param data_type output data type, If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
#' Otherwise output format will be inferred from R the datatype of the image.
#' @return Invisibly returns NULL after writing the volume to disk.
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
#' \donttest{
#' tmp1 <- tempfile(fileext = ".nii")
#' write_vol(bvol, tmp1)
#' unlink(tmp1)
#' }
#' @rdname write_vol-methods
setGeneric(name="write_vol",  def=function(x, file_name, format, data_type) standardGeneric("write_vol"))


#' Write a 4d image vector to disk
#'
#' @param x an image object, typically a \code{NeuroVec} instance.
#' @param file_name output file name.
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param data_type the numeric data type. If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
#' Otherwise output format will be inferred from R the datatype of the image.
#' @param ... extra args
#' @return Invisibly returns NULL after writing the vector to disk.
#' @export
#' @examples
#'
#' bvec <- NeuroVec(array(0, c(10,10,10,10)), NeuroSpace(c(10,10,10,10), c(1,1,1)))
#' \donttest{
#' # Create temporary files
#' tmp1 <- tempfile(fileext = ".nii")
#'
#' # Write vectors to temporary files
#' write_vec(bvec, tmp1)
#'
#' # Clean up
#' unlink(tmp1)
#' }
#' @rdname write_vec-methods
setGeneric(name="write_vec",  def=function(x, file_name, format, data_type, ...) standardGeneric("write_vec"))



#' Remap the grid-to-world coordinates mapping of an image.
#'
#' @param x the object
#' @param orient the orientation code indicating the "remapped" axes.
#' @return a reoriented image space
#' @rdname reorient-methods
#' @export
setGeneric(name="reorient", def=function(x, orient) standardGeneric("reorient"))


#' Convert 1d indices to n-dimensional grid coordinates
#'
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
#' @return A scaled version of the input 4D image, with each time series centered and/or scaled according to the parameters.
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
#' @return A sparse representation of the input object, containing only the elements specified by the mask.
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
#' @export
setGeneric(name="as.sparse", def=function(x, mask, ...) standardGeneric("as.sparse"))

#' Convert to dense representation
#'
#' @param x the object to densify
#' @return A dense representation of the input object.
#' @export
setGeneric(name="as.dense", def=function(x) standardGeneric("as.dense"))

#' Convert to a LogicalNeuroVol
#'
#' @param x the object to binarize
#' @param indices the indices to set to TRUE
#' @return A LogicalNeuroVol object with TRUE values at the specified indices.
#' @export
setGeneric(name="as.mask", def=function(x, indices) standardGeneric("as.mask"))



#' Generate a set of coordinate "patches" of fixed size from an image object.
#'
#' @param x the object to extract patches from
#' @param dims a vector indicating the dimensions of the patches
#' @param mask mask indicating the valid patch area
#' @param ... additional args
#' @return A list of coordinate patches, each representing a fixed-size region of the input object.
#' @export
setGeneric(name="patch_set", def=function(x, dims, mask, ...) standardGeneric("patch_set"))

#' Number of Clusters
#'
#' @param x the object to extract number of clusters
#' @return An integer representing the number of clusters in the object.
#' @export
setGeneric(name="num_clusters", def=function(x) standardGeneric("num_clusters"))


#' @title Extract coordinates from an object
#' @description This function extracts the coordinates from an input object.
#' @param x The object to extract coordinates from.
#' @param ... Additional arguments (not used in the generic function).
#' @return The extracted coordinates.
#' @export
#' @name coords
setGeneric(name="coords", def=function(x, ...) standardGeneric("coords"))

#' Extract indices
#'
#' @param x the object to extract indices
#' @return A vector of indices from the object.
#' @export
setGeneric(name="indices", def=function(x) standardGeneric("indices"))

#' Index Lookup operation
#' @param x the object to query
#' @param i the index to lookup
#' @param ... additional arguments
#' @return The value(s) at the specified index/indices.
#' @export
setGeneric(name="lookup", def=function(x, i, ...) standardGeneric("lookup"))


#' Extract time series from specific voxel coordinates and return as ROI object
#'
#' @description
#' Extracts time series data from a NeuroVec object at specified voxel coordinates
#' and returns it as an ROI object.
#'
#' @param x The NeuroVec object
#' @param i Numeric index for the first dimension
#' @param ... Additional arguments
#'
#' @return A ROIVec object containing the time series data for the specified coordinates
#'
#' @export
setGeneric("series_roi", function(x, i, ...) standardGeneric("series_roi"))

#' Extract one or more series from object
#' @param x the object
#' @param i the series indices
#' @param ... additional arguments
#' @return A list or array containing the extracted series.
#' @export
#' @rdname series-methods
setGeneric(name="series", def=function(x, i, ...) standardGeneric("series"))

#' Extract image slice
#'
#' Extract a 2D slice from an image volume
#'
#' @param x the object
#' @param zlevel coordinate (in voxel units) along the sliced axis
#' @param along the axis along which to slice
#' @param orientation the target orientation of the 2D slice
#' @param ... additional arguments
#' @return A 2D slice from the image volume.
#' @export
#' @rdname slice-methods
setGeneric(name="slice", def=function(x, zlevel, along, orientation, ...) standardGeneric("slice"))

#' Render an image to create a drawable image.
#'
#' Map an image intensities to an image with color values.
#'
#' @param x the object, e.g. an instance of type \code{NeuroSlice}
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @return A rendered image with the specified dimensions and color mapping.
#' @rdname render-methods
setGeneric(name="render", def=function(x, width, height, colmap,...) standardGeneric("render"))


#' Render a slice at z coordinate
#'
#' @param x the object, e.g. an instance of type \code{Layer} or \code{Overlay}
#' @param zpos the z coordinate to slice through.
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @return A rendered image of the specified slice with the given dimensions and color mapping.
#' @rdname render_slice-methods
setGeneric(name="render_slice", def=function(x, zpos, width, height, colmap,...) standardGeneric("render_slice"))




#' Extract permutation matrix associated with an image
#'
#' A permutation matrix defines how the native voxel coordinates can be transformed to standard (LPI) orientation.
#'
#' @param x the object
#' @param ... additional arguments
#' @export
#' @rdname perm_mat-methods
#' @return an N x N permutation matrix, where N is the dimensionality of the image.
#'
#' @details a permutation matrix can be used to convert between cardinal image orientations.
#' For example, if an image is stored in "RPI" (Right-Posterior-Inferior) format, a coordinate in this space
#' can be converted to LPI (Left-Posterior-Inferior) by multiplying a coordinate vector by the permutation matrix.
#'
#' @examples
#'
#' fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
#' vol <- read_vol(fname)
#' pmat <- perm_mat(space(vol))
#'
#' vox <- c(12,12,8)
#' pvox <- vox %*% perm_mat(space(vol))
#'
#' stopifnot(all(pvox == c(-12,12,8)))
setGeneric(name="perm_mat", def=function(x, ...) standardGeneric("perm_mat"))

#' Concatenate two objects in the time dimension
#'
#' @param x the first object, typically \code{NeuroVol} or \code{NeuroVec}
#' @param y the second object, typically \code{NeuroVol} or \code{NeuroVec}
#' @details The \code{x} and \code{y} images must have compatible dimensions. a \code{NeuroVol} can be concatenated to \code{NeuroVec}, and vice versa. See examples.
#' @param ... additional objects
#'
#' @return a temporally concatenated object.
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


#' Connected components
#'
#' Find connected components in an image
#'
#' @name conn_comp
#' @param x the image object
#' @param ... additional arguments
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


#' @name image
#' @title Generic Image Method for Creating Visual Representations
#'
#' @description Creates a visual representation (or image) from an object.
#'
#' @param x An object to be rendered as an image.
#' @param ... Additional arguments passed to methods.
#'
#' @export
if (!isGeneric("image"))
  setGeneric("image", function(x, ...) standardGeneric("image"))

#' @name as.raster
#' @title Generic Method for Converting Objects to Raster Format
#'
#' @description Converts an object to a raster (bitmap) representation.
#'
#' @param x An object to be converted.
#' @param ... Additional arguments passed to the conversion methods.
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

#' Generic function to read image meta info given a file
#' @param x descriptor instance
#' @param file_name file name containing meta information
#' @return A list containing the meta information read from the file.
#' @export
#' @rdname read_meta_info-methods
setGeneric(name="read_meta_info", def=function(x, file_name) standardGeneric("read_meta_info"))


#' Generic function to position kernel in a position in image space
#' @param x the kernel object
#' @param sp the space to embed the kernel
#' @param center_voxel the voxel marking the center of the kernel in the embedded space
#' @param ... extra args
#' @return An object representing the embedded kernel in the specified space.
#' @export
#' @rdname embed_kernel-methods
setGeneric("embed_kernel", def=function(x, sp, center_voxel, ...) standardGeneric("embed_kernel"))



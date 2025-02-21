#' @export
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' Generic Drop Method
#'
#' Provides a mechanism to remove dimensions or elements from an object.
#'
#' @param x An object.
#' @return An object of the same class as \code{x} with reduced dimensions or elements.
#'
#' @export
setGeneric("drop", function(x) standardGeneric("drop"))

#' Generic as.matrix Method
#'
#' Coerces an object to a matrix.
#'
#' @param x An object to be coerced to a matrix.
#' @param ... Additional arguments passed to methods.
#' @return A \code{matrix} representation of the input \code{x}.
#'
#' @export
setGeneric("as.matrix", function(x, ...) standardGeneric("as.matrix"))

#' Generic Scale Method
#'
#' Scales an object by (typically) subtracting the mean and dividing by the standard deviation.
#'
#' @param x The object to be scaled.
#' @param ... Additional arguments for scaling methods.
#' @return An object of the same class as \code{x}, scaled by the specified method.
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
#' @return An object representing the resampled \code{source} image, with the same spatial properties as \code{target}.
#'
#' @examples
#'
#' img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' rspace <- space(img)
#'
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
#' @return Invisibly returns \code{x}, called for its side effect of printing.
setGeneric(name="print_", def=function(x, ...) standardGeneric("print_"))

#' Extract Data Values of an Object
#'
#' @param x the object to get values from
#' @param ... additional arguments
#' @return A vector or array containing the values extracted from \code{x}.
#' @export
#' @rdname values-methods
#' @examples
#' x <- NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(rnorm(10 * 10 * 10), x)
#' values(vol)
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
#' @return A \code{vector} containing the values at the specified linear indices of \code{x}.
#' @export
#' @examples
#' # Create a sparse neuroimaging vector
#' bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
#' mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
#' mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
#' svec <- SparseNeuroVec(mat, bspace, mask)
#' 
#' # Extract values using linear indices
#' # Get values from first timepoint at voxels 1,2,3
#' indices <- c(1,2,3)
#' vals <- linear_access(svec, indices)
#' 
#' # Get values from multiple timepoints and voxels
#' # First voxel at timepoint 1, second voxel at timepoint 2
#' indices <- c(1, 1000 + 2) # 1000 = prod(10,10,10) 
#' vals <- linear_access(svec, indices)
setGeneric(name="linear_access", def=function(x, i, ...) standardGeneric("linear_access"))

#' Extract values from a 4D tensor using a matrix of time-space indices.
#'
#' This function efficiently extracts values from a 4D tensor (typically neuroimaging data) using a matrix of indices where each row contains
#' a time index in column 1 and a spatial index in column 2. The spatial index refers to the position
#' in the flattened spatial dimensions (x,y,z). This is primarily used internally by the \code{series()} 
#' method to efficiently access time series data for specific voxels.
#'
#' @param x a data source, typically a \code{SparseNeuroVec} object containing 4D neuroimaging data
#' @param i Either:
#'   \itemize{
#'     \item A matrix with 2 columns: [time_index, space_index] specifying which values to extract
#'     \item A numeric vector of spatial indices to extract all timepoints for those locations
#'   }
#' @param ... additional arguments to be passed to methods.
#' @return When \code{i} is a matrix, returns a numeric vector of values at the specified time-space coordinates.
#'         When \code{i} is a vector, returns a matrix where each column contains the full time series for each spatial index.
#' @rdname matricized_access-methods
#' @export
#' @examples 
#' # Create a sparse 4D neuroimaging vector
#' bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
#' mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
#' mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
#' svec <- SparseNeuroVec(mat, bspace, mask)
#' 
#' # Extract specific timepoint-voxel pairs
#' # Get value at timepoint 1, voxel 1 and timepoint 2, voxel 2
#' idx_mat <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
#' vals <- matricized_access(svec, idx_mat)
#' 
#' # Get full time series for voxels 1 and 2
#' ts_mat <- matricized_access(svec, c(1,2))
#' # Each column in ts_mat contains the full time series for that voxel
setGeneric(name="matricized_access", def=function(x, i, ...) standardGeneric("matricized_access"))

#' Read data from a data source.
#'
#' This function loads data from a data source and returns it in a format that is compatible with
#' other functions in the neuroim2 package. The format of the returned data depends on the type
#' of data source used.
#'
#' @param x a data source.
#' @param ... additional arguments to be passed to methods.
#' @return An R object containing loaded data, in a format compatible with the \pkg{neuroim2} package.
#' @keywords internal
#' @export
#' @examples
#' # Create a NeuroVolSource from a NIFTI file and load it
#' fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
#' src <- NeuroVolSource(fname)
#' vol <- load_data(src)
#' # The loaded volume is a DenseNeuroVol object
#' class(vol)
#' dim(vol)
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
#' @return The result of applying the mapping function to \code{x}.
#' @export
#' @examples
#' # Create a simple 3D volume
#' bspace <- NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), bspace)
#' 
#' # Create a 3x3x3 mean smoothing kernel
#' kern <- Kernel(c(3,3,3),  vdim=c(3,3,3))
#' 
#' # Apply the kernel to smooth the volume
#' smoothed_vol <- mapf(vol, kern)
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
#' @return A \code{list} containing the extracted 3D volumes from \code{x} in the same order as \code{indices}.
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
#' @return A \code{list} containing the extracted vectors from \code{x} in the same order as \code{subset}.
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
#' @return A \code{list} of sub-blocks, where each sub-block contains the elements from \code{x} corresponding to the matching \code{indices}.
#' @export
#' @examples
#' # Create a 4D neuroimaging vector with 20 timepoints
#' space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
#' vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)
#' 
#' # Split into 4 blocks by assigning timepoints to blocks 1-4 repeatedly
#' block_indices <- rep(1:4, length.out=20)
#' blocks <- split_blocks(vec, block_indices)
#'
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
#' @return a 3D \code{array} where each voxel is assigned to a cluster.
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
#' @return A \code{list} of sub-objects, where each sub-object corresponds to a unique cluster index.
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
#' @return A \code{list} of 2D \code{matrices}, each containing a slice from the input \code{x}.
#'
#' @export
#' @rdname slices-methods
#' @examples
#' # Create a simple 3D volume
#' space <- NeuroSpace(c(10,10,10), c(1,1,1))
#' vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), space)
#' 
#' # Get all slices along the z-axis
#' slc <- slices(vol)
#' 
#' # Number of slices equals the z dimension
#' length(slc) == dim(vol)[3]
#' 
#' # Each slice is a 2D matrix
#' dim(slc[[1]]) == c(10,10)
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
#' @return An integer representing the number of dimensions in \code{x}.
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
#' @return An integer representing the length of the specified \code{axis} of \code{x}.
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
#' @return An integer representing the dimension index of the specified \code{axis} for the object \code{x}.
#'
#' @export
#' @rdname which_dim-methods
#'
#' @examples
#'
#' x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
#' which_dim(x, x@axes@j) == 2
setGeneric(name="which_dim", def=function(x, axis) standardGeneric("which_dim"))

#' Add a Dimension to an Object
#'
#' This function adds a new dimension to a given object, such as a matrix or an array.
#'
#' @param x A dimensioned object, such as a matrix, an array, or a NeuroSpace object.
#' @param n An integer representing the size of the dimension to add.
#'
#' @return An object of the same class as \code{x} with the new dimension added.
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
#' @return An object of the same class as \code{x} with the specified dimension removed.
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
#' @return A \code{\linkS4class{NeuroSpace}} object representing the geometric space of \code{x}.
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
#' @return An object of the same class as \code{x}, with values replaced by the output of \code{FUN}.
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
#' @return An object of the same class as \code{x}, in which the original values have been replaced with the lookup table values.
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
#' @return An object of the same class as \code{x}, with row-subsets centered and/or scaled according to \code{f}.
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
#' @return A \code{matrix} (or matrix-like object) containing the summarized values after applying \code{FUN}.
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
#' @return A numeric vector specifying the voxel dimensions of \code{x}.
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
#' @return A numeric \code{matrix} with two columns specifying the min (column 1) and max (column 2) bounds of each dimension of \code{x}.
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
#' @return An object representing the \code{axes} of \code{x}.
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
#' @return A numeric vector giving the origin of \code{x}.
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
#' @return A numeric vector giving the centroid of \code{x}.
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
#' @return A numeric \code{matrix} where each row represents the coordinates of a centroid.
#' @export
#' @rdname centroids-methods
setGeneric(name="centroids", def=function(x, ...) standardGeneric("centroids"))

#' Extract image coordinate transformation
#'
#' @param x an object with a transformation
#' @return A numeric 4x4 matrix that maps from grid coordinates to real-world coordinates.
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
#' @return A numeric 4x4 \code{matrix} that maps from real-world coordinates back to grid coordinates.
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
#' @return A \code{vector} containing the elements read from \code{x}.
#' @export
#' @rdname read_elements-methods
#' @keywords internal
#' @examples
#' # Create a temporary binary file with test data
#' tmp <- tempfile()
#' con <- file(tmp, "wb")
#' test_data <- rnorm(100)
#' writeBin(test_data, con, size = 8)
#' close(con)
#' 
#' # Create a BinaryReader and read the data
#' reader <- BinaryReader(tmp, byte_offset = 0L,
#'                       data_type = "double", bytes_per_element = 8L)
#' data <- read_elements(reader, 100)
#' close(reader)
#' 
#' # Clean up
#' unlink(tmp)
setGeneric(name="read_elements", def=function(x, num_elements) standardGeneric("read_elements"))

#' Read a set of column vector from an input source (e.g. \code{ColumnReader})
#' @param x the input channel
#' @param column_indices the column indices
#' @return A numeric \code{matrix} consisting of the requested column vectors.
#' @export
#' @keywords internal
#' @rdname read_columns-methods
#' @examples
#' # Create a reader function that returns random data
#' reader_func <- function(cols) {
#'   matrix(rnorm(100 * length(cols)), 100, length(cols))
#' }
#' 
#' # Create a ColumnReader with 100 rows and 10 columns
#' col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)
#' 
#' # Read columns 1, 3, and 5
#' cols <- read_columns(col_reader, c(1L, 3L, 5L))
#' dim(cols) == c(100, 3)
#' 
setGeneric(name="read_columns", def=function(x, column_indices) standardGeneric("read_columns"))


#' Write a sequence of elements from an input source
#'
#' @param x the output channel
#' @param els the elements to write
#' @return Invisibly returns \code{NULL} after writing the elements.
#' @export
#' @rdname write_elements-methods
#' @examples
#' # Create a temporary binary file for writing
#' tmp <- tempfile()
#' writer <- BinaryWriter(tmp, byte_offset = 0L,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' 
#' # Write some random data
#' data <- rnorm(100)
#' write_elements(writer, data)
#' close(writer)
#' 
#' # Read back the data to verify
#' reader <- BinaryReader(tmp, byte_offset = 0L,
#'                       data_type = "double", bytes_per_element = 8L)
#' read_data <- read_elements(reader, 100)
#' close(reader)
#' 
#' # Verify data was written correctly
#' all.equal(data, read_data)
#' 
#' # Clean up
#' unlink(tmp)
setGeneric(name="write_elements", def=function(x, els) standardGeneric("write_elements"))


#' Write a 3d image volume to disk
#'
#' @param x an image object, typically a \code{\linkS4class{NeuroVol}} instance.
#' @param file_name output file name
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param data_type output data type, If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
#' Otherwise output format will be inferred from R the datatype of the image.
#' @return Invisibly returns \code{NULL} after writing the volume to disk.
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
#' @return Invisibly returns \code{NULL} after writing the vector to disk.
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
#' @return A reoriented version of \code{x}.
#' @details When \code{x} is a \code{NeuroSpace} object, the \code{orient} argument should be a character vector 
#' of length 3 specifying the desired anatomical orientation using single-letter codes. Each letter represents 
#' an anatomical direction:
#' \itemize{
#'   \item First position: "R" (Right) or "L" (Left)
#'   \item Second position: "A" (Anterior) or "P" (Posterior)
#'   \item Third position: "S" (Superior) or "I" (Inferior)
#' }
#' For example, \code{c("R", "A", "S")} specifies Right-Anterior-Superior orientation, while 
#' \code{c("L", "P", "I")} specifies Left-Posterior-Inferior orientation. The orientation codes 
#' determine how the voxel grid coordinates map to real-world anatomical space.
#' @rdname reorient-methods
#' @export
#' @examples
#' # Create a NeuroSpace object in LPI (Left-Posterior-Inferior) orientation
#' space <- NeuroSpace(c(64, 64, 40), c(2, 2, 2))
#' 
#' # Reorient to RAS (Right-Anterior-Superior) orientation
#' # Use individual axis codes: "R" for Right, "A" for Anterior, "S" for Superior
#' space_ras <- reorient(space, c("R", "A", "S"))
#' 
#' # The transformation matrix will be updated to reflect the new orientation
#' # Original and reoriented spaces will have different coordinate mappings
#' coords <- c(32, 32, 20)
#' orig_world <- grid_to_coord(space, coords)
#' new_world <- grid_to_coord(space_ras, coords)
setGeneric(name="reorient", def=function(x, orient) standardGeneric("reorient"))


#' Convert 1d indices to n-dimensional grid coordinates
#'
#' @param x the object
#' @param idx the 1d \code{vector} of indices
#' @return A numeric \code{matrix} of grid coordinates.
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
#' @return A numeric \code{matrix} of real-world coordinates.
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
#' @return An integer \code{vector} of 1D indices corresponding to \code{coords}.
#' @export
#' @examples
#' bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#' coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
#' idx <- coord_to_index(bvol, coords)
#' coords2 <- index_to_coord(bvol, idx)
#' all.equal(coords, coords2)
#' @rdname coord_to_index-methods
setGeneric(name="coord_to_index",   def=function(x, coords) standardGeneric("coord_to_index"))

#' convert n-dimensional real world coordinates to grid coordinates
#'
#' @param x the object
#' @param coords a matrix of real world coordinates
#' @return A numeric \code{matrix} of grid coordinates.
#' @export
#' @examples
#' # Create a simple 3D volume
#' bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#' coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
#' grid <- coord_to_grid(bvol, coords)
#' world <- grid_to_coord(bvol, grid)
#' all.equal(coords, world)
#' @rdname coord_to_grid-methods
setGeneric(name="coord_to_grid",   def=function(x, coords) standardGeneric("coord_to_grid"))

#' Generic function to convert N-dimensional grid coordinates to real world coordinates
#'
#' @param x the object
#' @param coords a matrix of grid coordinates
#' @return A numeric \code{matrix} of real-world coordinates.
#' @export
#' @examples
#' # Create a simple 3D volume
#' bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
#' grid_coords <- matrix(c(1.5,1.5,1.5, 5.5,5.5,5.5), ncol=3, byrow=TRUE)
#' world <- grid_to_coord(bvol, grid_coords)
#' grid <- coord_to_grid(bvol, world)
#' all.equal(grid_coords, grid)
#' @rdname grid_to_coord-methods
setGeneric(name="grid_to_coord",   def=function(x, coords) standardGeneric("grid_to_coord"))

#' Generic function to convert voxel coordinates in the reference space (LPI) to native array space.
#'
#' @param x the object
#' @param vox a matrix of LPI voxel coordinates
#' @return A numeric \code{matrix} of native voxel coordinates.
#' @export
#' @examples
#' # Create a simple 3D volume in LPI orientation
#' space <- NeuroSpace(c(10,10,10), c(2,2,2))
#' 
#' # Create a reoriented space in RAS orientation
#' space_ras <- reorient(space, c("R", "A", "S"))
#' 
#' # Convert coordinates between orientations
#' voxel_coords <- matrix(rbind(c(1,1,1)))
#' new_coords <- grid_to_grid(space_ras, voxel_coords)
#' 
#' @rdname grid_to_grid-methods
setGeneric(name="grid_to_grid",   def=function(x, vox) standardGeneric("grid_to_grid"))


#' Generic function to convert N-dimensional grid coordinates to 1D indices
#' @param x the object, typically a \code{NeuroVol} or \code{NeuroSpace} instance.
#' @param coords a matrix where each row is a coordinate or a vector of length equal to \code{ndim(x)}
#' @return An integer \code{vector} of 1D indices corresponding to \code{coords}.
#' @export
#' @examples
#' # Create a 2D space (10x10)
#' space_2d <- NeuroSpace(c(10,10), c(1,1))
#' 
#' # Convert 2D grid coordinates to linear indices
#' coords_2d <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
#' idx_2d <- grid_to_index(space_2d, coords_2d)
#' # First coordinate (1,1) maps to index 1
#' # Second coordinate (2,2) maps to index 12 (= 2 + (2-1)*10)
#' 
#' # Create a 3D space (10x10x10)
#' space_3d <- NeuroSpace(c(10,10,10), c(1,1,1))
#' 
#' # Convert 3D grid coordinates to linear indices
#' coords_3d <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)
#' idx_3d <- grid_to_index(space_3d, coords_3d)
#' 
#' # Single coordinate can also be converted
#' idx <- grid_to_index(space_3d, c(1,1,1))
#' 
#' @rdname grid_to_index-methods
setGeneric(name="grid_to_index",   def=function(x, coords) standardGeneric("grid_to_index"))



#' Generic function to extract a sub-vector from a \code{NeuroVec} object.
#' @param x four-dimensional image
#' @param i the indices of the volume(s) to extract
#' @param ... additional arguments
#' @return A \code{NeuroVec} object that is a sub-sequence of the supplied object.
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
#' @return An object of the same class as \code{x}, with each time series centered and/or scaled.
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
#' @return A sparse representation of the input object, containing only the elements specified by \code{mask}.
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
#' @return A \code{LogicalNeuroVol} object with \code{TRUE} values at the specified \code{indices}.
#' @export
setGeneric(name="as.mask", def=function(x, indices) standardGeneric("as.mask"))



#' Generate a set of coordinate "patches" of fixed size from an image object.
#'
#' @param x the object to extract patches from
#' @param dims a vector indicating the dimensions of the patches
#' @param mask mask indicating the valid patch area
#' @param ... additional args
#' @return A \code{list} of coordinate patches, each representing a fixed-size region of the input object.
#' @export
setGeneric(name="patch_set", def=function(x, dims, mask, ...) standardGeneric("patch_set"))

#' Number of Clusters
#'
#' @param x the object to extract number of clusters
#' @return An \code{integer} representing the number of clusters in \code{x}.
#' @export
setGeneric(name="num_clusters", def=function(x) standardGeneric("num_clusters"))


#' @title Extract coordinates from an object
#' @description This function extracts the coordinates from an input object.
#' @param x The object to extract coordinates from.
#' @param ... Additional arguments (not used in the generic function).
#' @return A numeric \code{matrix} or \code{vector} containing the coordinates of \code{x}.
#' @export
#' @name coords
setGeneric(name="coords", def=function(x, ...) standardGeneric("coords"))

#' Extract indices
#'
#' @param x the object to extract indices
#' @return A \code{vector} of indices from \code{x}.
#' @export
setGeneric(name="indices", def=function(x) standardGeneric("indices"))

#' Index Lookup operation
#' @param x the object to query
#' @param i the index to lookup
#' @param ... additional arguments
#' @return The value(s) at the specified index/indices of \code{x}.
#' @export
setGeneric(name="lookup", def=function(x, i, ...) standardGeneric("lookup"))


#' Extract time series from specific voxel coordinates and return as ROI object
#'
#' @description
#' Extracts time series data from a \code{NeuroVec} object at specified voxel coordinates
#' and returns it as an ROI object.
#'
#' @param x The \code{NeuroVec} object
#' @param i Numeric index for the first dimension
#' @param ... Additional arguments
#'
#' @return A \code{ROIVec} object containing the time series data for the specified coordinates.
#'
#' @export
setGeneric("series_roi", function(x, i, ...) standardGeneric("series_roi"))

#' Extract one or more series from object
#' @param x the object
#' @param i the series indices
#' @param ... additional arguments
#' @return A \code{list} or \code{array} containing the extracted series.
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
#' Map image intensities to an image with color values.
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
#' @return A numeric N x N \code{matrix} representing the permutation transform, where N is the dimensionality of the image.
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
#' @details The \code{x} and \code{y} images must have compatible dimensions. A \code{NeuroVol} can be concatenated to \code{NeuroVec}, and vice versa. See examples.
#' @param ... additional objects
#'
#' @return A temporally concatenated object.
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
#' @return An object representing the connected components of \code{x}.
setGeneric(name="conn_comp", def=function(x, ...) standardGeneric("conn_comp"))

#' extract voxel coordinates
#'
#' @param x the object to extract voxels from
#' @param ... additional arguments to function
#' @export
#' @rdname voxels-methods
#' @return A \code{matrix} or \code{vector} representing voxel coordinates from \code{x}.
setGeneric(name="voxels", def=function(x, ...) standardGeneric("voxels"))


#' @name image
#' @title Generic Image Method for Creating Visual Representations
#'
#' @description Creates a visual representation (or image) from an object.
#'
#' @param x An object to be rendered as an image.
#' @param ... Additional arguments passed to methods.
#'
#' @return An image object representing \code{x}.
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
#' @return A \code{raster} object representing \code{x}.
#' @rdname as.raster
if (!isGeneric("as.raster"))
  setGeneric("as.raster", function(x, ...) standardGeneric("as.raster"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{FileFormat}} instance.
#' Will test for match to either header file or data file
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return \code{TRUE} for match, \code{FALSE} otherwise.
#' @rdname file_matches-methods
#' @export
setGeneric(name="file_matches", def=function(x, file_name) standardGeneric("file_matches"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{FileFormat}} instance.
#' Will test for match to header file only
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return \code{TRUE} for match, \code{FALSE} otherwise.
#' @export
#' @rdname header_file_matches-methods
setGeneric(name="header_file_matches", def=function(x, file_name) standardGeneric("header_file_matches"))

#' Generic function to test whether a file name conforms to the given a \code{\linkS4class{FileFormat}} instance.
#' Will test for match to data file only
#' @param x object for which the file name is to matched to
#' @param file_name file name to be matched
#' @return \code{TRUE} for match, \code{FALSE} otherwise.
#' @export
#' @rdname data_file_matches-methods
setGeneric(name="data_file_matches", def=function(x, file_name) standardGeneric("data_file_matches"))

#' Generic function to get the name of the header file, given a file name and a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return The correct header file name as a \code{character} string.
#' @rdname header_file-methods
#' @export
setGeneric(name="header_file", def=function(x, file_name) standardGeneric("header_file"))

#' Generic function to get the name of the data file, given a file name and a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return The correct data file name as a \code{character} string.
#' @export
#' @rdname data_file-methods
setGeneric(name="data_file", def=function(x, file_name) standardGeneric("data_file"))

#' Generic function to strip extension from file name, given a \code{\linkS4class{FileFormat}} instance.
#' @param x descriptor instance
#' @param file_name file name to be stripped of its extension
#' @return A \code{character} string \code{file_name} without its extension.
#' @export
#' @rdname strip_extension-methods
setGeneric(name="strip_extension", def=function(x, file_name) standardGeneric("strip_extension"))

#' Generic function to read image meta info given a file
#' @param x descriptor instance
#' @param file_name file name containing meta information
#' @return A \code{list} containing the meta information read from the file.
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

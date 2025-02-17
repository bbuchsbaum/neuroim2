#' @include all_generic.R
#' @import methods
#' @import Matrix
NULL

setOldClass(c("file", "connection"))
setOldClass(c("gzfile", "connection"))
setOldClass("environment")
setOldClass("mmap")
#setOldClass("FBM")

#' ArrayLike5D Class
#'
#' A virtual class for representing five-dimensional array-like objects.
#' This class serves as an interface for objects that mimic 5D arrays.
#'
#' @export
setClass("ArrayLike5D")

#' ArrayLike4D Class
#'
#' A virtual class for representing four-dimensional array-like objects.
#' It is intended to serve as a base class for 4D array representations.
#'
#' @export
setClass("ArrayLike4D")

#' ArrayLike3D Class
#'
#' A virtual class for representing three-dimensional array-like objects.
#' It provides a common interface for 3D array operations.
#'
#' @export
setClass("ArrayLike3D")

#' NamedAxis
#'
#' This class represents an axis with a name attribute
#'
#' @rdname NamedAxis-class
#' @slot axis the name of the axis
#' @slot direction of axis (-1,+1)
#' @export
setClass("NamedAxis", representation(axis="character", direction="numeric"))

#'
#' AxisSet
#'
#' Virtual base class representing an ordered set of named axes.
#' @rdname AxisSet-class
#' @slot ndim the number of axes (or dimensions)
#' @export
setClass("AxisSet", representation(ndim="integer"))

#' AxisSet1D
#'
#' A one-dimensional axis set
#' @rdname AxisSet1D-class
#' @slot i the first axis
#' @export
setClass("AxisSet1D", representation(i="NamedAxis"), contains=c("AxisSet"))

#' AxisSet2D
#'
#' A two-dimensional axis set representing an ordered pair of named axes.
#'
#' @rdname AxisSet2D-class
#' @slot i The first axis, inherited from AxisSet1D
#' @slot j The second axis, of class "NamedAxis"
#' @seealso \code{\link{AxisSet1D-class}}, \code{\link{AxisSet3D-class}}
#' @examples
#' # Create an AxisSet2D object
#' axis1 <- new("NamedAxis", axis = "x", direction = 1)
#' axis2 <- new("NamedAxis", axis = "y", direction = 1)
#' axisSet2D <- new("AxisSet2D", i = axis1, j = axis2, ndim = 2L)
#' @export
setClass("AxisSet2D", representation(j="NamedAxis"), contains=c("AxisSet1D"))

#' AxisSet3D Class
#'
#' @description
#' A class representing a three-dimensional axis set, extending the AxisSet2D class
#' with an additional third axis.
#'
#' @slot k A \code{NamedAxis} object representing the third axis.
#'
#' @seealso \code{\link{AxisSet2D-class}}, \code{\link{NamedAxis-class}}
#'
#' @examples
#' # Create NamedAxis objects for each dimension
#' x_axis <- new("NamedAxis", axis = "x", direction = 1)
#' y_axis <- new("NamedAxis", axis = "y", direction = 1)
#' z_axis <- new("NamedAxis", axis = "z", direction = 1)
#'
#' # Create an AxisSet3D object
#' axis_set_3d <- new("AxisSet3D", i = x_axis, j = y_axis, k = z_axis, ndim = 3L)
#'
#' @export
#' @rdname AxisSet3D-class
setClass("AxisSet3D", representation(k="NamedAxis"), contains=c("AxisSet2D"))

#' AxisSet4D Class
#'
#' @description
#' A class representing a four-dimensional axis set, extending the AxisSet3D class
#' with an additional fourth axis.
#'
#' @slot l A \code{NamedAxis} object representing the fourth axis.
#'
#' @seealso \code{\link{AxisSet3D-class}}, \code{\link{NamedAxis-class}}
#'
#' @examples
#' # Create NamedAxis objects for each dimension
#' x_axis <- new("NamedAxis", axis = "x", direction = 1)
#' y_axis <- new("NamedAxis", axis = "y", direction = 1)
#' z_axis <- new("NamedAxis", axis = "z", direction = 1)
#' t_axis <- new("NamedAxis", axis = "t", direction = 1)
#'
#' # Create an AxisSet4D object
#' axis_set_4d <- new("AxisSet4D", i = x_axis, j = y_axis, k = z_axis,
#'                    l = t_axis, ndim = 4L)
#'
#' @export
#' @rdname AxisSet4D-class
setClass("AxisSet4D", representation(l="NamedAxis"), contains=c("AxisSet3D"))

#' AxisSet5D Class
#'
#' @description
#' A class representing a five-dimensional axis set, extending the AxisSet4D class
#' with an additional fifth axis.
#'
#' @slot m A \code{NamedAxis} object representing the fifth axis.
#'
#' @seealso \code{\link{AxisSet4D-class}}, \code{\link{NamedAxis-class}}
#'
#' @examples
#' # Create NamedAxis objects for each dimension
#' x_axis <- new("NamedAxis", axis = "x", direction = 1)
#' y_axis <- new("NamedAxis", axis = "y", direction = 1)
#' z_axis <- new("NamedAxis", axis = "z", direction = 1)
#' t_axis <- new("NamedAxis", axis = "t", direction = 1)
#' v_axis <- new("NamedAxis", axis = "v", direction = 1)
#'
#' # Create an AxisSet5D object
#' axis_set_5d <- new("AxisSet5D", i = x_axis, j = y_axis, k = z_axis,
#'                    l = t_axis, m = v_axis, ndim = 5L)
#'
#' @export
#' @rdname AxisSet5D-class
setClass("AxisSet5D", representation(m="NamedAxis"), contains=c("AxisSet4D"))

#' FileFormat Class
#'
#' @description
#' This class represents a neuroimaging file format descriptor, containing
#' information about the file format, encoding, and extensions for both header
#' and data components.
#'
#' @slot file_format A \code{character} string specifying the name of the file format (e.g., "NIfTI").
#' @slot header_encoding A \code{character} string specifying the file encoding of the header file
#'   (e.g., "raw" for binary, "gzip" for gz compressed).
#' @slot header_extension A \code{character} string specifying the file extension for the header file
#'   (e.g., "nii" for NIfTI single files).
#' @slot data_encoding A \code{character} string specifying the file encoding for the data file.
#' @slot data_extension A \code{character} string specifying the file extension for the data file
#'   (e.g., "nii" for NIfTI single files).
#'
#' @examples
#' # Create a FileFormat object for NIfTI format
#' nifti_format <- new("FileFormat",
#'                     file_format = "NIfTI",
#'                     header_encoding = "raw",
#'                     header_extension = "nii",
#'                     data_encoding = "raw",
#'                     data_extension = "nii")
#'
#'
#' @export
#' @rdname FileFormat-class
setClass("FileFormat",
         representation = representation(
           file_format = "character",
           header_encoding = "character",
           header_extension = "character",
           data_encoding = "character",
           data_extension = "character"
         )
)

#' NIFTIFormat Class
#'
#' @description
#' This class represents the NIFTI (Neuroimaging Informatics Technology Initiative) file format.
#' It extends the FileFormat class with NIFTI-specific attributes.
#'
#'
#' @keywords internal
#' @noRd
setClass("NIFTIFormat", contains = c("FileFormat"))

#' AFNIFormat Class
#'
#' @description
#' This class represents the AFNI (Analysis of Functional NeuroImages) file format.
#' It extends the FileFormat class with AFNI-specific attributes.
#'
#'
#' @keywords internal
#' @noRd
setClass("AFNIFormat", contains = c("FileFormat"))


#' MetaInfo Class
#'
#' @description
#' This class encapsulates meta information for neuroimaging data types,
#' including spatial and temporal characteristics, data type, and labeling.
#'
#' @slot data_type A \code{character} string specifying the data type code (e.g., "FLOAT", "INT").
#' @slot dims A \code{numeric} vector representing image dimensions.
#' @slot spatial_axes An \code{\linkS4class{AxisSet3D}} object representing image axes for spatial dimensions (x, y, z).
#' @slot additional_axes An \code{\linkS4class{AxisSet}} object representing axes for dimensions beyond spatial (e.g., time, color band, direction).
#' @slot spacing A \code{numeric} vector representing voxel dimensions in real-world units.
#' @slot origin A \code{numeric} vector representing the coordinate origin.
#' @slot label A \code{character} vector containing name(s) of images or data series.
#'
#' @details
#' The MetaInfo class provides a structured way to store and access essential
#' metadata for neuroimaging data. This includes information about the data type,
#' spatial and temporal dimensions, voxel spacing, and coordinate system origin.
#'
#' @seealso \code{\link{FileMetaInfo-class}}, \code{\link{AxisSet3D-class}}, \code{\link{AxisSet-class}}
#'
#' @examples
#' # Create a MetaInfo object
#' meta_info <- new("MetaInfo",
#'                  data_type = "FLOAT",
#'                  dims = c(64, 64, 32, 100),
#'                  spatial_axes = new("AxisSet3D"),
#'                  additional_axes = new("AxisSet"),
#'                  spacing = c(3, 3, 4),
#'                  origin = c(0, 0, 0),
#'                  label = "fMRI_run1")
#'
#' @export
#' @rdname MetaInfo-class
setClass("MetaInfo",
         representation = representation(
           data_type = "character",
           dims = "numeric",
           spatial_axes = "AxisSet3D",
           additional_axes = "AxisSet",
           spacing = "numeric",
           origin = "numeric",
           label = "character"
         ))

#' FileMetaInfo Class
#'
#' @description
#' This class extends MetaInfo to include file-specific metadata for neuroimaging data files.
#'
#' @slot header_file A \code{character} string specifying the name of the file containing meta information.
#' @slot data_file A \code{character} string specifying the name of the file containing image data.
#' @slot descriptor A \code{\linkS4class{FileFormat}} object describing the image file format.
#' @slot endian A \code{character} string specifying the byte order of data ('little' or 'big').
#' @slot data_offset A \code{numeric} value indicating the number of bytes preceding the start of image data in the data file.
#' @slot bytes_per_element An \code{integer} specifying the number of bytes per data element.
#' @slot intercept A \code{numeric} vector of constant values added to image data (one per sub-image).
#' @slot slope A \code{numeric} vector of multipliers for image data (one per sub-image).
#' @slot header A \code{list} of format-specific attributes.
#'
#' @seealso \code{\link{MetaInfo-class}}, \code{\link{NIFTIMetaInfo-class}}, \code{\link{AFNIMetaInfo-class}}
#'
#' @export
#' @rdname FileMetaInfo-class
setClass("FileMetaInfo",
         representation(header_file="character",
                        data_file="character",
                        descriptor="FileFormat",
                        endian="character",
                        data_offset="numeric",
                        bytes_per_element="integer",
                        intercept="numeric",
                        slope="numeric",
                        header="list"),
         contains="MetaInfo")

#' NIFTIMetaInfo Class
#'
#' @description
#' This class extends FileMetaInfo with NIfTI-specific metadata.
#'
#' @slot nifti_header A \code{list} of attributes specific to the NIfTI file format.
#'
#' @seealso \code{\link{FileMetaInfo-class}}
#'
#' @export
#' @rdname FileMetaInfo-class
setClass("NIFTIMetaInfo",
         representation(nifti_header="list"),
         contains=c("FileMetaInfo"))

#' AFNIMetaInfo Class
#'
#' @description
#' This class extends FileMetaInfo with AFNI-specific metadata.
#'
#' @slot afni_header A \code{list} of attributes specific to the AFNI file format.
#'
#' @seealso \code{\link{FileMetaInfo-class}}
#'
#' @export
#' @rdname FileMetaInfo-class
setClass("AFNIMetaInfo",
         representation(afni_header="list"),
         contains=c("FileMetaInfo"))

#' FileSource Class
#'
#' @description
#' Base class for representing a data source for images. The purpose of this class is to provide a layer in between
#' low level IO and image loading functionality.
#'
#' @slot meta_info An object of class \code{\linkS4class{FileMetaInfo}} containing meta information for the data source.
#'
#'
#' @export
#' @rdname FileSource-class
setClass("FileSource", representation(meta_info="FileMetaInfo"))

#' NeuroVolSource Class
#'
#' @description
#' A class used to produce a \code{\linkS4class{NeuroVol}} instance.
#'
#' @slot index An \code{integer} representing the index of the volume to be read (must be of length 1).
#'
#' @seealso \code{\link{FileSource-class}}, \code{\link{NeuroVol-class}}
#'
#' @keywords internal
#' @noRd
setClass("NeuroVolSource", representation(index="integer"), contains="FileSource")

#' NeuroVecSource Class
#'
#' @description
#' A class used to produce a \code{\linkS4class{NeuroVec}} instance.
#'
#' @slot indices An \code{integer} vector representing the indices of the volumes to be loaded.
#'
#' @seealso \code{\link{FileSource-class}}, \code{\link{NeuroVec-class}}
#'
#' @export
setClass("NeuroVecSource", representation(indices="integer"), contains="FileSource")


#' BinaryReader Class
#'
#' @description
#' Class supporting reading of bulk binary data from a connection
#'
#' @slot input The binary input connection
#' @slot byte_offset The number of bytes to skip at the start of input
#' @slot data_type The data type of the binary elements
#' @slot bytes_per_element The number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian The endianness of the binary input connection
#' @slot signed Logical indicating whether the data type is signed
#'
#' @export
#' @rdname BinaryReader-class
setClass("BinaryReader",
           representation(input="connection",
                          byte_offset="numeric",
                          data_type="character",
                          bytes_per_element="integer",
                          endian="character",
                          signed="logical"))

#' BinaryWriter Class
#'
#' @description
#' This class supports writing of bulk binary data to a connection
#'
#' @slot output The binary output connection
#' @slot byte_offset The number of bytes to skip at the start of input
#' @slot data_type The data type of the binary elements
#' @slot bytes_per_element The number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian The endianness of the binary output connection
#'
#' @export
#' @rdname BinaryWriter-class
setClass("BinaryWriter",
           representation(output="connection",
                          byte_offset="numeric",
                          data_type="character",
                          bytes_per_element="integer",
                          endian="character"))


#' NeuroSpace Class
#'
#' @description
#' The NeuroSpace class represents the geometric properties of a brain image,
#' including its dimensions, origin, spacing, axes, and coordinate transformations.
#' It provides a comprehensive framework for handling spatial information in
#' neuroimaging data analysis.
#'
#' @slot dim An integer vector representing the grid dimensions of the image.
#' @slot origin A numeric vector representing the coordinates of the spatial origin.
#' @slot spacing A numeric vector representing the dimensions (in mm) of the grid units (voxels).
#' @slot axes A named \code{\linkS4class{AxisSet}} object representing the set of spatial axes in the untransformed native grid space.
#' @slot trans A matrix representing an affine transformation that converts grid coordinates to real-world coordinates.
#' @slot inverse A matrix representing an inverse transformation that converts real-world coordinates to grid coordinates.
#'
#' @section Validity:
#' A \code{NeuroSpace} object is considered valid if:
#' \itemize{
#'   \item The length of the \code{dim} slot is equal to the lengths of the \code{spacing}, \code{origin}, and number of axes in the \code{axes} slots.
#'   \item The \code{dim} slot contains only non-negative values.
#' }
#'
#' @section Methods:
#' The following methods are available for \code{NeuroSpace} objects:
#' \itemize{
#'   \item \code{\link{dim}}: Get the dimensions of the space.
#'   \item \code{\link{origin}}: Get or set the origin of the space.
#'   \item \code{\link{spacing}}: Get or set the spacing of the space.
#'   \item \code{\link{axes}}: Get the axes of the space.
#'   \item \code{\link{trans}}: Apply the affine transformation to coordinates.
#' }
#'
#' @section Usage:
#' The \code{NeuroSpace} class is fundamental in representing and manipulating
#' the spatial properties of neuroimaging data. It is used extensively throughout
#' the package for operations that require spatial information, such as image
#' registration, resampling, and coordinate transformations.
#'
#' @examples
#' # Create a NeuroSpace object
#' space <- NeuroSpace(dim = c(64L, 64L, 64L),
#'                     origin = c(0, 0, 0),
#'                     spacing = c(1, 1, 1))
#'
#' # Get the dimensions
#' dim(space)
#'
#'
#'
#' @seealso
#' \code{\link{AxisSet-class}} for details on the axis set representation.
#' \code{\link{NeuroVol-class}} and \code{\link{NeuroVec-class}} for classes that use NeuroSpace.
#'
#' @references
#' For more information on spatial transformations in neuroimaging:
#' Brett, M., Johnsrude, I. S., & Owen, A. M. (2002). The problem of functional localization in the human brain. Nature Reviews Neuroscience, 3(3), 243-249.
#'
#' @export
#' @rdname NeuroSpace-class
setClass("NeuroSpace",
        representation(dim = "integer", origin = "numeric", spacing = "numeric",
                       axes = "AxisSet", trans = "matrix", inverse = "matrix"),

         validity = function(object) {
           dim <- object@dim
           if (any(dim < 0)) {
             return("'dim' slot must contain non-negative values")
           }
           TRUE
         })

#' NeuroObj Class
#'
#' @description
#' Base class for all neuroimaging data objects with a Cartesian spatial representation.
#' This class provides a foundation for more specific neuroimaging data structures.
#'
#' @slot space An object of class \code{\linkS4class{NeuroSpace}} representing the geometry of the image object.
#'
#' @seealso \code{\link{NeuroSpace-class}}, \code{\link{NeuroSlice-class}}, \code{\link{NeuroVol-class}}
#'
#' @export
#' @rdname NeuroObj-class
setClass("NeuroObj", representation(space="NeuroSpace"))

#' NeuroSlice Class
#'
#' @description
#' Represents a two-dimensional brain image slice. This class extends both the \code{array}
#' class for data storage and the \code{\linkS4class{NeuroObj}} class for spatial information.
#'
#' @details
#' NeuroSlice objects are typically used to represent individual slices of 3D brain volumes
#' or 2D projections of 3D data. They inherit the spatial properties from NeuroObj and
#' the data storage capabilities from array.
#'
#' @seealso \code{\link{NeuroObj-class}}, \code{\link{NeuroVol-class}}
#'
#' @examples
#' # Create a simple 2D brain slice
#' slice_data <- matrix(rnorm(64*64), 64, 64)
#' slice_space <- NeuroSpace(dim=c(64L, 64L), origin=c(0, 0), spacing=c(1, 1))
#' brain_slice <- new("NeuroSlice", .Data=slice_data, space=slice_space)
#'
#' @export
#' @rdname NeuroSlice-class
setClass("NeuroSlice", contains=c("array", "NeuroObj"))

#' NeuroVol Class
#'
#' @description
#' Base class for representing 3D volumetric neuroimaging data. This class extends
#' \code{\linkS4class{NeuroObj}} to provide a foundation for various types of 3D brain images.
#'
#' @details
#' NeuroVol serves as an abstract base class for more specific 3D neuroimaging data structures.
#' It inherits spatial properties from NeuroObj but does not specify a particular data storage method.
#'
#' @seealso \code{\link{NeuroObj-class}}, \code{\link{DenseNeuroVol-class}}
#'
#' @export
#' @rdname NeuroVol-class
setClass("NeuroVol", contains="NeuroObj")

#' DenseNeuroVol Class
#'
#' @description
#' Represents a three-dimensional brain image backed by a dense array. This class
#' combines the spatial properties of \code{\linkS4class{NeuroVol}} with the data
#' storage capabilities of an array.
#'
#' @details
#' DenseNeuroVol objects are used for 3D brain images where most or all voxels contain
#' meaningful data. They provide efficient access to individual voxel values and are
#' suitable for operations that require frequent random access to voxel data.
#'
#' @seealso \code{\link{NeuroVol-class}}, \code{\link{SparseNeuroVol-class}}
#'
#' @examples
#' # Create a simple 3D brain volume
#' vol_data <- array(rnorm(64*64*64), c(64, 64, 64))
#' vol_space <- NeuroSpace(dim=c(64L, 64L, 64L), origin=c(0, 0, 0), spacing=c(1, 1, 1))
#' brain_vol <- new("DenseNeuroVol", .Data=vol_data, space=vol_space)
#'
#' @export
#' @rdname DenseNeuroVol-class
setClass("DenseNeuroVol", contains=c("NeuroVol", "array"))

#' SparseNeuroVol Class
#'
#' @description
#' This class represents a three-dimensional brain image using a sparse data
#' representation. It is particularly useful for large brain images with a high
#' proportion of zero or missing values, offering efficient storage and processing.
#'
#' @details
#' The SparseNeuroVol class extends the \code{\linkS4class{NeuroVol}} class and
#' implements the ArrayLike3D interface. It uses a \code{sparseVector}
#' from the Matrix package to store the image data, which allows for memory-efficient
#' representation of sparse 3D neuroimaging data.
#'
#' @slot data A \code{sparseVector} object from the Matrix package,
#'   storing the image volume data in a sparse format.
#'
#' @seealso
#' \code{\link{NeuroVol-class}} for the base volumetric image class.
#' \code{\link{DenseNeuroVol-class}} for a dense representation of 3D brain images.
#'
#' @examples
#' \dontrun{
#' # Create a sparse 3D brain image
#' dim <- c(64L, 64L, 64L)
#' space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
#' sparse_data <- Matrix::sparseVector(x = c(1, 2, 3),
#'                                     i = c(100, 1000, 10000),
#'                                     length = prod(dim))
#' sparse_vol <- new("SparseNeuroVol", space = space, data = sparse_data)
#' }
#'
#' @references
#' Bates, D., & Maechler, M. (2019). Matrix: Sparse and Dense Matrix Classes
#' and Methods. R package version 1.2-18.
#' https://CRAN.R-project.org/package=Matrix
#'
#' @importFrom Matrix sparseVector
#'
#' @export
#' @rdname SparseNeuroVol-class
setClass("SparseNeuroVol",
         representation = representation(data = "sparseVector"),
         contains = c("NeuroVol", "ArrayLike3D"))

#' LogicalNeuroVol Class
#'
#' @description
#' This class represents a three-dimensional brain image where all values are
#' either TRUE or FALSE. It is particularly useful for creating and managing
#' binary masks for brain images.
#'
#' @details
#' The LogicalNeuroVol class extends the \code{\linkS4class{DenseNeuroVol}} class,
#' inheriting its spatial properties and array-based storage. However, it
#' constrains the values to be logical (TRUE or FALSE), making it ideal for
#' representing binary masks, regions of interest (ROIs), or segmentation results
#' in neuroimaging analyses.
#'
#' @slot .Data A logical array containing the binary volume data.
#' @slot space A \code{\linkS4class{NeuroSpace}} object defining the spatial properties of the volume.
#'
#' @section Methods:
#' This class inherits methods from \code{\linkS4class{DenseNeuroVol}}. Additional
#' methods specific to logical operations may be available.
#'
#' @seealso
#' \code{\link{DenseNeuroVol-class}} for the parent class.
#' \code{\link{NeuroVol-class}} for the base volumetric image class.
#'
#' @examples
#' # Create a simple logical brain volume (e.g., a mask)
#' dim <- c(64L, 64L, 64L)
#' mask_data <- array(sample(c(TRUE, FALSE), prod(dim), replace = TRUE), dim)
#' mask_space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
#' brain_mask <- new("LogicalNeuroVol", .Data = mask_data, space = mask_space)
#'
#' # Check the proportion of TRUE voxels
#' true_proportion <- sum(brain_mask) / prod(dim(brain_mask))
#' print(paste("Proportion of TRUE voxels:", true_proportion))
#'
#' @export
#' @rdname LogicalNeuroVol-class
setClass("LogicalNeuroVol", contains = c("DenseNeuroVol"))

#' ClusteredNeuroVol Class
#'
#' @description
#' This class represents a three-dimensional brain image divided into N disjoint
#' partitions or clusters. It extends the \code{\linkS4class{SparseNeuroVol}} class
#' to provide efficient storage and manipulation of clustered neuroimaging data.
#'
#' @slot mask A \code{\linkS4class{LogicalNeuroVol}} object representing the logical
#'   mask indicating the spatial domain of the set of clusters.
#' @slot clusters An integer vector representing the cluster number for each voxel
#'   in the mask.
#' @slot label_map A named list where each element represents a cluster and its name.
#' @slot cluster_map An \code{environment} object that maps from cluster id to the
#'   set of 1D spatial indices belonging to that cluster.
#'
#' @details
#' The ClusteredNeuroVol class is designed for efficient representation and
#' manipulation of brain images with distinct, non-overlapping regions or clusters.
#' It combines the memory efficiency of sparse representations with additional
#' structures for managing cluster information.
#'
#' @section Methods:
#' This class inherits methods from the \code{\linkS4class{SparseNeuroVol}} class.
#' Additional methods specific to cluster operations may be available.
#'
#' @section Usage:
#' ClusteredNeuroVol objects are particularly useful for:
#' \itemize{
#'   \item Representing parcellated brain images
#'   \item Storing results of clustering algorithms applied to neuroimaging data
#'   \item Efficient manipulation and analysis of region-based neuroimaging data
#' }
#'
#' @seealso
#' \code{\link{SparseNeuroVol-class}} for the parent sparse volume class.
#' \code{\link{LogicalNeuroVol-class}} for the mask representation.
#'
#' @examples
#' \dontrun{
#' # Create a simple clustered brain volume
#' dim <- c(10L, 10L, 10L)
#' mask_data <- array(rep(c(TRUE, FALSE), 500), dim)
#' mask <- new("LogicalNeuroVol", .Data = mask_data,
#'             space = NeuroSpace(dim = dim, origin = c(0,0,0), spacing = c(1,1,1)))
#'
#' clusters <- as.integer(runif(sum(mask_data)) * 5)
#' label_map <- list("Cluster1" = 1, "Cluster2" = 2, "Cluster3" = 3,
#'                   "Cluster4" = 4, "Cluster5" = 5)
#'
#' cluster_map <- new.env()
#' for (i in 1:5) {
#'   cluster_map[[as.character(i)]] <- which(clusters == i)
#' }
#'
#' clustered_vol <- new("ClusteredNeuroVol",
#'                      mask = mask,
#'                      clusters = clusters,
#'                      label_map = label_map,
#'                      cluster_map = cluster_map)
#' }
#'
#' @export
#' @rdname ClusteredNeuroVol-class
setClass("ClusteredNeuroVol",
         representation = representation(
           mask = "LogicalNeuroVol",
           clusters = "integer",
           label_map = "list",
           cluster_map = "environment"),
         contains = c("SparseNeuroVol"))

#' IndexLookupVol Class
#'
#' @description
#' A three-dimensional brain image class that serves as a map between 1D grid indices
#' and a table of values. This class is primarily used in conjunction with the
#' \code{\linkS4class{SparseNeuroVec}} class to efficiently represent and access
#' sparse neuroimaging data.
#'
#' @slot space A \code{\linkS4class{NeuroSpace}} object representing the 3D space of the brain image.
#' @slot indices An integer vector containing the 1D indices of the non-zero voxels in the grid.
#' @slot map An integer vector containing the mapping between the 1D indices and the table of values.
#'
#' @details
#' The IndexLookupVol class extends \code{\linkS4class{NeuroVol}} and provides a
#' mechanism for efficient lookup and mapping of sparse 3D neuroimaging data. It
#' stores only the indices of non-zero voxels and their corresponding mappings,
#' allowing for memory-efficient representation of large, sparse brain images.
#'
#' @section Methods:
#' This class inherits methods from \code{\linkS4class{NeuroVol}}. Additional
#' methods specific to index lookup and mapping operations may be available.
#'
#' @seealso
#' \code{\link{SparseNeuroVec-class}} for the primary class that utilizes IndexLookupVol.
#' \code{\link{NeuroVol-class}} for the base volumetric image class.
#'
#' @examples
#' # Create a NeuroSpace object
#' space <- NeuroSpace(dim = c(2L, 2L, 2L), origin = c(0, 0, 0), spacing = c(1, 1, 1))
#'
#' # Create a 3D mask
#' mask <- array(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE), dim = c(2, 2, 2))
#'
#' # Create indices and map for the IndexLookupVol
#' indices <- which(mask)
#' map <- seq_along(indices)
#'
#' # Create an IndexLookupVol object
#' ilv <- IndexLookupVol(space = space, indices = as.integer(indices))
#'
#' # Access the indices
#' print(ilv@indices)
#'
#' # Access the map
#' print(ilv@map)
#'
#' @export
#' @rdname IndexLookupVol-class
setClass("IndexLookupVol",
         representation = representation(
           space = "NeuroSpace",
           indices = "integer",
           map = "integer"
         ),
         contains = c("NeuroVol"))

#' NeuroVec Class
#'
#' @description
#' This S4 class represents a four-dimensional brain image, which is used to store
#' and process time series neuroimaging data such as fMRI or 4D functional
#' connectivity maps. The class extends the basic functionality of \code{\linkS4class{NeuroObj}}.
#'
#' @details
#' NeuroVec objects are designed to handle 4D neuroimaging data, where the first
#' three dimensions represent spatial coordinates, and the fourth dimension typically
#' represents time or another series dimension. This structure is particularly
#' useful for storing and analyzing functional MRI data, time series of brain
#' states, or multiple 3D volumes in a single object.
#'
#' @section Slots:
#' \describe{
#'   \item{space}{A \code{\linkS4class{NeuroSpace}} object defining the spatial properties of the image.}
#'   \item{label}{A character string providing a label for the NeuroVec object.}
#' }
#'
#' @section Methods:
#' Methods specific to NeuroVec objects may include operations for time series
#' analysis, 4D data manipulation, and extraction of 3D volumes or time courses.
#'
#' @section Usage:
#' To create a NeuroVec object, use the constructor function \code{NeuroVec()}.
#' This function should handle the appropriate initialization of the 4D data
#' structure and associated spatial information.
#'
#' @seealso
#' \code{\link{NeuroObj-class}} for the parent class.
#' \code{\link{DenseNeuroVec-class}} and \code{\link{SparseNeuroVec-class}} for specific implementations.
#'
#' @examples
#' \dontrun{
#' # Load an example 4D brain image
#' example_4d_image <- read_vec(system.file("extdata", "example_4d.nii", package = "neuroim2"))
#'
#' # Create a NeuroVec object
#' neuro_vec <- NeuroVec(data = array(rnorm(64*64*32*100), dim = c(64, 64, 32, 100)),
#'                       space = NeuroSpace(dim = c(64, 64, 32),
#'                       origin = c(0, 0, 0),
#'                       spacing = c(3, 3, 4))
#'
#' # Access the dimensions of the 4D image
#' dim(neuro_vec)
#'
#' # Extract a single 3D volume (e.g., the first time point)
#' first_volume <- neuro_vec[,,,1]
#' }
#'
#' @export
#' @rdname NeuroVec-class
setClass("NeuroVec",
         slots = c(label = "character"),
         prototype = list(label = ""),
         contains = c("NeuroObj"))

#' DenseNeuroVec Class
#'
#' @description
#' A class representing a four-dimensional brain image, backed by a dense array.
#' This class is designed for neuroimaging data where most voxels contain non-zero values.
#'
#' @details
#' DenseNeuroVec objects store their data in a dense array format, which is efficient
#' for operations that require frequent access to all voxels. This class inherits from
#' both \code{\linkS4class{NeuroVec}} and \code{array} classes, combining spatial
#' information with array-based storage.
#'
#' @section Validity:
#' A DenseNeuroVec object is considered valid if:
#' \itemize{
#'   \item The underlying data is a four-dimensional array.
#' }
#'
#' @seealso
#' \code{\link{NeuroVec-class}} for the parent class.
#' \code{\link{SparseNeuroVec-class}} for a sparse representation alternative.
#'
#' @examples
#' \dontrun{
#' # Create a simple 4D brain image
#' data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
#' space <- NeuroSpace(dim = c(64, 64, 32,10), origin = c(0, 0, 0), spacing = c(3, 3, 4))
#' dense_vec <- new("DenseNeuroVec", .Data = data, space = space)
#'
#' # Access dimensions
#' dim(dense_vec)
#'
#' # Extract a single 3D volume
#' first_volume <- dense_vec[[1]]
#' }
#'
#' @export
#' @rdname DenseNeuroVec-class
setClass("DenseNeuroVec", contains = c("array", "NeuroVec"))

#' @keywords internal
#' @noRd
setValidity("DenseNeuroVec", function(object) {
  if (length(dim(object)) != 4) {
    "data must be four-dimensional array"
  } else {
    TRUE
  }
})



#' MappedNeuroVec Class
#'
#' @description
#' A class representing a four-dimensional brain image backed by a memory-mapped file.
#' This class provides efficient access to large brain images without loading the
#' entire dataset into memory.
#'
#' @slot filemap An object of class \code{mmap} representing the memory-mapped file
#'   containing the brain image data.
#' @slot offset An integer representing the byte offset within the memory-mapped file
#'   where the brain image data starts.
#'
#' @details
#' MappedNeuroVec objects use memory-mapped files to store and access large 4D brain
#' images efficiently. This approach allows for rapid access to specific portions of
#' the data without requiring the entire dataset to be loaded into memory at once.
#'
#' @section Methods:
#' This class inherits methods from \code{\linkS4class{NeuroVec}} and implements
#' the \code{ArrayLike4D} interface. Additional methods specific to memory-mapped
#' operations may be available.
#'
#' @seealso
#' \code{\link{NeuroVec-class}} for the parent class.
#' \code{\link[mmap]{mmap}} for details on memory-mapped file objects.
#'
#' @examples
#' \dontrun{
#' # Create a MappedNeuroVec object (pseudo-code)
#' file_path <- "/path/to/large/brain/image.dat"
#' mapped_vec <- MappedNeuroVec(file_path, dim = c(64, 64, 32, 100))
#'
#' # Access a subset of the data
#' subset <- mapped_vec[,,, 1:10]
#' }
#'
#' @importFrom mmap mmap
#' @export
#' @rdname MappedNeuroVec-class
setClass("MappedNeuroVec",
         slots = c(
           filemap = "mmap",
           offset = "integer"
         ),
         prototype = list(label = ""),
         contains = c("NeuroVec", "ArrayLike4D"))

#' AbstractSparseNeuroVec Class
#'
#' @description
#' An abstract base class for sparse four-dimensional brain image representations.
#' This class provides the foundation for efficient storage and manipulation of
#' large, sparse neuroimaging data.
#'
#' @slot mask An object of class \code{\linkS4class{LogicalNeuroVol}} defining the
#'   sparse domain of the brain image. This mask indicates which voxels contain
#'   non-zero data.
#' @slot map An object of class \code{\linkS4class{IndexLookupVol}} used to map
#'   between spatial coordinates and index/row coordinates in the sparse representation.
#'
#' @details
#' The AbstractSparseNeuroVec class serves as a template for implementing various
#' sparse representations of 4D brain images. It combines the spatial properties of
#' \code{\linkS4class{NeuroVec}} with the efficiency of sparse data structures.
#'
#' @section Subclasses:
#' Concrete implementations of this abstract class should provide specific data
#' storage mechanisms and methods for efficient access and manipulation of sparse
#' 4D brain image data.
#'
#' @seealso
#' \code{\link{NeuroVec-class}} for the parent class.
#' \code{\link{LogicalNeuroVol-class}} for the mask representation.
#' \code{\link{IndexLookupVol-class}} for the spatial-to-index mapping.
#'
#' @examples
#' \dontrun{
#' # This is an abstract class and should not be instantiated directly.
#' # Instead, use one of its concrete subclasses, such as SparseNeuroVec.
#' }
#'
#' @export
#' @rdname AbstractSparseNeuroVec-class
setClass("AbstractSparseNeuroVec",
         slots = c(
           mask = "LogicalNeuroVol",
           map = "IndexLookupVol"
         ),
         contains = c("NeuroVec", "ArrayLike4D"))

#' SparseNeuroVec Class
#'
#' @description
#' A class representing a sparse four-dimensional brain image, optimized for
#' efficient storage and access of large, sparse neuroimaging data.
#'
#' @slot data A \code{matrix} where each column represents a non-zero vector
#'   spanning the fourth dimension (e.g., time series for each voxel). Rows
#'   correspond to voxels in the sparse domain defined by the mask.
#'
#' @details
#' SparseNeuroVec objects store data in a compressed format, where only non-zero
#' values are retained. This approach significantly reduces memory usage for
#' sparse brain images. The class leverages the mask and mapping from its parent
#' class \code{\linkS4class{AbstractSparseNeuroVec}} to efficiently manage the
#' spatial structure of the data.
#'
#' @section Inheritance:
#' \code{SparseNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\linkS4class{NeuroVec}}: Base class for 4D brain images
#'   \item \code{\linkS4class{AbstractSparseNeuroVec}}: Provides sparse representation framework
#'   \item \code{ArrayLike4D}: Interface for 4D array-like operations
#' }
#'
#' @seealso
#' \code{\link{AbstractSparseNeuroVec-class}} for the parent sparse representation class.
#' \code{\link{NeuroVec-class}} for the base 4D brain image class.
#'
#' @examples
#' \dontrun{
#' # Create a sparse 4D brain image
#' mask <- LogicalNeuroVol(array(runif(64*64*32) > 0.7, dim=c(64,64,32)))
#' data <- matrix(rnorm(sum(mask) * 100), nrow=sum(mask), ncol=100)
#' sparse_vec <- SparseNeuroVec(data=data, mask=mask, space=NeuroSpace(dim=c(64,64,32)))
#'
#' # Access a subset of the data
#' subset <- sparse_vec[,,, 1:10]
#' }
#'
#' @export
#' @rdname SparseNeuroVec-class
setClass("SparseNeuroVec",
         slots = c(
           data = "matrix"
         ),
         contains = c("NeuroVec", "AbstractSparseNeuroVec", "ArrayLike4D"))



#' NeuroHyperVec Class
#'
#' @description
#' A class representing a five-dimensional brain image, where the first three dimensions are spatial,
#' the fourth dimension is typically time or trials, and the fifth dimension represents features within a trial.
#'
#' @slot mask An object of class \code{\linkS4class{LogicalNeuroVol}} defining the sparse spatial domain of the brain image.
#' @slot data A 3D array with dimensions [features x trials x voxels] containing the neuroimaging data.
#' @slot space A \code{\linkS4class{NeuroSpace}} object representing the dimensions and voxel spacing of the neuroimaging data.
#' @slot lookup_map An integer vector for O(1) spatial index lookups.
#'
#' @export
#' @rdname NeuroHyperVec-class
setClass(
  "NeuroHyperVec",
  slots = c(
    mask = "LogicalNeuroVol",
    data = "array",  # Dimensions: [features x trials x voxels]
    space = "NeuroSpace",
    lookup_map = "integer"
  ),
  contains = c( "ArrayLike5D"),
  validity = function(object) {
    # Validate that data is a 3D array
    if (!is.array(object@data) || length(dim(object@data)) != 3) {
      return("Data must be a 3D array with dimensions [features x trials x voxels]")
    }

    # Get expected dimensions
    num_voxels <- sum(object@mask@.Data)
    num_trials <- dim(object@space)[4]
    num_features <- dim(object@space)[5]

    # Validate array dimensions
    expected_dims <- c(num_features, num_trials, num_voxels)
    if (!identical(dim(object@data), expected_dims)) {
      return(sprintf(
        "Data array dimensions [%s] do not match expected [%d x %d x %d] (features x trials x voxels)",
        paste(dim(object@data), collapse=" x "),
        num_features, num_trials, num_voxels))
    }

    # Validate that mask is consistent with spatial dimensions
    mask_dims <- dim(object@mask)
    if (length(mask_dims) != 3) {
      return("Mask must be a 3D volume")
    }

    # Validate lookup_map
    if (length(object@lookup_map) != prod(dim(object@mask))) {
      return("lookup_map length must match the total number of voxels in the mask")
    }

    TRUE
  }
)



#' BigNeuroVec Class
#'
#' @description
#' A class representing a sparse four-dimensional brain image backed by a disk-based
#' big matrix. BigNeuroVec objects are designed for efficient handling of large-scale
#' brain imaging data that exceeds available memory.
#'
#' @slot data An instance of class \code{FBM} from the \code{bigstatsr} package,
#'   containing time-series data. The FBM (File-Backed Big Matrix) is a matrix-like
#'   structure stored on disk, enabling efficient handling of large-scale data.
#'
#' @details
#' BigNeuroVec leverages file-backed storage to manage large 4D neuroimaging datasets
#' that would typically exceed available RAM. It combines the sparse representation
#' framework of \code{\linkS4class{AbstractSparseNeuroVec}} with the disk-based
#' storage capabilities of \code{FBM}, allowing for out-of-core computations on
#' massive datasets.
#'
#' @section Inheritance:
#' \code{BigNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\linkS4class{NeuroVec}}: Base class for 4D brain images
#'   \item \code{\linkS4class{AbstractSparseNeuroVec}}: Provides sparse representation framework
#'   \item \code{ArrayLike4D}: Interface for 4D array-like operations
#' }
#'
#' @seealso
#' \code{\link{AbstractSparseNeuroVec-class}} for the parent sparse representation class.
#' \code{\link{NeuroVec-class}} for the base 4D brain image class.
#' \code{\link[bigstatsr]{FBM}} for details on File-Backed Big Matrix objects.
#'
#' @examples
#' \dontrun{
#' # Create a BigNeuroVec object
#' library(bigstatsr)
#'
#' # Create a file-backed big matrix
#' fbm <- FBM(10000, 1000, init = rnorm(10000000))
#'
#' # Create a mask for sparse representation
#' mask <- LogicalNeuroVol(array(runif(100*100*100) > 0.7, dim=c(100,100,100)))
#'
#' # Create a BigNeuroVec object
#' big_vec <- BigNeuroVec(data = fbm, mask = mask, space = NeuroSpace(dim=c(100,100,100)))
#'
#' # Access a subset of the data
#' subset <- big_vec[,,, 1:10]
#' }
#'
#' @importFrom bigstatsr FBM
#' @export
#' @rdname BigNeuroVec-class
setClass("BigNeuroVec",
         slots = c(
           data = "FBM"
         ),
         contains = c("NeuroVec", "AbstractSparseNeuroVec", "ArrayLike4D"))

#' FileBackedNeuroVec Class
#'
#' @description
#' A class representing a four-dimensional brain image that uses on-demand loading
#' through memory-mapped file access. This approach enables efficient handling of
#' large-scale brain imaging data by loading only the required portions of the data
#' into memory when needed.
#'
#' @slot meta An instance of class \code{\linkS4class{FileMetaInfo}} containing
#'   file metadata such as file path, format, and other associated information.
#'
#' @details
#' FileBackedNeuroVec objects provide a memory-efficient solution for working with
#' large 4D neuroimaging datasets. By utilizing memory-mapped file access, this class
#' allows users to work with datasets that exceed available RAM, only loading the
#' necessary data segments into memory as they are accessed.
#'
#' @section Inheritance:
#' \code{FileBackedNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\linkS4class{NeuroVec}}: Base class for 4D brain images
#'   \item \code{ArrayLike4D}: Interface for 4D array-like operations
#' }
#'
#' @seealso
#' \code{\link{NeuroVec-class}} for the base 4D brain image class.
#' \code{\link{FileMetaInfo-class}} for details on file metadata representation.
#'
#' @examples
#' \dontrun{
#' # Create a FileBackedNeuroVec object
#' file_path <- "/path/to/large/brain/image.nii.gz"
#' file_meta <- FileMetaInfo(file_path, format = "NIFTI")
#' file_backed_vec <- FileBackedNeuroVec(meta = file_meta)
#'
#' # Access a subset of the data (this will load only the required portion)
#' subset <- file_backed_vec[,,, 1:10]
#' }
#'
#' @export
#' @rdname FileBackedNeuroVec-class
setClass("FileBackedNeuroVec",
         slots = c(
           meta = "FileMetaInfo"
         ),
         contains = c("NeuroVec", "ArrayLike4D"))

#' NeuroVecSeq Class
#'
#' @description
#' A concatenated sequence of \code{\linkS4class{NeuroVec}} instances.
#'
#' @slot vecs The sequences of \code{NeuroVec} instances
#' @slot lens The number of volumes in each \code{NeuroVec} sequence
#'
#' @export
#' @rdname NeuroVecSeq-class
setClass("NeuroVecSeq",
         representation(vecs="list", lens="numeric"),
         contains=c("NeuroVec", "ArrayLike4D"),

         validity = function(object) {
           assert_that(all(purrr::map_lgl(object@vecs, ~ inherits(., "NeuroVec"))))
           dimlist <- purrr::map(object@vecs, ~ dim(.)[1:3])
           splist <- purrr::map(object@vecs, ~ spacing(.))
           assert_that(all(purrr::map_lgl(dimlist, ~ all(dimlist[[1]] == .))))
           assert_that(all(purrr::map_lgl(splist, ~ all(splist[[1]] == .))))

         })





#' SparseNeuroVecSource Class
#'
#' @description
#' A class used to produce a \code{\linkS4class{SparseNeuroVec}} instance. It encapsulates
#' the necessary information to create a sparse representation of a 4D neuroimaging dataset.
#'
#' @slot mask An object of class \code{\linkS4class{LogicalNeuroVol}} representing the subset
#'   of voxels that will be stored in memory. This mask defines the sparse structure of the
#'   resulting SparseNeuroVec.
#'
#' @details
#' SparseNeuroVecSource acts as a factory for SparseNeuroVec objects. It holds the
#' spatial mask that determines which voxels will be included in the sparse representation.
#' This class is typically used in data loading or preprocessing pipelines where the
#' sparse structure of the data is known or determined before the full dataset is loaded.
#'
#' @section Inheritance:
#' \code{SparseNeuroVecSource} inherits from:
#' \itemize{
#'   \item \code{\linkS4class{NeuroVecSource}}: Base class for NeuroVec source objects
#' }
#'
#' @seealso
#' \code{\link{SparseNeuroVec-class}} for the resulting sparse 4D neuroimaging data class.
#' \code{\link{LogicalNeuroVol-class}} for the mask representation.
#'
#' @examples
#' \dontrun{
#' # Create a simple mask
#' mask_data <- array(runif(64*64*32) > 0.7, dim = c(64, 64, 32))
#' mask <- LogicalNeuroVol(mask_data, space = NeuroSpace(dim = c(64, 64, 32)))
#'
#' # Create a SparseNeuroVecSource
#' sparse_source <- new("SparseNeuroVecSource", mask = mask)
#'
#' # Use the source to create a SparseNeuroVec (pseudo-code)
#' # sparse_vec <- create_sparse_neuro_vec(sparse_source, data)
#' }
#'
#' @export
#' @rdname SparseNeuroVecSource-class
setClass("SparseNeuroVecSource",
         slots = c(mask = "LogicalNeuroVol"),
         contains = "NeuroVecSource")

#' MappedNeuroVecSource Class
#'
#' @description
#' A class used to produce a \code{\linkS4class{MappedNeuroVec}} instance. It encapsulates
#' the necessary information to create a memory-mapped representation of a 4D neuroimaging dataset.
#'
#' @details
#' MappedNeuroVecSource acts as a factory for MappedNeuroVec objects. While it doesn't
#' have any additional slots beyond its parent class, it specifies the intent to create
#' a memory-mapped representation of the neuroimaging data. This class is typically used
#' in data loading pipelines where large datasets need to be accessed efficiently without
#' loading the entire dataset into memory.
#'
#' @section Inheritance:
#' \code{MappedNeuroVecSource} inherits from:
#' \itemize{
#'   \item \code{\linkS4class{NeuroVecSource}}: Base class for NeuroVec source objects
#' }
#'
#' @examples
#' \dontrun{
#' # Create a MappedNeuroVecSource
#' mapped_source <- new("MappedNeuroVecSource")
#'
#' # Use the source to create a MappedNeuroVec (pseudo-code)
#' # file_path <- "/path/to/large/brain/image.nii"
#' # mapped_vec <- create_mapped_neuro_vec(mapped_source, file_path)
#' }
#'
#' @export
#' @rdname MappedNeuroVecSource-class
setClass("MappedNeuroVecSource", contains = "NeuroVecSource")

#' numericOrMatrix Union
#'
#' A class union that includes both numeric vectors and matrices.
#'
#' @export
setClassUnion("numericOrMatrix", c("numeric", "matrix"))

#' ROI
#'
#' Base marker class for a region of interest (ROI)
#'
#' @export
setClass("ROI", contains="VIRTUAL")


#' ROICoords
#'
#' A class representing a region of interest (ROI) in a brain image, defined by a set of coordinates. This class stores the geometric space of the image and the coordinates of the voxels within the ROI.
#'
#' @slot space An instance of class \code{\linkS4class{NeuroSpace}} representing the geometric space of the image data.
#' @slot coords A \code{matrix} containing the coordinates of the voxels within the ROI.
#' Each row represents a coordinate as, e.g. (i,   j,  k).
#'
#'
#' @name ROICoords-class
#' @export
setClass("ROICoords",
         representation=representation(space="NeuroSpace", coords="matrix"),
         contains=c("ROI"))

#' ROIVol
#'
#' A class representing a volumetric region of interest (ROI) in a brain image, defined by a set of coordinates and associated data values.
#'
#' @slot coords A \code{matrix} containing the 3D coordinates of the voxels within the ROI. Each row represents a voxel coordinate as (x, y, z).
#' @slot .Data A \code{numeric} vector containing the data values associated with each voxel in the ROI. The length of this vector should match the number of rows in the \code{coords} matrix.
#'
#'
#' @section Validity:
#' An object of class \code{ROIVol} is considered valid if:
#' - The \code{coords} slot is a matrix with 3 columns.
#' - The \code{.Data} slot is a numeric vector.
#' - The length of the \code{.Data} vector is equal to the number of rows in the \code{coords} matrix.
#'
#' @name ROIVol-class
#' @export
setClass("ROIVol",
         contains=c("ROICoords", "numeric"),
         validity = function(object) {
           if (ncol(object@coords) != 3) {
             stop("coords slot must be a matrix with 3 columns")
           }
           if (!is.vector(object@.Data)) {
             stop("'data' must be a vector")
           }
           if (length(object@.Data) != nrow(object@coords)) {

             stop("length of data vector must equal 'nrow(coords)'")
           }
         })


#' ROIVolWindow
#'
#' A class representing a spatially windowed volumetric region of interest (ROI) in a brain image, derived from a larger parent ROI.
#'
#' @slot parent_index An \code{integer} representing the 1D index of the center voxel in the parent space.
#' @slot center_index An \code{integer} representing the location in the coordinate matrix of the center voxel in the window.
#' @slot coords A \code{matrix} containing the 3D coordinates of the voxels within the ROI. Each row represents a voxel coordinate as (x, y, z).
#' @slot .Data A \code{numeric} vector containing the data values associated with each voxel in the ROI. The length of this vector should match the number of rows in the \code{coords} matrix.
#'
#'
#' @section Validity:
#' An object of class \code{ROIVolWindow} is considered valid if:
#' - The \code{coords} slot is a matrix with 3 columns.
#' - The \code{.Data} slot is a numeric vector.
#' - The length of the \code{.Data} vector is equal to the number of rows in the \code{coords} matrix.
#'
#' @name ROIVolWindow-class
#' @export
setClass("ROIVolWindow",
         representation=representation(parent_index="integer", center_index="integer"),
         contains=c("ROIVol"),
         validity = function(object) {
           if (ncol(object@coords) != 3) {
             stop("coords slot must be a matrix with 3 columns")
           }
           if (!is.vector(object@.Data)) {
             stop("'data' must be a vector")
           }
           if (length(object@.Data) != nrow(object@coords)) {
             stop("length of data vector must equal 'nrow(coords)'")
           }
         })




#' ROIVec
#'
#' A class representing a vector-valued volumetric region of interest (ROI) in a brain image.
#'
#' @slot coords A \code{matrix} containing the 3D coordinates of the voxels within the ROI. Each row represents a voxel coordinate as (x, y, z).
#' @slot .Data A \code{matrix} containing the data values associated with each voxel in the ROI. Each row corresponds to a unique vector value, and the number of rows should match the number of rows in the \code{coords} matrix.
#'
#'
#' @section Validity:
#' An object of class \code{ROIVec} is considered valid if:
#' - The \code{coords} slot is a matrix with 3 columns.
#' - The \code{.Data} slot is a matrix.
#' - The number of rows in the \code{.Data} matrix is equal to the number of rows in the \code{coords} matrix.
#'
#' @name ROIVec-class
#' @export
setClass("ROIVec",
         contains=c("ROICoords", "matrix"),
         validity = function(object) {
           if (ncol(object@coords) != 3) {
             stop("coords slot must be a matrix with 3 columns")
           }

           if (ncol(object@.Data) != nrow(object@coords)) {
             stop("'ncol(object)' must equal 'nrow(coords)'")
           }
         })


#' ROIVecWindow
#'
#' A class representing a spatially windowed, vector-valued volumetric region of interest (ROI) in a brain image.
#'
#' @slot coords A \code{matrix} containing the 3D coordinates of the voxels within the ROI. Each row represents a voxel coordinate as (x, y, z).
#' @slot .Data A \code{matrix} containing the data values associated with each voxel in the ROI. Each row corresponds to a unique vector value, and the number of rows should match the number of rows in the \code{coords} matrix.
#' @slot parent_index An \code{integer} representing the 1D index of the center voxel in the parent space.
#' @slot center_index An \code{integer} representing the location in the coordinate matrix of the center voxel in the window.
#'
#'
#' @section Validity:
#' An object of class \code{ROIVecWindow} is considered valid if:
#' - The \code{coords} slot is a matrix with 3 columns.
#' - The \code{.Data} slot is a matrix.
#' - The number of rows in the \code{.Data} matrix is equal to the number of rows in the \code{coords} matrix.
#'
#' @name ROIVecWindow-class
#' @export
setClass("ROIVecWindow",
         representation=representation(parent_index="integer", center_index="integer"),
         contains=c("ROIVec"),
         validity = function(object) {
           if (ncol(object@coords) != 3) {
             stop("coords slot must be a matrix with 3 columns")
           }

           if (ncol(object@.Data) != nrow(object@coords)) {
             stop("'ncol(object)' must equal 'nrow(coords)'")
           }
         })


#' Kernel
#'
#' A class representing an image kernel for image processing, such as convolution or filtering operations in brain images.
#'
#' @slot width A \code{numeric} value representing the width of the kernel in voxels. The width is typically an odd number to maintain symmetry.
#' @slot weights A \code{numeric} vector containing the weights associated with each voxel in the kernel.
#' @slot voxels A \code{matrix} containing the relative voxel coordinates of the kernel. Each row represents a voxel coordinate as (x, y, z).
#' @slot coords A \code{matrix} containing the relative real-world coordinates of the kernel, corresponding to the voxel coordinates.
#'
#'
#' @name Kernel-class
#' @export
setClass("Kernel", representation(width="numeric", weights="numeric", voxels="matrix", coords="matrix"))



#' NeuroBucket
#'
#' a four-dimensional image that consists of a sequence of labeled image volumes backed by a list
#'
#' @rdname NeuroBucket-class
#' @slot labels the names of the sub-volumes contained in the bucket
#' @slot data a list of \code{\linkS4class{NeuroVol}} instances with names corresponding to volume labels
#' @export
#' @importFrom purrr map_lgl
setClass("NeuroBucket",
         representation=representation(labels="character", data="list"),
         validity = function(object) {
           if (any(map_lgl(object@data, function(obj) !is(obj, "NeuroVol")))) {
             stop("all elements of data list must be of type `NeuroVol`")
           } else {
             TRUE
           }
         },
         contains=c("NeuroVec"))


#' ColumnReader
#'
#' A class that supports reading of data from a matrix-like storage format, such as a file or a database, in a column-wise manner.
#'
#' @slot nrow An \code{integer} representing the number of rows in the matrix-like storage.
#' @slot ncol An \code{integer} representing the number of columns in the matrix-like storage.
#' @slot reader A \code{function} that takes a set of column indices as input and returns a \code{matrix} containing the requested columns from the storage.
#'
#' @name ColumnReader-class
#' @export
setClass("ColumnReader", representation=
          representation(nrow="integer", ncol="integer", reader="function"))
          #contains=c("function")


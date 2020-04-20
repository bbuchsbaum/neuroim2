#' @include all_generic.R
#' @import methods
#' @import Matrix
NULL

setOldClass(c("file", "connection"))
setOldClass(c("gzfile", "connection"))
setOldClass("environment")
setOldClass("mmap")
setOldClass("H5File")
#setOldClass("FBM")

setClass("ArrayLike5D")
setClass("ArrayLike4D")
setClass("ArrayLike3D")

#' NamedAxis
#'
#' This class represents an axis with a name attribute
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
#' A two-dimensional axis set
#'
#' @rdname AxisSet2D-class
#' @slot j the second axis
#' @export
setClass("AxisSet2D", representation(j="NamedAxis"), contains=c("AxisSet1D"))

#' AxisSet3D
#'
#' A three-dimensional axis set
#'
#' @rdname AxisSet3D-class
#' @slot k the third axis
#' @export
setClass("AxisSet3D", representation(k="NamedAxis"),contains=c("AxisSet2D"))

#' AxisSet4D
#'
#' A four-dimensional axis set
#' @name AxisSet4D-class
#' @slot l the fourth axis
#' @export
setClass("AxisSet4D", representation(l="NamedAxis"), contains=c("AxisSet3D"))

#' AxisSet5D
#'
#' A five-dimensional axis set
#'
#' @name AxisSet5D-class
#' @slot m the fifth axis
#' @export
setClass("AxisSet5D", representation(m="NamedAxis"), contains=c("AxisSet4D"))

#' FileFormatDescriptor
#'
#' This class represents a neuroimaging file format
#'
#' @rdname FileFormat-class
#'
#' @slot file_format the name of the file format (e.g. NIfTI)
#' @slot header_encoding the file encoding of the header file (e.g. 'raw' for binary, 'gzip' for gz compressed')
#' @slot header_extension the file extension for the header file (e.g. 'nii' for NIfTI single files)
#' @slot data_encoding the file encoding for the data file
#' @slot data_extension the file extension for the data file (e.g. 'nii' for NIfTI single files)
#' @exportClass FileFormat
setClass("FileFormat",
         representation=
           representation(file_format="character",
                          header_encoding="character",
                          header_extension="character",
                          data_encoding="character",
                          data_extension="character")
)

#' NIFTIFormat
#'
#' This class supports the NIFTI file format
#'
#' @rdname NIFTIFileDescriptor-class
#' @export
setClass("NIFTIFormat", contains=c("FileFormat"))


#' AFNIFormat
#'
#' This class supports the AFNI file format
#' @rdname AFNIDescriptor-class
#' @export
setClass("AFNIFormat", contains=c("FileFormat"))

#' H5Format
#'
#' This class supports the AFNI file format
#' @keyword internal
setClass("H5Format", contains=c("FileFormat"))


#' MetaInfo
#'
#' This class contains meta information for an image data type
#' @rdname MetaInfo-class
#' @slot data_type the data type code, e.g. FLOAT
#' @slot dims image dimensions
#' @slot spatial_axes image axes for spatial dimensions (x,y,z)
#' @slot additional_axes axes for dimensions > 3 (e.g. time, color band, direction)
#' @slot spacing voxel dimensions
#' @slot origin coordinate origin
#' @slot label name(s) of images
#' @export
setClass("MetaInfo",
         representation=
           representation(
             data_type="character",
             dims="numeric",
             spatial_axes="AxisSet3D",
             additional_axes="AxisSet",
             spacing="numeric",
             origin="numeric",
             label="character"))

#' FileMetaInfo
#'
#' This class contains meta information from an image data file
#'
#' @rdname FileMetaInfo-class
#' @slot header_file name of the file containing meta information
#' @slot data_file name of the file containing data
#' @slot descriptor descriptor of image file format
#' @slot endian byte order of data ('little' or 'big')
#' @slot data_offset the number of bytes preceding the start of image data in data file
#' @slot bytes_per_element number of bytes per element
#' @slot intercept constant value added to image -- multiple values allowed (must equal numer of sub-images)
#' @slot slope image multiplier -- multiple values allowed (must equal numer of sub-images)
#' @slot header a list of format specific attributes
#' @export
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

#' NIFTIMetaInfo
#'
#' This class contains meta information for a NIfTI image file
#' @rdname FileMetaInfo-class
#' @slot nifti_header a \code{list} of attributes specific to the NIfTI file format
#' @export
setClass("NIFTIMetaInfo",
         representation=(nifti_header="list"),
         contains=c("FileMetaInfo"))

#' AFNIMetaInfo
#'
#' This class contains meta information for a AFNI image file.
#'
#' @rdname FileMetaInfo-class
#' @slot afni_header a list of attributes specific to the AFNI file format
#' @slot afni_header a \code{list} of attributes specific to the AFNI file format
#' @export
setClass("AFNIMetaInfo",
         representation(afni_header="list"),
         contains=c("FileMetaInfo"))



#' FileSource
#'
#' Base class for representing a data source for images. The purpose of this class is to provide a layer in between
#' low level IO and image loading functionality.
#'
#' @rdname FileSource-class
#' @slot meta_info meta information for the data source
#' @exportClass FileSource
setClass("FileSource", representation(meta_info="FileMetaInfo"))



#' NeuroVolSource
#'
#' A class is used to produce a \code{\linkS4class{NeuroVol}} instance
#' @rdname NeuroVolSource-class
#' @slot index the index of the volume to be read -- must be of length 1.
#' @exportClass NeuroVolSource
setClass("NeuroVolSource", representation(index="integer"), contains="FileSource")



#' NeuroVecSource
#'
#' A class that is used to produce a \code{\linkS4class{NeuroVec}} instance
#' @rdname NeuroVecSource-class
#' @slot indices the index vector of the volumes to be loaded
#' @export
setClass("NeuroVecSource", representation(indices="integer"), contains="FileSource")


#' H5NeuroVecSource
#'
#' A class that is used to produce a \code{\linkS4class{H5NeuroVecSource}} instance
#'
#' @rdname H5NeuroVecSource-class
#' @slot file_name the name of the hdf5 file.
setClass("H5NeuroVecSource", representation(file_name="character"))


#' LatentNeuroVecSource
#'
#' A class that is used to produce a \code{\linkS4class{LatentNeuroVecSource}} instance
#'
#' @rdname LatentNeuroVecSource-class
#' @slot file_name the name of the hdf5 file.
setClass("LatentNeuroVecSource", representation(file_name="character"))



#' BinaryReader
#'
#' This class supports reading of bulk binary data from a connection
#' @rdname BinaryReader-class
#' @slot input the binary input connection
#' @slot byte_offset the number of bytes to skip at the start of input
#' @slot data_type the dataType of the binary Elements
#' @slot bytes_per_element number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian endianness of binary input connection
#' @export
setClass("BinaryReader",
           representation(input="connection",
                          byte_offset="numeric",
                          data_type="character",
                          bytes_per_element="integer",
                          endian="character"))



#' BinaryWriter
#'
#' This class supports writing of bulk binary data to a connection
#' @rdname BinaryWriter-class
#' @slot output the binary output connection
#' @slot byte_offset the number of bytes to skip at the start of input
#' @slot data_type the dataType of the binary Elements
#' @slot bytes_per_element number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian endianness of binary output connection
#' @export
setClass("BinaryWriter",
           representation(output="connection",
                          byte_offset="numeric",
                          data_type="character",
                          bytes_per_element="integer",
                          endian="character"))

#' NeuroSpace
#'
#' This class represents the geometry of a brain image
#' @rdname NeuroSpace-class
#' @slot dim the grid dimensions of the image
#' @slot origin the coordinates of the spatial origin
#' @slot spacing the dimensions (in mm) of the grid units (voxels)
#' @slot axes the set of named spatial axes in the untransformed native grid space.
#' @slot trans an affine transformation matrix that moves from grid -> real world coordinates
#' @slot inverse an inverse matrix that moves from real world -> grid coordinates
#' @export
#'
# TODO add 'ref_space' e.g. the name of the coordinate reference space (e.g. LPI)?
setClass("NeuroSpace",
        representation(dim = "integer", origin = "numeric", spacing = "numeric",
                       axes="AxisSet", trans="matrix", inverse="matrix"),

         validity = function(object) {
           dim <- object@dim
           if (length(dim) < length(object@spacing)) {
             return("dim slot must be of same length as spacing slot")
           }
           if (length(dim) < length(object@origin)) {
             return("@dim slot must be of same length as origin slot")
           }
           if (length(dim) < ndim(object@axes)) {
             return("@dim slot must be of same length as number of axes in AxisSet")
           }

           if (any(dim) < 0) {
             return("@dim slot must contain non-negative values")
           }
         })

#' NeuroObj
#' Base class for all data objects with a cartesion spatial represenetation
#'
#' @slot space the geometry of the image object
setClass("NeuroObj", representation(space="NeuroSpace"))

#' NeuroSlice
#'
#' A two-dimensional brain image.
#'
#' @rdname NeuroSlice-class
#' @export
setClass("NeuroSlice", contains=c("array", "NeuroObj"))


#' NeuroVol
#'
#' Base class for image representing 3D volumetric data.
#' @rdname NeuroVol-class
#' @export
setClass("NeuroVol", contains="NeuroObj")



#' DenseNeuroVol
#'
#' Three-dimensional brain image, backed by an \code{array}
#' @rdname DenseNeuroVol-class
#' @export
setClass("DenseNeuroVol", contains=c("NeuroVol", "array"))


#' SparseNeuroVol
#'
#' Three-dimensional brain image, backed by a \code{sparseVector} for \code{Matrix} package
#' @slot data a \code{sparseVector} instance
#' @importFrom Matrix sparseVector
#' @rdname SparseNeuroVol-class
#' @export
setClass("SparseNeuroVol",
         representation=representation(data="sparseVector"),
         contains=c("NeuroVol", "ArrayLike3D"))


#' LogicalNeuroVol
#'
#' Three-dimensional brain image where all values are either TRUE or FALSE
#' @rdname LogicalNeuroVol-class
#' @export
setClass("LogicalNeuroVol", contains=c("DenseNeuroVol"))

#' ClusteredNeuroVol
#'
#' Three-dimensional brain image that is divided into N disjoint partitions
#'
#' @slot mask the \code{logical} mask indicating the spatial domain of the set of clusters
#' @slot clusters an integer index indicating the cluster number for each voxel in the mask
#' @slot label_map a list mapping from name to cluster id
#' @slot cluster_map an \code{environment} mapping from cluster id to the set of 1D spatial indices
#' @rdname ClusteredNeuroVol-class
#' @export
setClass("ClusteredNeuroVol",
         representation=representation(
                                       mask="LogicalNeuroVol",
                                       clusters="integer",
                                       label_map="list",
                                       cluster_map="environment"),
         contains=c("SparseNeuroVol"))


#' IndexLookupVol
#'
#' Three-dimensional brain image that can be used as a map between 1D grid indices and a table of values
#' Currently used in the \code{\linkS4class{SparseNeuroVec}} class.
#'
#' @rdname IndexLookupVol-class
#' @export
setClass("IndexLookupVol",
         representation=
           representation(space="NeuroSpace", indices="integer", map="integer"),
         contains=c("NeuroVol"))



#' NeuroVec
#'
#' Four-dimensional brain image
#'
#' @rdname NeuroVec-class
#' @export
setClass("NeuroVec", contains="NeuroObj")


#' NeuroHyperVec
#'
#' Five-dimensional brain image
#'
#' @rdname NeuroHyperVec-class
#' @export
setClass("NeuroHyperVec", contains="NeuroObj",
         representation(vecs="list"),

    validity = function(object) {
      assert_that(all(purrr::map_lgl(object@vecs, ~ inherits(., "NeuroVec"))))
      dimlist <- purrr::map(object@vecs, ~ dim(.)[1:3])
      d4 <- purrr::map(object@vecs, ~ dim(.)[4])
      splist <- purrr::map(object@vecs, ~ spacing(.))
      assert_that(all(purrr::map_lgl(dimlist, ~ all(dimlist[[1]] == .))))
      assert_that(all(purrr::map_lgl(splist, ~ all(splist[[1]] == .))))
      assert_that(all(purrr::map_lgl(d4, ~ all(d4[[1]] == .))))
    }
)


#' DenseNeuroVec
#'
#' Four-dimensional brain image, backed by an array
#' @name DenseNeuroVec-class
#' @export
setClass("DenseNeuroVec",  contains=c("NeuroVec", "array"))



#' MappedNeuroVec
#'
#' Four-dimensional brain image, backed by an memory-mapped file.
#'
#' @name MappedNeuroVec-class
#' @export
setClass("MappedNeuroVec",  representation(filemap="mmap", offset="integer"), contains=c("NeuroVec", "ArrayLike4D"))


#' AbstractSparseNeuroVec
#'
#' @slot mask the mask defining the sparse domain
#' @slot map instance of class \code{\linkS4class{IndexLookupVol}} is used to map between spatial and index/row coordinates
setClass("AbstractSparseNeuroVec",
         representation(mask="LogicalNeuroVol", map="IndexLookupVol"),
         contains=c("NeuroVec", "ArrayLike4D"))


#' SparseNeuroVec
#'
#' a sparse four-dimensional brain image, backed by a \code{matrix}, where each column represents
#' a non-zero vector spanning the fourth dimension (e.g. time), and defined by a volumetric mask.
#'
#' @slot data the matrix of series, where rows span across voxel space and columns span the fourth dimensions
#'
#' @rdname SparseNeuroVec-class
setClass("SparseNeuroVec",
         representation(data="matrix"),
         contains=c("NeuroVec", "AbstractSparseNeuroVec", "ArrayLike4D"))

#' BigNeuroVec
#'
#' a sparse four-dimensional brain image that is backed by a disk-based big-matrix
#'
#' @rdname BigNeuroVec-class
#'
#' @slot data an instance of class \code{FBM} from \code{bigstatsr} package which contains time-series data.
#' @importFrom bigstatsr FBM
setClass("BigNeuroVec",
         representation(data="FBM"),
         contains=c("NeuroVec", "AbstractSparseNeuroVec", "ArrayLike4D"))





#' FileBackedNeuroVec
#'
#' a four-dimensional brain image that is read in to memory "on demand" using memory-mapped file access.
#'
#' @rdname FileBackedNeuroVec-class
#' @slot meta the file meta information of type \code{\linkS4class{FileMetaInfo}}
setClass("FileBackedNeuroVec",
         representation(meta="FileMetaInfo"),
         contains=c("NeuroVec", "ArrayLike4D"))



#' H5NeuroVol
#'
#' a three-dimensional brain image backed by an HDF5 dataset
#'
#' @rdname H5NeuroVol-class
# @importClassesFrom hdf5r H5File
setClass("H5NeuroVol",
         representation(h5obj="H5File"),
         contains=c("NeuroVol", "ArrayLike3D"))


setClass("H5NeuroVec",
         representation(obj="H5File"),
         contains=c("NeuroVec", "ArrayLike4D"),
         validity = function(object) {
            TRUE
         })



#' NeuroVecSeq
#'
#' A concatenated sequence of \code{\linkS4class{NeuroVec}} instances.
#'
#' @rdname NeuroVecSeq-class
#' @slot vecs the sequences of \code{NeuroVec} instances
#' @slot lens the number of volumes in each \code{NeuroVec} sequence
#' @export
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


#' LatentNeuroVec
#'
#' a class that stores a represents a 4-dimensional array as a set of basis functions (dictionary) and
#' corresponding set of loadings (coefficients).
#'
#' @rdname LatentNeuroVec-class
#'
#' @slot basis the matrix of bases, were each column is a basis vector.
#' @slot loadings the \code{sparseMatrix} of loadings
#' @slot offset an offset vector
#' @export
setClass("LatentNeuroVec",
         representation(basis="Matrix",
                        loadings="Matrix",
                        offset="numeric"),
         contains=c("NeuroVec", "AbstractSparseNeuroVec"))



#' SparseNeuroVecSource
#'
#' A class that is used to produce a \code{\linkS4class{SparseNeuroVec}} instance
#'
#' @rdname SparseNeuroVecSource-class
#' @slot mask the subset of voxels that will be stored in memory
#' @export
setClass("SparseNeuroVecSource", representation(mask="LogicalNeuroVol"), contains=c("NeuroVecSource"))



#' MappedNeuroVecSource
#'
#' A class that is used to produce a \code{\linkS4class{MappedNeuroVec}} instance
#'
#' @rdname MappedNeuroVecSource-class
#' @export
setClass("MappedNeuroVecSource", contains=c("NeuroVecSource"))



setClassUnion("numericOrMatrix", c("numeric", "matrix"))

#' ROI
#'
#' Base marker class for a region of interest (ROI)
#'
#' @export
setClass("ROI", contains="VIRTUAL")

setClass("ROICoords",
         representation=representation(space="NeuroSpace", coords="matrix"),
         contains=c("ROI"))

#' ROIVol
#'
#' A class that represents a volumetric region of interest
#'
#' @rdname ROIVol-class
#' @exportClass ROIVol
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
#' A class that represents a spatially windowed volumetric region of interest
#'
#' @rdname ROIVolWindow-class
#' @slot parent_index the 1D index of the center voxel in the parent space.
#' @slot center_index the location in the coordinate matrix of the center voxel in the window
#' @export
setClass("ROIVolWindow",
         contains=c("ROIVol"),
         representation=representation(parent_index="integer", center_index="integer"),
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
#' A class that represents a vector-valued volumetric region of interest
#'
#' @rdname ROIVec-class
#' @exportClass ROIVec
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
#' A class that represents a vector-valued volumetric region of interest
#'
#' @slot parent_index the 1D index of the center voxel in the parent space.
#' @slot center_index the location in the coordinate matrix of the center voxel in the window
#' @rdname ROIVecWindow-class
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
#' A class representing an image kernel
#'
#' @rdname Kernel-class
#' @slot width the width in voxels of the kernel
#' @slot weights the kernel weights
#' @slot voxels the relative voxel coordinates of the kernel
#' @slot coords the relative real coordinates of the kernel
#' @export
setClass("Kernel", representation(width="numeric", weights="numeric", voxels="matrix", coords="matrix"))



#' NeuroBucket
#'
#' a four-dimensional image that conists of a sequence of labeled image volumes backed by a list
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
#' This class supports reading of data froma matrix-like storage format
#' @rdname ColumnReader-class
#' @slot nrow the number of rows
#' @slot ncol the number of columns
#' @slot reader a function that takes a set of column indices and returns a \code{matrix}
#' @export
setClass("ColumnReader", representation=
          representation(nrow="integer", ncol="integer", reader="function"))
          #contains=c("function")


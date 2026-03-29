pkgname <- "neuroim2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('neuroim2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AxisSet2D-class")
### * AxisSet2D-class

flush(stderr()); flush(stdout())

### Name: AxisSet2D-class
### Title: AxisSet2D
### Aliases: AxisSet2D-class

### ** Examples

# Create an AxisSet2D object
axis1 <- new("NamedAxis", axis = "x", direction = 1)
axis2 <- new("NamedAxis", axis = "y", direction = 1)
axisSet2D <- new("AxisSet2D", i = axis1, j = axis2, ndim = 2L)



cleanEx()
nameEx("AxisSet3D-class")
### * AxisSet3D-class

flush(stderr()); flush(stdout())

### Name: AxisSet3D-class
### Title: AxisSet3D Class
### Aliases: AxisSet3D-class

### ** Examples

# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)

# Create an AxisSet3D object
axis_set_3d <- new("AxisSet3D", i = x_axis, j = y_axis, k = z_axis, ndim = 3L)




cleanEx()
nameEx("AxisSet4D-class")
### * AxisSet4D-class

flush(stderr()); flush(stdout())

### Name: AxisSet4D-class
### Title: AxisSet4D Class
### Aliases: AxisSet4D-class

### ** Examples

# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)
t_axis <- new("NamedAxis", axis = "t", direction = 1)

# Create an AxisSet4D object
axis_set_4d <- new("AxisSet4D", i = x_axis, j = y_axis, k = z_axis,
                   l = t_axis, ndim = 4L)




cleanEx()
nameEx("AxisSet5D-class")
### * AxisSet5D-class

flush(stderr()); flush(stdout())

### Name: AxisSet5D-class
### Title: AxisSet5D Class
### Aliases: AxisSet5D-class

### ** Examples

# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)
t_axis <- new("NamedAxis", axis = "t", direction = 1)
v_axis <- new("NamedAxis", axis = "v", direction = 1)

# Create an AxisSet5D object
axis_set_5d <- new("AxisSet5D", i = x_axis, j = y_axis, k = z_axis,
                   l = t_axis, m = v_axis, ndim = 5L)




cleanEx()
nameEx("BigNeuroVec-methods")
### * BigNeuroVec-methods

flush(stderr()); flush(stdout())

### Name: BigNeuroVec
### Title: Create a Memory-Mapped Neuroimaging Vector
### Aliases: BigNeuroVec

### ** Examples





cleanEx()
nameEx("BinaryReader")
### * BinaryReader

flush(stderr()); flush(stdout())

### Name: BinaryReader
### Title: Create Binary Reader Object
### Aliases: BinaryReader

### ** Examples




cleanEx()
nameEx("BinaryWriter")
### * BinaryWriter

flush(stderr()); flush(stdout())

### Name: BinaryWriter
### Title: Create Binary Writer Object
### Aliases: BinaryWriter

### ** Examples




cleanEx()
nameEx("ClusteredNeuroVec")
### * ClusteredNeuroVec

flush(stderr()); flush(stdout())

### Name: ClusteredNeuroVec
### Title: ClusteredNeuroVec: Cluster-aware 4D neuroimaging data
### Aliases: ClusteredNeuroVec

### ** Examples

# Create synthetic 4D data (10x10x10 volume, 20 timepoints)
sp4 <- NeuroSpace(c(10,10,10,20), c(1,1,1))
arr <- array(rnorm(10*10*10*20), dim=c(10,10,10,20))
vec <- NeuroVec(arr, sp4)

# Create a mask covering the central region
sp3 <- NeuroSpace(c(10,10,10), c(1,1,1))
mask_arr <- array(FALSE, dim=c(10,10,10))
mask_arr[3:8, 3:8, 3:8] <- TRUE
mask <- LogicalNeuroVol(mask_arr, sp3)

# Assign voxels to 5 random clusters
n_voxels <- sum(mask_arr)
clusters <- sample(1:5, n_voxels, replace=TRUE)
cvol <- ClusteredNeuroVol(mask, clusters)

# Create clustered representation
cv <- ClusteredNeuroVec(vec, cvol)

# Access like a regular NeuroVec
vol_t1 <- cv[,,,1]  # 3D volume at time 1
ts <- series(cv, 5, 5, 5)  # time-series at voxel (5,5,5)

# Get cluster time-series matrix
cluster_matrix <- as.matrix(cv)  # T x K matrix
dim(cluster_matrix)  # 20 x 5



cleanEx()
nameEx("ClusteredNeuroVol-class")
### * ClusteredNeuroVol-class

flush(stderr()); flush(stdout())

### Name: ClusteredNeuroVol-class
### Title: ClusteredNeuroVol Class
### Aliases: ClusteredNeuroVol-class ClusteredNeuroVol

### ** Examples


# Create a simple clustered brain volume
dim <- c(10L, 10L, 10L)
mask_data <- array(rep(c(TRUE, FALSE), 500), dim)
mask <- new("LogicalNeuroVol", .Data = mask_data,
            space = NeuroSpace(dim = dim, origin = c(0,0,0), spacing = c(1,1,1)))

clusters <- as.integer(runif(sum(mask_data)) * 5)+1
label_map <- list("Cluster1" = 1, "Cluster2" = 2, "Cluster3" = 3,
                  "Cluster4" = 4, "Cluster5" = 5)

cluster_map <- list()
for (i in 1:5) {
  cluster_map[[as.character(i)]] <- which(clusters == i)
}

clustered_vol <- ClusteredNeuroVol(
                     mask = mask,
                     clusters = clusters,
                     label_map = label_map)



# Create a simple space and volume
space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
vol <- NeuroVol(vol_data, space)

# Create a binary mask (e.g., values > 0)
mask_data <- vol_data > 0
mask_vol <- LogicalNeuroVol(mask_data, space)

# Get coordinates of masked voxels
mask_idx <- which(mask_data)
coords <- index_to_coord(mask_vol, mask_idx)

# Cluster the coordinates into 10 groups
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 10)

# Create the clustered volume
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Print information about the clusters
print(clustered_vol)



cleanEx()
nameEx("ColumnReader")
### * ColumnReader

flush(stderr()); flush(stdout())

### Name: ColumnReader
### Title: Create Column Reader Object
### Aliases: ColumnReader

### ** Examples


reader_func <- function(cols) {
  matrix(rnorm(100 * length(cols)), 100, length(cols))
}
col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)




cleanEx()
nameEx("DenseNeuroVec-class")
### * DenseNeuroVec-class

flush(stderr()); flush(stdout())

### Name: DenseNeuroVec-class
### Title: DenseNeuroVec Class
### Aliases: DenseNeuroVec-class DenseNeuroVec

### ** Examples


# Create a simple 4D brain image
data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
space <- NeuroSpace(dim = c(64, 64, 32,10), origin = c(0, 0, 0), spacing = c(3, 3, 4))
dense_vec <- new("DenseNeuroVec", .Data = data, space = space)

# Access dimensions
dim(dense_vec)

# Extract a single 3D volume
first_volume <- dense_vec[[1]]


# Create a simple 4D brain image
dim <- c(64, 64, 32, 10)  # 64x64x32 volume with 10 time points
data <- array(rnorm(prod(dim)), dim)
space <- NeuroSpace(dim, spacing = c(3, 3, 4))

# Create a DenseNeuroVec object
dense_vec <- DenseNeuroVec(data = data, space = space, label = "Example")
print(dense_vec)

# Create from a matrix (voxels x time)
mat_data <- matrix(rnorm(prod(dim)), nrow = prod(dim[1:3]), ncol = dim[4])
dense_vec_mat <- DenseNeuroVec(data = mat_data, space = space)
print(dense_vec_mat)




cleanEx()
nameEx("DenseNeuroVol-class")
### * DenseNeuroVol-class

flush(stderr()); flush(stdout())

### Name: DenseNeuroVol-class
### Title: DenseNeuroVol Class
### Aliases: DenseNeuroVol-class DenseNeuroVol

### ** Examples

# Create a simple 3D brain volume
vol_data <- array(rnorm(64*64*64), c(64, 64, 64))
vol_space <- NeuroSpace(dim=c(64L, 64L, 64L), origin=c(0, 0, 0), spacing=c(1, 1, 1))
brain_vol <- new("DenseNeuroVol", .Data=vol_data, space=vol_space)




cleanEx()
nameEx("FileBackedNeuroVec-class")
### * FileBackedNeuroVec-class

flush(stderr()); flush(stdout())

### Name: FileBackedNeuroVec-class
### Title: FileBackedNeuroVec Class
### Aliases: FileBackedNeuroVec-class

### ** Examples

# Load example 4D image file included with package
file_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
fbvec <- FileBackedNeuroVec(file_path)

# Get dimensions of the image
dim(fbvec)

# Extract first volume
vol1 <- sub_vector(fbvec, 1)

# Extract multiple volumes
vols <- sub_vector(fbvec, 1:2)




cleanEx()
nameEx("FileBackedNeuroVec")
### * FileBackedNeuroVec

flush(stderr()); flush(stdout())

### Name: FileBackedNeuroVec
### Title: Create a File-Backed Neuroimaging Vector
### Aliases: FileBackedNeuroVec

### ** Examples


# Create a file-backed vector from a NIFTI file
fbvec <- FileBackedNeuroVec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Access specific volumes without loading entire dataset
first_vol <- sub_vector(fbvec, 1)





cleanEx()
nameEx("FileFormat-class")
### * FileFormat-class

flush(stderr()); flush(stdout())

### Name: FileFormat-class
### Title: FileFormat Class
### Aliases: FileFormat-class

### ** Examples

# Create a FileFormat object for NIfTI format
nifti_format <- new("FileFormat",
                    file_format = "NIfTI",
                    header_encoding = "raw",
                    header_extension = "nii",
                    data_encoding = "raw",
                    data_extension = "nii")





cleanEx()
nameEx("IndexLookupVol-class")
### * IndexLookupVol-class

flush(stderr()); flush(stdout())

### Name: IndexLookupVol-class
### Title: IndexLookupVol Class
### Aliases: IndexLookupVol-class IndexLookupVol

### ** Examples

# Create a NeuroSpace object
space <- NeuroSpace(dim = c(2L, 2L, 2L), origin = c(0, 0, 0), spacing = c(1, 1, 1))

# Create a 3D mask
mask <- array(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE), dim = c(2, 2, 2))

# Create indices and map for the IndexLookupVol
indices <- which(mask)
map <- seq_along(indices)

# Create an IndexLookupVol object
ilv <- IndexLookupVol(space = space, indices = as.integer(indices))

# Access the indices
print(ilv@indices)

# Access the map
print(ilv@map)


# Create a 64x64x64 space
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))

# Create a lookup volume with random indices
indices <- sample(1:262144, 10000)  # Select 10000 random voxels
ilv <- IndexLookupVol(space, indices)

# Look up coordinates for specific indices
coords <- coords(ilv, indices[1:10])





cleanEx()
nameEx("Kernel")
### * Kernel

flush(stderr()); flush(stdout())

### Name: Kernel
### Title: Create a Kernel object from a function of distance from kernel
###   center
### Aliases: Kernel

### ** Examples

kdim <- c(3, 3, 3)
vdim <- c(1, 1, 1)
k <- Kernel(kerndim = kdim, vdim = vdim, FUN = dnorm, sd = 1)



cleanEx()
nameEx("Logic-methods")
### * Logic-methods

flush(stderr()); flush(stdout())

### Name: Logic-methods
### Title: Logic Operations for Neuroimaging Volumes
### Aliases: Logic-methods Logic,DenseNeuroVol,DenseNeuroVol-method
###   Logic,SparseNeuroVol,SparseNeuroVol-method
###   Logic,SparseNeuroVol,NeuroVol-method
###   Logic,NeuroVol,SparseNeuroVol-method Logic,NeuroVol,logical-method
###   Logic,logical,NeuroVol-method

### ** Examples

sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
v1 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
v2 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
intersection <- v1 & v2
union_mask  <- v1 | v2




cleanEx()
nameEx("LogicalNeuroVol-class")
### * LogicalNeuroVol-class

flush(stderr()); flush(stdout())

### Name: LogicalNeuroVol-class
### Title: LogicalNeuroVol Class
### Aliases: LogicalNeuroVol-class LogicalNeuroVol

### ** Examples

# Create a simple logical brain volume (e.g., a mask)
dim <- c(64L, 64L, 64L)
mask_data <- array(sample(c(TRUE, FALSE), prod(dim), replace = TRUE), dim)
mask_space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
brain_mask <- new("LogicalNeuroVol", .Data = mask_data, space = mask_space)

# Check the proportion of TRUE voxels
true_proportion <- sum(brain_mask) / prod(dim(brain_mask))
print(paste("Proportion of TRUE voxels:", true_proportion))

# Load an example brain mask
brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Convert the brain mask to a LogicalNeuroVol
logical_vol <- LogicalNeuroVol(brain_mask, space(brain_mask))




cleanEx()
nameEx("MappedNeuroVec-class")
### * MappedNeuroVec-class

flush(stderr()); flush(stdout())

### Name: MappedNeuroVec-class
### Title: MappedNeuroVec Class
### Aliases: MappedNeuroVec-class MappedNeuroVec

### ** Examples


# Create a MappedNeuroVec object (pseudo-code)
file_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
mapped_vec <- MappedNeuroVec(file_path)

# Access a subset of the data
subset <- mapped_vec[,,, 1:2]


# Create mapped vector from NIFTI file
mvec <- MappedNeuroVec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Extract first volume
vol1 <- mvec[[1]]

# Get dimensions
dim(mvec)

# Access specific timepoint
timepoint <- mvec[, , , 2]





cleanEx()
nameEx("MappedNeuroVecSource-class")
### * MappedNeuroVecSource-class

flush(stderr()); flush(stdout())

### Name: MappedNeuroVecSource-class
### Title: MappedNeuroVecSource Class
### Aliases: MappedNeuroVecSource-class MappedNeuroVecSource

### ** Examples

# Create a MappedNeuroVecSource
mapped_source <- new("MappedNeuroVecSource")





cleanEx()
nameEx("MetaInfo-class")
### * MetaInfo-class

flush(stderr()); flush(stdout())

### Name: MetaInfo-class
### Title: MetaInfo Class
### Aliases: MetaInfo-class

### ** Examples

# Create a MetaInfo object
meta_info <- new("MetaInfo",
                 data_type = "FLOAT",
                 dims = c(64, 64, 32, 100),
                 spatial_axes = new("AxisSet3D"),
                 additional_axes = new("AxisSet"),
                 spacing = c(3, 3, 4),
                 origin = c(0, 0, 0),
                 label = "fMRI_run1")




cleanEx()
nameEx("MetaInfo")
### * MetaInfo

flush(stderr()); flush(stdout())

### Name: MetaInfo
### Title: Create Neuroimaging Metadata Object
### Aliases: MetaInfo

### ** Examples

# Create metadata for 3D structural MRI
meta <- MetaInfo(
  Dim = c(256, 256, 180),
  spacing = c(1, 1, 1),
  data_type = "FLOAT",
  label = "T1w"
)

# Get image dimensions
dim(meta)

# Get transformation matrix
trans(meta)




cleanEx()
nameEx("NIFTIMetaInfo")
### * NIFTIMetaInfo

flush(stderr()); flush(stdout())

### Name: NIFTIMetaInfo
### Title: Create NIFTI Format Metadata Object
### Aliases: NIFTIMetaInfo

### ** Examples


# Read NIFTI header
header <- read_header(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Create format descriptor
fmt <- new("NIFTIFormat",
           file_format = "NIFTI",
           header_encoding = "raw",
           header_extension = "nii",
           data_encoding = "raw",
           data_extension = "nii")

# Create metadata
meta <- NIFTIMetaInfo(fmt, header@header)

# Check dimensions
dim(meta)





cleanEx()
nameEx("NeuroHyperVec-class")
### * NeuroHyperVec-class

flush(stderr()); flush(stdout())

### Name: NeuroHyperVec-class
### Title: NeuroHyperVec Class
### Aliases: NeuroHyperVec-class [,NeuroHyperVec,ANY,ANY,ANY-method
###   [.NeuroHyperVec

### ** Examples


# Create a simple 5D dataset (10x10x10 spatial, 5 trials, 3 features)
dims <- c(10, 10, 10)
space <- NeuroSpace(c(dims, 5, 3))

# Create a sparse mask (20% of voxels)
mask_data <- array(runif(prod(dims)) < 0.2, dims)
mask <- LogicalNeuroVol(mask_data, NeuroSpace(dims))

# Generate random data for active voxels
n_voxels <- sum(mask_data)
data <- array(rnorm(3 * 5 * n_voxels), dim = c(3, 5, n_voxels))  # [features x trials x voxels]

# Create NeuroHyperVec object
hvec <- NeuroHyperVec(data, space, mask)

# Access operations
# Get data for specific voxel across all trials/features
series(hvec, 5, 5, 5)

# Extract a 3D volume for specific trial and feature
hvec[,,,2,1]





cleanEx()
nameEx("NeuroHyperVec")
### * NeuroHyperVec

flush(stderr()); flush(stdout())

### Name: NeuroHyperVec
### Title: Constructor for NeuroHyperVec class
### Aliases: NeuroHyperVec

### ** Examples

# Create a 5D space (10x10x10 spatial, 2 trials, 2 features)
space <- NeuroSpace(c(10,10,10,2,2))

# Create a mask for the spatial dimensions
space3d <- NeuroSpace(c(10,10,10))
mask_data <- array(TRUE, dim=c(10,10,10))  # All voxels active
mask <- LogicalNeuroVol(mask_data, space3d)

# Create data in the format expected by NeuroHyperVec:
# 3D array with dimensions [features x trials x voxels]
n_features <- 2
n_trials <- 2
n_voxels <- sum(mask_data)  # 1000 voxels
data_array <- array(rnorm(n_features * n_trials * n_voxels),
                   dim = c(n_features, n_trials, n_voxels))

# Create the NeuroHyperVec object
hvec <- NeuroHyperVec(data_array, space, mask)




cleanEx()
nameEx("NeuroSlice-class")
### * NeuroSlice-class

flush(stderr()); flush(stdout())

### Name: NeuroSlice-class
### Title: NeuroSlice Class
### Aliases: NeuroSlice-class

### ** Examples

# Create a simple 2D brain slice
slice_data <- matrix(rnorm(64*64), 64, 64)
slice_space <- NeuroSpace(dim=c(64L, 64L), origin=c(0, 0), spacing=c(1, 1))
brain_slice <- new("NeuroSlice", .Data=slice_data, space=slice_space)




cleanEx()
nameEx("NeuroSlice")
### * NeuroSlice

flush(stderr()); flush(stdout())

### Name: NeuroSlice
### Title: NeuroSlice: 2D Neuroimaging Data Container
### Aliases: NeuroSlice

### ** Examples

# Create a 64x64 slice space
slice_space <- NeuroSpace(c(64, 64), spacing = c(2, 2))

# Example 1: Dense slice from matrix
slice_data <- matrix(rnorm(64*64), 64, 64)
dense_slice <- NeuroSlice(slice_data, slice_space)

# Example 2: Dense slice from vector
vec_data <- rnorm(64*64)
vec_slice <- NeuroSlice(vec_data, slice_space)

# Example 3: Sparse slice with specific values
n_points <- 100
sparse_data <- rnorm(n_points)
sparse_indices <- sample(1:(64*64), n_points)
sparse_slice <- NeuroSlice(sparse_data, slice_space, indices = sparse_indices)




cleanEx()
nameEx("NeuroSpace-class")
### * NeuroSpace-class

flush(stderr()); flush(stdout())

### Name: NeuroSpace-class
### Title: NeuroSpace Class
### Aliases: NeuroSpace-class

### ** Examples

# Create a NeuroSpace object
space <- NeuroSpace(dim = c(64L, 64L, 64L),
                    origin = c(0, 0, 0),
                    spacing = c(1, 1, 1))

# Get the dimensions
dim(space)






cleanEx()
nameEx("NeuroSpace")
### * NeuroSpace

flush(stderr()); flush(stdout())

### Name: NeuroSpace
### Title: NeuroSpace: Spatial Reference System for Neuroimaging Data
### Aliases: NeuroSpace

### ** Examples

# Create a standard 3D space (64x64x40 voxels, 2mm isotropic)
space_3d <- NeuroSpace(
  dim = c(64L, 64L, 40L),
  spacing = c(2, 2, 2),
  origin = c(-90, -126, -72)
)

# Check properties
dim(space_3d)           # Image dimensions
spacing(space_3d)       # Voxel sizes
origin(space_3d)        # World-space origin

# Create a 2D slice space
space_2d <- NeuroSpace(
  dim = c(128L, 128L),
  spacing = c(1.5, 1.5),
  origin = c(-96, -96)
)

# Convert between coordinate systems
world_coords <- c(0, 0, 0)
vox_idx <- coord_to_index(space_3d, world_coords)
back_to_world <- index_to_coord(space_3d, vox_idx)




cleanEx()
nameEx("NeuroVec-class")
### * NeuroVec-class

flush(stderr()); flush(stdout())

### Name: NeuroVec-class
### Title: NeuroVec Class
### Aliases: NeuroVec-class NeuroVec

### ** Examples


# Load an example 4D brain image
example_4d_image <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Create a NeuroVec object
neuro_vec <- NeuroVec(data = array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10)),
                      space = NeuroSpace(dim = c(64, 64, 32,10),
                      origin = c(0, 0, 0),
                      spacing = c(3, 3, 4)))


dim(neuro_vec)

# Extract a single 3D volume (e.g., the first time point)
first_volume <- neuro_vec[[1]]


# Load an example 4D brain image
example_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
example_4d_image <- read_vec(example_file)

# Create a DenseNeuroVec object
dense_vec <- NeuroVec(data = example_4d_image@.Data,
                      space = space(example_4d_image))
print(dense_vec)

# Create a SparseNeuroVec object with a mask
mask <- array(runif(prod(dim(example_4d_image)[1:3])) > 0.5,
              dim = dim(example_4d_image)[1:3])
sparse_vec <- NeuroVec(data = example_4d_image@.Data,
                       space = space(example_4d_image),
                       mask = mask)
print(sparse_vec)




cleanEx()
nameEx("NeuroVecSeq")
### * NeuroVecSeq

flush(stderr()); flush(stdout())

### Name: NeuroVecSeq
### Title: NeuroVecSeq: A Container for Sequential NeuroVec Objects
### Aliases: NeuroVecSeq

### ** Examples

# Create some example NeuroVec objects
v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
               space = NeuroSpace(dim = c(5, 5, 5, 2)))
v2 <- NeuroVec(array(1, c(5, 5, 5, 4)),
               space = NeuroSpace(dim = c(5, 5, 5, 4)))
v3 <- NeuroVec(array(2, c(5, 5, 5, 6)),
               space = NeuroSpace(dim = c(5, 5, 5, 6)))

# Combine them into a sequence
vs <- NeuroVecSeq(v1, v2, v3)

# Access properties
length(vs)  # Total time points
vs[[5]]     # Get the 5th volume

# Extract a subsequence
sub_seq <- sub_vector(vs, 1:5)


# Create sample vectors
v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
               space = NeuroSpace(dim = c(5, 5, 5, 2)))
v2 <- NeuroVec(array(0, c(5, 5, 5, 4)),
               space = NeuroSpace(dim = c(5, 5, 5, 4)))

# Combine into sequence
vs <- NeuroVecSeq(v1, v2)
print(vs)





cleanEx()
nameEx("NeuroVol")
### * NeuroVol

flush(stderr()); flush(stdout())

### Name: NeuroVol
### Title: NeuroVol: 3D Neuroimaging Volume Class
### Aliases: NeuroVol

### ** Examples

bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
dat <- array(rnorm(64*64*64), c(64,64,64))
bvol <- NeuroVol(dat,bspace, label="test")



cleanEx()
nameEx("NiftiExtension-class")
### * NiftiExtension-class

flush(stderr()); flush(stdout())

### Name: NiftiExtension-class
### Title: NiftiExtension Class
### Aliases: NiftiExtension-class show,NiftiExtension-method

### ** Examples

# Create a simple comment extension
comment_text <- "This is a test comment"
ext <- NiftiExtension(ecode = 6L, data = comment_text)

# Access the extension code
ext@ecode




cleanEx()
nameEx("NiftiExtension")
### * NiftiExtension

flush(stderr()); flush(stdout())

### Name: NiftiExtension
### Title: Create a NIfTI Extension
### Aliases: NiftiExtension

### ** Examples

# Create a comment extension
ext <- NiftiExtension(ecode = 6L, data = "This is a comment")
ext@ecode
ext@esize

# Create an AFNI extension with XML data
afni_xml <- '<?xml version="1.0"?><AFNI_attributes></AFNI_attributes>'
afni_ext <- NiftiExtension(ecode = 4L, data = afni_xml)




cleanEx()
nameEx("NiftiExtensionCodes")
### * NiftiExtensionCodes

flush(stderr()); flush(stdout())

### Name: NiftiExtensionCodes
### Title: Known NIfTI Extension Codes
### Aliases: NiftiExtensionCodes
### Keywords: datasets

### ** Examples

# Get the code for AFNI extensions
NiftiExtensionCodes["AFNI"]  # Returns 4

# Get the name for a code
names(NiftiExtensionCodes)[NiftiExtensionCodes == 4]  # Returns "AFNI"




cleanEx()
nameEx("NiftiExtensionList-class")
### * NiftiExtensionList-class

flush(stderr()); flush(stdout())

### Name: NiftiExtensionList-class
### Title: NiftiExtensionList Class
### Aliases: NiftiExtensionList-class show,NiftiExtensionList-method

### ** Examples

# Create an empty extension list
ext_list <- new("NiftiExtensionList")

# Create a list with extensions
ext1 <- NiftiExtension(ecode = 6L, data = "Comment 1")
ext2 <- NiftiExtension(ecode = 6L, data = "Comment 2")
ext_list <- new("NiftiExtensionList", list(ext1, ext2))




cleanEx()
nameEx("ROICoords")
### * ROICoords

flush(stderr()); flush(stdout())

### Name: ROICoords
### Title: Create ROI Coordinates Object
### Aliases: ROICoords

### ** Examples

coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
roi_coords <- ROICoords(coords)




cleanEx()
nameEx("ROIVec")
### * ROIVec

flush(stderr()); flush(stdout())

### Name: ROIVec
### Title: Create an instance of class 'ROIVec'
### Aliases: ROIVec

### ** Examples

# Create a NeuroSpace object
vspace <- NeuroSpace(dim = c(5, 5, 5, 10), spacing = c(1, 1, 1))

# Define voxel coordinates for the ROI
coords <- matrix(c(1, 2, 3, 2, 2, 2, 3, 3, 3), ncol = 3)

# Create a data matrix for the ROI
data <- matrix(rnorm(30), nrow = 10, ncol = 3)

# Create a ROIVec object
roi_vec <- ROIVec(vspace, coords, data)



cleanEx()
nameEx("ROIVol")
### * ROIVol

flush(stderr()); flush(stdout())

### Name: ROIVol
### Title: Create ROI Volume Object
### Aliases: ROIVol

### ** Examples

space <- NeuroSpace(c(64,64,64))
coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
data <- c(1.5, 2.5)
roi_vol <- ROIVol(space, coords, data)




cleanEx()
nameEx("SparseNeuroVec-class")
### * SparseNeuroVec-class

flush(stderr()); flush(stdout())

### Name: SparseNeuroVec-class
### Title: SparseNeuroVec Class
### Aliases: SparseNeuroVec-class SparseNeuroVec

### ** Examples


# Create a sparse 4D brain image
mask <- LogicalNeuroVol(array(runif(64*64*32) > 0.7, c(64,64,32)), NeuroSpace(c(64,64,32)))
data <- matrix(rnorm(sum(mask) * 100), nrow=sum(mask), ncol=100)
sparse_vec <- SparseNeuroVec(data=data, mask=mask, space=NeuroSpace(dim=c(64,64,32,100)))

# Access a subset of the data
subset <- sparse_vec[,,, 1:10]


bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)
length(indices(svec)) == sum(mask)



cleanEx()
nameEx("SparseNeuroVecSource-class")
### * SparseNeuroVecSource-class

flush(stderr()); flush(stdout())

### Name: SparseNeuroVecSource-class
### Title: SparseNeuroVecSource Class
### Aliases: SparseNeuroVecSource-class

### ** Examples

# Create a simple mask
mask_data <- array(runif(64*64*32) > 0.7, dim = c(64, 64, 32))
mask <- LogicalNeuroVol(mask_data, space = NeuroSpace(dim = c(64, 64, 32)))

# Create a SparseNeuroVecSource
sparse_source <- new("SparseNeuroVecSource", mask = mask)





cleanEx()
nameEx("SparseNeuroVol-class")
### * SparseNeuroVol-class

flush(stderr()); flush(stdout())

### Name: SparseNeuroVol-class
### Title: SparseNeuroVol Class
### Aliases: SparseNeuroVol-class SparseNeuroVol

### ** Examples


# Create a sparse 3D brain image
dim <- c(64L, 64L, 64L)
space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
sparse_data <- Matrix::sparseVector(x = c(1, 2, 3),
                                    i = c(100, 1000, 10000),
                                    length = prod(dim))
sparse_vol <- new("SparseNeuroVol", space = space, data = sparse_data)
sparse_vol[1000] == 1

data <- 1:10
indices <- seq(1,1000, length.out=10)
bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
sparsevol <- SparseNeuroVol(data,bspace,indices=indices)
densevol <- NeuroVol(data,bspace,indices=indices)
sum(sparsevol) == sum(densevol)




cleanEx()
nameEx("Summary-methods")
### * Summary-methods

flush(stderr()); flush(stdout())

### Name: Summary-methods
### Title: Summary Methods for Neuroimaging Objects
### Aliases: Summary-methods Summary,SparseNeuroVec-method
###   Summary,SparseNeuroVol-method Summary,DenseNeuroVol-method

### ** Examples

# Create a simple volume
vol <- DenseNeuroVol(array(1:27, c(3,3,3)),
                     NeuroSpace(c(3L,3L,3L), c(1,1,1)))
sum(vol)
range(vol)




cleanEx()
nameEx("add_dim-methods")
### * add_dim-methods

flush(stderr()); flush(stdout())

### Name: add_dim
### Title: Add a Dimension to an Object
### Aliases: add_dim add_dim,NeuroSpace,numeric-method

### ** Examples

# Create a NeuroSpace object
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Add a new dimension with size 10
x1 <- add_dim(x, 10)

# Check the new dimension
ndim(x1) == 4
dim(x1)[4] == 10




cleanEx()
nameEx("affine_utils")
### * affine_utils

flush(stderr()); flush(stdout())

### Name: affine_utils
### Title: Affine utility functions
### Aliases: affine_utils apply_affine to_matvec from_matvec append_diag
###   dot_reduce voxel_sizes obliquity rescale_affine

### ** Examples

aff <- diag(c(2, 3, 4, 1))
aff[1:3, 4] <- c(10, 20, 30)

pts <- rbind(c(1, 2, 3), c(4, 5, 6))
apply_affine(aff, pts)

mv <- to_matvec(aff)
from_matvec(mv$matrix, mv$vector)



cleanEx()
nameEx("as-ClusteredNeuroVol-DenseNeuroVol")
### * as-ClusteredNeuroVol-DenseNeuroVol

flush(stderr()); flush(stdout())

### Name: as-ClusteredNeuroVol-DenseNeuroVol
### Title: Convert ClusteredNeuroVol to DenseNeuroVol
### Aliases: as-ClusteredNeuroVol-DenseNeuroVol
###   coerce,ClusteredNeuroVol,DenseNeuroVol-method

### ** Examples


# Create a clustered volume
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
clusters <- rep(1:5, length.out=sum(mask))
cvol <- ClusteredNeuroVol(mask, clusters)

# Convert to DenseNeuroVol
dvol <- as(cvol, "DenseNeuroVol")




cleanEx()
nameEx("as.dense")
### * as.dense

flush(stderr()); flush(stdout())

### Name: as.dense
### Title: Convert to dense representation
### Aliases: as.dense

### ** Examples

# Create a sparse representation
space <- NeuroSpace(c(10,10,10,4), c(1,1,1))
mask <- array(runif(10*10*10) > 0.8, c(10,10,10))  # ~20% of voxels active
data <- matrix(rnorm(sum(mask) * 4), 4, sum(mask))  # Random data for active voxels
sparse_vec <- SparseNeuroVec(data, space, mask)

# Convert to dense representation
dense_vec <- as.dense(sparse_vec)
# The dense representation has the same dimensions but stores all voxels
identical(dim(sparse_vec), dim(dense_vec))




cleanEx()
nameEx("as.mask")
### * as.mask

flush(stderr()); flush(stdout())

### Name: as.mask
### Title: Convert to a LogicalNeuroVol
### Aliases: as.mask

### ** Examples

# Create a simple 3D volume with random values
space <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
vol <- NeuroVol(array(runif(1000), c(10,10,10)), space)

# Create a mask by thresholding (values > 0.5 become TRUE)
mask1 <- as.mask(vol > 0.5)

# Create a mask by specifying indices
indices <- which(vol > 0.8)  # get indices of high values
mask2 <- as.mask(vol, indices)

# Both masks are LogicalNeuroVol objects
identical(class(mask1), class(mask2))



cleanEx()
nameEx("as.sparse")
### * as.sparse

flush(stderr()); flush(stdout())

### Name: as.sparse
### Title: Convert to from dense to sparse representation
### Aliases: as.sparse

### ** Examples

bvol <- NeuroVol(array(runif(24*24*24), c(24,24,24)), NeuroSpace(c(24,24,24), c(1,1,1)))
indmask <- sort(sample(1:(24*24*24), 100))
svol <- as.sparse(bvol, indmask)


mask <- LogicalNeuroVol(runif(length(indmask)), space=space(bvol), indices=indmask)
sum(mask) == 100



cleanEx()
nameEx("axes-methods")
### * axes-methods

flush(stderr()); flush(stdout())

### Name: axes
### Title: Extract Image Axes
### Aliases: axes axes,NeuroSpace-method

### ** Examples

x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
class(axes(x)) == "AxisSet3D"




cleanEx()
nameEx("bilateral_filter")
### * bilateral_filter

flush(stderr()); flush(stdout())

### Name: bilateral_filter
### Title: Apply a bilateral filter to a volumetric image
### Aliases: bilateral_filter

### ** Examples

brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Apply bilateral filtering to the brain volume
filtered_vol <- bilateral_filter(brain_mask, brain_mask, spatial_sigma = 2,
intensity_sigma = 25, window = 1)




cleanEx()
nameEx("bilateral_filter_4d")
### * bilateral_filter_4d

flush(stderr()); flush(stdout())

### Name: bilateral_filter_4d
### Title: Apply a 4D bilateral filter to a NeuroVec
### Aliases: bilateral_filter_4d

### ** Examples





cleanEx()
nameEx("bounds-methods")
### * bounds-methods

flush(stderr()); flush(stdout())

### Name: bounds
### Title: Extract Spatial Bounds of an Image
### Aliases: bounds bounds,NeuroSpace-method

### ** Examples

bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
b <- bounds(bspace)
nrow(b) == ndim(bspace)
ncol(b) == 2




cleanEx()
nameEx("centroid-methods")
### * centroid-methods

flush(stderr()); flush(stdout())

### Name: centroid
### Title: return the centroid of an object
### Aliases: centroid centroid,NeuroSpace-method centroid,ROICoords-method

### ** Examples


bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
centroid(bspace)




cleanEx()
nameEx("cgb_filter")
### * cgb_filter

flush(stderr()); flush(stdout())

### Name: cgb_filter
### Title: Correlation-guided bilateral filtering (convenience wrapper)
### Aliases: cgb_filter

### ** Examples





cleanEx()
nameEx("cgb_make_graph")
### * cgb_make_graph

flush(stderr()); flush(stdout())

### Name: cgb_make_graph
### Title: Build a correlation-guided bilateral (CGB) graph
### Aliases: cgb_make_graph

### ** Examples




cleanEx()
nameEx("close-methods")
### * close-methods

flush(stderr()); flush(stdout())

### Name: close,BinaryReader-method
### Title: Close a BinaryReader or BinaryWriter
### Aliases: close,BinaryReader-method close,BinaryWriter-method

### ** Examples





cleanEx()
nameEx("cluster_searchlight_series")
### * cluster_searchlight_series

flush(stderr()); flush(stdout())

### Name: cluster_searchlight_series
### Title: Cluster-centroid searchlight over cluster time-series
### Aliases: cluster_searchlight_series

### ** Examples

# Create synthetic 4D data (8x8x8 volume, 10 timepoints)
sp4 <- NeuroSpace(c(8,8,8,10), c(1,1,1))
arr <- array(rnorm(8*8*8*10), dim=c(8,8,8,10))
vec <- NeuroVec(arr, sp4)

# Create a mask covering most of the volume
sp3 <- NeuroSpace(c(8,8,8), c(1,1,1))
mask_arr <- array(FALSE, dim=c(8,8,8))
mask_arr[2:7, 2:7, 2:7] <- TRUE
mask <- LogicalNeuroVol(mask_arr, sp3)

# Assign voxels to 10 clusters
n_voxels <- sum(mask_arr)
clusters <- sample(1:10, n_voxels, replace=TRUE)
cvol <- ClusteredNeuroVol(mask, clusters)

# Create clustered representation
cv <- ClusteredNeuroVec(vec, cvol)

# Get cluster searchlight with 3 nearest neighbors
windows <- cluster_searchlight_series(cv, k = 3)
length(windows)  # 10 windows (one per cluster)

# Check first window
roi1 <- windows[[1]]
dim(values(roi1))  # 10 x 3 (timepoints x neighbors)

# Use radius-based neighborhoods (5mm radius)
windows_radius <- cluster_searchlight_series(cv, radius = 5)
# Each window may have different number of neighbors



cleanEx()
nameEx("clustered_searchlight")
### * clustered_searchlight

flush(stderr()); flush(stdout())

### Name: clustered_searchlight
### Title: Create a clustered searchlight iterator
### Aliases: clustered_searchlight

### ** Examples

# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate a clustered searchlight iterator with 5 clusters
clust_searchlight <- clustered_searchlight(mask, csize = 5)





cleanEx()
nameEx("concat-methods")
### * concat-methods

flush(stderr()); flush(stdout())

### Name: concat
### Title: Concatenate two objects in the time dimension
### Aliases: concat concat,NeuroVec,NeuroVol-method
###   concat,NeuroVol,NeuroVec-method concat,NeuroVec,NeuroVec-method
###   concat,ROIVec,ROIVec-method concat,DenseNeuroVol,missing-method
###   concat,DenseNeuroVol,DenseNeuroVol-method
###   concat,AbstractSparseNeuroVec,missing-method
###   concat,SparseNeuroVec,SparseNeuroVec-method

### ** Examples

bv1 <- NeuroVol(rep(1,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
bv2 <- NeuroVol(rep(2,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
bv3 <- concat(bv1,bv2)
inherits(bv3, "NeuroVec")

bv4 <- concat(bv3, bv1)
dim(bv4)[4] == 3
bv5 <- concat(bv1, bv3)
dim(bv4)[4] == 3

bv6 <- concat(bv4,bv5)
dim(bv6)[4] == 6




cleanEx()
nameEx("conn_comp-methods")
### * conn_comp-methods

flush(stderr()); flush(stdout())

### Name: conn_comp
### Title: Connected components
### Aliases: conn_comp conn_comp,NeuroVol-method

### ** Examples

# Create a simple 3D volume with two distinct regions
space <- NeuroSpace(c(10,10,10), c(1,1,1))
vol_data <- array(0, c(10,10,10))

# Create first cluster in corner (2x2x2)
vol_data[1:2, 1:2, 1:2] <- 1

# Create second cluster in opposite corner (2x2x2)
vol_data[8:9, 8:9, 8:9] <- 1

# Create NeuroVol object
vol <- NeuroVol(vol_data, space)

# Find connected components with default 26-connectivity
# Returns components above threshold 0
comps <- conn_comp(vol, threshold=0)

# Access results
max(comps$index) == 2  # Should have 2 clusters
all(comps$size >= 0)    # All clusters should have >= 0

# Get cluster statistics
comps <- conn_comp(vol, threshold=0, cluster_table=TRUE)
# cluster_table contains: index, x, y, z, N (size), Area, value

# Find local maxima within clusters
comps <- conn_comp(vol, threshold=0, local_maxima=TRUE,
                  local_maxima_dist=2)
# local_maxima contains: index, x, y, z, value




cleanEx()
nameEx("conn_comp_3D")
### * conn_comp_3D

flush(stderr()); flush(stdout())

### Name: conn_comp_3D
### Title: Extract Connected Components from a 3D Binary Mask
### Aliases: conn_comp_3D

### ** Examples

# Create a simple 3D binary mask with two disconnected components
mask <- array(FALSE, c(4, 4, 4))
mask[1:2, 1:2, 1:2] <- TRUE  # First component
mask[3:4, 3:4, 3:4] <- TRUE  # Second component

# Extract components using different connectivity patterns
comps <- conn_comp_3D(mask, connect = "6-connect")

# Number of components
max_comps <- max(comps$index)
cat("Found", max_comps, "components\n")

# Size of each component
unique_sizes <- unique(comps$size[comps$size > 0])
cat("Component sizes:", paste(unique_sizes, collapse=", "), "\n")

# Try with different connectivity
comps_26 <- conn_comp_3D(mask, connect = "26-connect")
cat("Number of components with 26-connectivity:", max(comps_26$index), "\n")




cleanEx()
nameEx("coord_to_grid-methods")
### * coord_to_grid-methods

flush(stderr()); flush(stdout())

### Name: coord_to_grid
### Title: convert n-dimensional real world coordinates to grid coordinates
### Aliases: coord_to_grid coord_to_grid,NeuroSpace,matrix-method
###   coord_to_grid,NeuroSpace,numeric-method
###   coord_to_grid,NeuroVol,matrix-method
###   coord_to_grid,NeuroVol,numeric-method

### ** Examples

# Create a simple 3D volume
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
grid <- coord_to_grid(bvol, coords)
world <- grid_to_coord(bvol, grid)
all.equal(coords, world)



cleanEx()
nameEx("coord_to_index-methods")
### * coord_to_index-methods

flush(stderr()); flush(stdout())

### Name: coord_to_index
### Title: convert n-dimensional real world coordinates to 1D indices
### Aliases: coord_to_index coord_to_index,NeuroSpace,matrix-method
###   coord_to_index,NeuroSpace,numeric-method
###   coord_to_index,NeuroVol,matrix-method

### ** Examples

bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
idx <- coord_to_index(bvol, coords)
coords2 <- index_to_coord(bvol, idx)
all.equal(coords, coords2)



cleanEx()
nameEx("coords-methods")
### * coords-methods

flush(stderr()); flush(stdout())

### Name: coords,IndexLookupVol-method
### Title: Extract Coordinates from an IndexLookupVol Object
### Aliases: coords,IndexLookupVol-method coords,ROIVol-method
###   coords,ROICoords-method coords,AbstractSparseNeuroVec-method

### ** Examples


space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
coords(ilv, 1)  # Extract coordinates for index 1





cleanEx()
nameEx("coords")
### * coords

flush(stderr()); flush(stdout())

### Name: coords
### Title: Extract coordinates from an object
### Aliases: coords

### ** Examples

# Create a NeuroSpace object with 3mm voxels
space <- NeuroSpace(c(10,10,10), spacing=c(3,3,3))

# Create ROI coordinates in voxel space
coords <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)
roi_coords <- ROICoords(coords)

# Get coordinates in voxel space
vox_coords <- coords(roi_coords)
# First coordinate is (1,1,1)

# Get coordinates
cds <- coords(roi_coords)
nrow(cds) == 2



cleanEx()
nameEx("cuboid_roi")
### * cuboid_roi

flush(stderr()); flush(stdout())

### Name: cuboid_roi
### Title: Create A Cuboid Region of Interest
### Aliases: cuboid_roi

### ** Examples

 sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
 cube <- cuboid_roi(sp1, c(5,5,5), 3)
 vox <- coords(cube)
 cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=5)





cleanEx()
nameEx("data_file-methods")
### * data_file-methods

flush(stderr()); flush(stdout())

### Name: data_file
### Title: Generic function to get the name of the data file, given a file
###   name and a 'FileFormat' instance.
### Aliases: data_file data_file,FileFormat,character-method

### ** Examples


fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
data_file(fmt, "brain_scan.img")  # Returns "brain_scan.img"
data_file(fmt, "brain_scan.hdr")  # Also Returns "brain_scan.img"





cleanEx()
nameEx("data_file_matches-methods")
### * data_file_matches-methods

flush(stderr()); flush(stdout())

### Name: data_file_matches
### Title: Generic function to test whether a file name conforms to the
###   given a 'FileFormat' instance. Will test for match to data file only
### Aliases: data_file_matches
###   data_file_matches,FileFormat,character-method

### ** Examples


fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
data_file_matches(fmt, "brain_scan.img")  # TRUE
data_file_matches(fmt, "brain_scan.hdr")  # FALSE
data_file_matches(fmt, "brain.img.gz")    # FALSE





cleanEx()
nameEx("data_reader")
### * data_reader

flush(stderr()); flush(stdout())

### Name: data_reader
### Title: Create a Data Reader
### Aliases: data_reader

### ** Examples


# Create reader for NIFTI file
meta <- read_header(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
reader <- data_reader(meta, offset = 0)

# Read first 100 voxels
data <- read_elements(reader, 100)





cleanEx()
nameEx("deoblique")
### * deoblique

flush(stderr()); flush(stdout())

### Name: deoblique
### Title: Deoblique a Neuroimaging Space or Volume
### Aliases: deoblique

### ** Examples

sp <- NeuroSpace(c(32, 32, 20), spacing = c(2, 2, 3))
tx <- trans(sp)
tx[1, 2] <- 0.15
sp_obl <- NeuroSpace(dim(sp), spacing = spacing(sp), trans = tx)

# Build deobliqued target space (minimum spacing default)
sp_deob <- deoblique(sp_obl)

# Resample a volume to deobliqued space
vol <- NeuroVol(array(rnorm(prod(dim(sp_obl))), dim(sp_obl)), sp_obl)




cleanEx()
nameEx("dim_of-methods")
### * dim_of-methods

flush(stderr()); flush(stdout())

### Name: dim_of
### Title: Get the length of a given dimension of an object
### Aliases: dim_of dim_of,NeuroSpace,NamedAxis-method

### ** Examples


x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
stopifnot(dim_of(x, x@axes@i) == 10)



cleanEx()
nameEx("downsample-methods")
### * downsample-methods

flush(stderr()); flush(stdout())

### Name: downsample
### Title: Downsample an Image
### Aliases: downsample downsample,DenseNeuroVec-method
###   downsample,SparseNeuroVec-method downsample,NeuroVec-method
###   downsample,DenseNeuroVol-method downsample,NeuroVol-method

### ** Examples

# Create a sample 4D image
data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
space <- NeuroSpace(dim = c(64, 64, 32, 10), 
                    origin = c(0, 0, 0),
                    spacing = c(2, 2, 2))
nvec <- DenseNeuroVec(data, space)

# Downsample by factor
nvec_down1 <- downsample(nvec, factor = 0.5)

# Downsample to target spacing
nvec_down2 <- downsample(nvec, spacing = c(4, 4, 4))

# Downsample to target dimensions
nvec_down3 <- downsample(nvec, outdim = c(32, 32, 16))

# Create a sample 3D volume
data <- array(rnorm(64*64*32), dim = c(64, 64, 32))
space <- NeuroSpace(dim = c(64, 64, 32), 
                    origin = c(0, 0, 0),
                    spacing = c(2, 2, 2))
vol <- DenseNeuroVol(data, space)

# Downsample by factor
vol_down1 <- downsample(vol, factor = 0.5)

# Downsample to target spacing
vol_down2 <- downsample(vol, spacing = c(4, 4, 4))

# Downsample to target dimensions
vol_down3 <- downsample(vol, outdim = c(32, 32, 16))




cleanEx()
nameEx("drop_dim-methods")
### * drop_dim-methods

flush(stderr()); flush(stdout())

### Name: drop_dim
### Title: Drop a Dimension from an Object
### Aliases: drop_dim drop_dim,AxisSet2D,numeric-method
###   drop_dim,AxisSet2D,missing-method drop_dim,AxisSet3D,numeric-method
###   drop_dim,AxisSet3D,missing-method drop_dim,NeuroSpace,numeric-method
###   drop_dim,NeuroSpace,missing-method

### ** Examples

# Create a NeuroSpace object with dimensions (10, 10, 10)
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Drop the first dimension
x1 <- drop_dim(x, 1)

# Check the new dimensions
ndim(x1) == 2
dim(x1)[1] == 10




cleanEx()
nameEx("ecode_name")
### * ecode_name

flush(stderr()); flush(stdout())

### Name: ecode_name
### Title: Get Extension Code Name
### Aliases: ecode_name

### ** Examples

ecode_name(4L)   # Returns "AFNI"
ecode_name(6L)   # Returns "comment"
ecode_name(999L) # Returns "unknown"




cleanEx()
nameEx("embed_kernel-methods")
### * embed_kernel-methods

flush(stderr()); flush(stdout())

### Name: embed_kernel
### Title: Generic function to position kernel in a position in image space
### Aliases: embed_kernel embed_kernel,Kernel,NeuroSpace,numeric-method

### ** Examples

# Create a 3D Gaussian kernel with dimensions 3x3x3 and voxel size 1x1x1
kern <- Kernel(kerndim = c(3,3,3), vdim = c(1,1,1), FUN = dnorm, sd = 1)

# Create a NeuroSpace object to embed the kernel in
space <- NeuroSpace(c(10,10,10), c(1,1,1))

# Embed the kernel at the center of the space (position 5,5,5)
embedded_kern <- embed_kernel(kern, space, c(5,5,5))

# The result is a SparseNeuroVol with kernel weights centered at (5,5,5)
# We can also scale the kernel weights by using the weight parameter
embedded_kern_scaled <- embed_kernel(kern, space, c(5,5,5), weight = 2)

# The scaled kernel has weights twice as large as the original
max(values(embedded_kern_scaled)) == 2 * max(values(embedded_kern))




cleanEx()
nameEx("file_matches-methods")
### * file_matches-methods

flush(stderr()); flush(stdout())

### Name: file_matches
### Title: Generic function to test whether a file name conforms to the
###   given 'FileFormat' instance. Will test for match to either header
###   file or data file
### Aliases: file_matches file_matches,FileFormat,character-method

### ** Examples


# Create a FileFormat for NIFTI format




cleanEx()
nameEx("findAnatomy3D")
### * findAnatomy3D

flush(stderr()); flush(stdout())

### Name: findAnatomy3D
### Title: Find 3D anatomical orientation from axis abbreviations
### Aliases: findAnatomy3D

### ** Examples

# Create orientation with default LPI axes
orient <- findAnatomy3D()
# Create orientation with custom axes
orient <- findAnatomy3D("R", "A", "S")



cleanEx()
nameEx("gaussian_blur")
### * gaussian_blur

flush(stderr()); flush(stdout())

### Name: gaussian_blur
### Title: Gaussian Blur for Volumetric Images
### Aliases: gaussian_blur

### ** Examples

# Load a sample brain mask
brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Apply Gaussian blurring to the brain volume
blurred_vol <- gaussian_blur(brain_mask, brain_mask, sigma = 2, window = 1)

# View a slice of the original and blurred volumes
image(brain_mask[,,12])
image(blurred_vol[,,12])




cleanEx()
nameEx("get_afni_attribute")
### * get_afni_attribute

flush(stderr()); flush(stdout())

### Name: get_afni_attribute
### Title: Get AFNI Attribute from Extension
### Aliases: get_afni_attribute

### ** Examples

## Not run: 
##D # Get the history note from an AFNI extension
##D history <- get_afni_attribute(afni_ext, "HISTORY_NOTE")
## End(Not run)




cleanEx()
nameEx("grid_to_coord-methods")
### * grid_to_coord-methods

flush(stderr()); flush(stdout())

### Name: grid_to_coord
### Title: Generic function to convert N-dimensional grid coordinates to
###   real world coordinates
### Aliases: grid_to_coord grid_to_coord,NeuroSpace,matrix-method
###   grid_to_coord,NeuroSpace,numeric-method
###   grid_to_coord,NeuroVol,matrix-method

### ** Examples

# Create a simple 3D volume
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
grid_coords <- matrix(c(1.5,1.5,1.5, 5.5,5.5,5.5), ncol=3, byrow=TRUE)
world <- grid_to_coord(bvol, grid_coords)
grid <- coord_to_grid(bvol, world)
all.equal(grid_coords, grid)



cleanEx()
nameEx("grid_to_grid-methods")
### * grid_to_grid-methods

flush(stderr()); flush(stdout())

### Name: grid_to_grid
### Title: Generic function to convert voxel coordinates in the reference
###   space (LPI) to native array space.
### Aliases: grid_to_grid grid_to_grid,NeuroSpace,matrix-method
###   grid_to_grid,matrix,matrix-method

### ** Examples

# Create a simple 3D volume in LPI orientation
space <- NeuroSpace(c(10,10,10), c(2,2,2))

# Create a reoriented space in RAS orientation
space_ras <- reorient(space, c("R", "A", "S"))

# Convert coordinates between orientations
voxel_coords <- t(matrix(c(1,1,1)))
new_coords <- grid_to_grid(space_ras, voxel_coords)
print(new_coords)



cleanEx()
nameEx("grid_to_index-methods")
### * grid_to_index-methods

flush(stderr()); flush(stdout())

### Name: grid_to_index
### Title: Generic function to convert N-dimensional grid coordinates to 1D
###   indices
### Aliases: grid_to_index grid_to_index,NeuroSlice,matrix-method
###   grid_to_index,NeuroSlice,numeric-method
###   grid_to_index,NeuroSpace,matrix-method
###   grid_to_index,NeuroSpace,numeric-method
###   grid_to_index,NeuroVol,matrix-method
###   grid_to_index,NeuroVol,numeric-method

### ** Examples

# Create a 2D space (10x10)
space_2d <- NeuroSpace(c(10,10), c(1,1))

# Convert 2D grid coordinates to linear indices
coords_2d <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
idx_2d <- grid_to_index(space_2d, coords_2d)
# First coordinate (1,1) maps to index 1
# Second coordinate (2,2) maps to index 12 (= 2 + (2-1)*10)

# Create a 3D space (10x10x10)
space_3d <- NeuroSpace(c(10,10,10), c(1,1,1))

# Convert 3D grid coordinates to linear indices
coords_3d <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)
idx_3d <- grid_to_index(space_3d, coords_3d)

# Single coordinate can also be converted
idx <- grid_to_index(space_3d, c(1,1,1))

slice_space <- NeuroSpace(c(10, 10))
slice_data <- matrix(1:100, 10, 10)
slice <- NeuroSlice(slice_data, slice_space)

# Convert single coordinate
idx <- grid_to_index(slice, c(5, 5))

# Convert multiple coordinates
coords <- matrix(c(1,1, 2,2, 3,3), ncol=2, byrow=TRUE)
indices <- grid_to_index(slice, coords)




cleanEx()
nameEx("guided_filter")
### * guided_filter

flush(stderr()); flush(stdout())

### Name: guided_filter
### Title: Edge-Preserving Guided Filter for Volumetric Images
### Aliases: guided_filter

### ** Examples

# Load an example brain volume
brain_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Apply guided filtering to the brain volume




cleanEx()
nameEx("header_file-methods")
### * header_file-methods

flush(stderr()); flush(stdout())

### Name: header_file
### Title: Generic function to get the name of the header file, given a
###   file name and a 'FileFormat' instance.
### Aliases: header_file header_file,FileFormat,character-method

### ** Examples


fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file(fmt, "brain_scan.hdr")  # Returns "brain_scan.hdr"
header_file(fmt, "brain_scan.img")  # Returns "brain_scan.hdr"





cleanEx()
nameEx("header_file_matches-methods")
### * header_file_matches-methods

flush(stderr()); flush(stdout())

### Name: header_file_matches
### Title: Generic function to test whether a file name conforms to the
###   given 'FileFormat' instance. Will test for match to header file only
### Aliases: header_file_matches
###   header_file_matches,FileFormat,character-method

### ** Examples


fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file_matches(fmt, "brain_scan.hdr")  # TRUE
header_file_matches(fmt, "brain_scan.img")  # FALSE
header_file_matches(fmt, "brain.hdr.gz")    # FALSE





cleanEx()
nameEx("index_to_coord-methods")
### * index_to_coord-methods

flush(stderr()); flush(stdout())

### Name: index_to_coord
### Title: convert 1d indices to n-dimensional real world coordinates
### Aliases: index_to_coord index_to_coord,NeuroSpace,numeric-method
###   index_to_coord,NeuroSpace,integer-method
###   index_to_coord,NeuroVol,integer-method
###   index_to_coord,NeuroVec,integer-method

### ** Examples

bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
idx <- 1:10
g <- index_to_coord(bvol, idx)
idx2 <- coord_to_index(bvol, g)
all.equal(idx, idx2)



cleanEx()
nameEx("index_to_grid-methods")
### * index_to_grid-methods

flush(stderr()); flush(stdout())

### Name: index_to_grid
### Title: Convert 1d indices to n-dimensional grid coordinates
### Aliases: index_to_grid index_to_grid,NeuroSlice,numeric-method
###   index_to_grid,NeuroSpace,numeric-method
###   index_to_grid,NeuroVec,index-method
###   index_to_grid,NeuroVec,integer-method
###   index_to_grid,NeuroVol,index-method
###   index_to_grid,NeuroVol,integer-method

### ** Examples


 bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
 idx <- 1:10
 g <- index_to_grid(bvol, idx)
 bvol[g]

slice_space <- NeuroSpace(c(10, 10))
slice_data <- matrix(1:100, 10, 10)
slice <- NeuroSlice(slice_data, slice_space)

# Convert single index
coords <- index_to_grid(slice, 55)

# Convert multiple indices
indices <- c(1, 25, 50, 75, 100)
coords_mat <- index_to_grid(slice, indices)




cleanEx()
nameEx("indices-methods")
### * indices-methods

flush(stderr()); flush(stdout())

### Name: indices,IndexLookupVol-method
### Title: Get Indices from an IndexLookupVol Object
### Aliases: indices,IndexLookupVol-method indices,ROIVol-method
###   indices,ROIVec-method indices,AbstractSparseNeuroVec-method

### ** Examples


space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
idx <- indices(ilv)  # Get included indices





cleanEx()
nameEx("indices")
### * indices

flush(stderr()); flush(stdout())

### Name: indices
### Title: Extract indices
### Aliases: indices

### ** Examples

# Create a NeuroSpace object with 3mm voxels
space <- NeuroSpace(c(10,10,10), spacing=c(3,3,3))

# Create ROI coordinates in voxel space
coords <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)

# Create ROI volume
roi_vol <- ROIVol(space, coords, data=c(1,2))

# Get linear indices of ROI voxels
idx <- indices(roi_vol)
# These indices can be used to index into a 3D array of size 10x10x10



cleanEx()
nameEx("inverse_trans-methods")
### * inverse_trans-methods

flush(stderr()); flush(stdout())

### Name: inverse_trans
### Title: Extract inverse image coordinate transformation
### Aliases: inverse_trans inverse_trans,NeuroSpace-method

### ** Examples

bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
itrans <- inverse_trans(bspace)
identical(trans(bspace) %*% inverse_trans(bspace), diag(4))



cleanEx()
nameEx("linear_access-methods")
### * linear_access-methods

flush(stderr()); flush(stdout())

### Name: linear_access,DenseNeuroVol,numeric-method
### Title: Linear Access Method for FileBackedNeuroVec
### Aliases: linear_access,DenseNeuroVol,numeric-method
###   linear_access,DenseNeuroVec,numeric-method
###   linear_access,DenseNeuroVol,integer-method
###   linear_access,DenseNeuroVec,integer-method
###   linear_access,FileBackedNeuroVec,numeric-method
###   linear_access,MappedNeuroVec,numeric-method
###   linear_access,NeuroHyperVec,ANY-method
###   linear_access,NeuroVecSeq,numeric-method
###   linear_access,SparseNeuroVol,numeric-method
###   linear_access,AbstractSparseNeuroVec,numeric-method

### ** Examples





cleanEx()
nameEx("linear_access")
### * linear_access

flush(stderr()); flush(stdout())

### Name: linear_access
### Title: Extract values from an array-like object using linear indexing.
### Aliases: linear_access

### ** Examples

# Create a sparse neuroimaging vector
bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)

# Extract values using linear indices
# Get values from first timepoint at voxels 1,2,3
indices <- c(1,2,3)
vals <- linear_access(svec, indices)

# Get values from multiple timepoints and voxels
# First voxel at timepoint 1, second voxel at timepoint 2
indices <- c(1, 1000 + 2) # 1000 = prod(10,10,10)
vals <- linear_access(svec, indices)



cleanEx()
nameEx("list_afni_attributes")
### * list_afni_attributes

flush(stderr()); flush(stdout())

### Name: list_afni_attributes
### Title: List AFNI Attributes in Extension
### Aliases: list_afni_attributes

### ** Examples

## Not run: 
##D # List all attributes in an AFNI extension
##D attrs <- list_afni_attributes(afni_ext)
##D print(attrs)
## End(Not run)




cleanEx()
nameEx("load_data")
### * load_data

flush(stderr()); flush(stdout())

### Name: load_data
### Title: Read data from a data source.
### Aliases: load_data
### Keywords: internal

### ** Examples

# Create a NeuroVolSource from a NIFTI file and load it
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
src <- NeuroVolSource(fname)
vol <- load_data(src)
# The loaded volume is a DenseNeuroVol object
class(vol)
dim(vol)



cleanEx()
nameEx("lookup-methods")
### * lookup-methods

flush(stderr()); flush(stdout())

### Name: lookup,IndexLookupVol,numeric-method
### Title: Lookup Values in an IndexLookupVol Object
### Aliases: lookup,IndexLookupVol,numeric-method
###   lookup,AbstractSparseNeuroVec,numeric-method

### ** Examples


space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
lookup(ilv, c(1, 2, 3))  # Look up values for indices 1, 2, and 3





cleanEx()
nameEx("lookup")
### * lookup

flush(stderr()); flush(stdout())

### Name: lookup
### Title: Index Lookup operation
### Aliases: lookup

### ** Examples

# Create a 64x64x64 space
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))

# Create a lookup volume with first 100 indices
ilv <- IndexLookupVol(space, 1:100)

# Look up values for indices 1, 2, and 3
# Returns their positions in the sparse representation
lookup(ilv, c(1, 2, 3))

# Look up values outside the included indices
# Returns 0 for indices not in the lookup volume
lookup(ilv, c(101, 102))



cleanEx()
nameEx("map-methods")
### * map-methods

flush(stderr()); flush(stdout())

### Name: mapf
### Title: Apply a function to an object.
### Aliases: mapf mapf,NeuroVol,Kernel-method

### ** Examples

# Create a simple 3D volume
bspace <- NeuroSpace(c(10,10,10), c(1,1,1))
vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), bspace)

# Create a 3x3x3 mean smoothing kernel
kern <- Kernel(c(3,3,3),  vdim=c(3,3,3))

# Apply the kernel to smooth the volume
smoothed_vol <- mapf(vol, kern)



cleanEx()
nameEx("map_values-methods")
### * map_values-methods

flush(stderr()); flush(stdout())

### Name: map_values
### Title: Map Values from One Set to Another Using a User-supplied Lookup
###   Table
### Aliases: map_values map_values,NeuroVol,list-method
###   map_values,NeuroVol,matrix-method

### ** Examples

x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
vol <- NeuroVol(sample(1:10, 10 * 10 * 10, replace = TRUE), x)

## Lookup table is a list
lookup <- lapply(1:10, function(i) i * 10)
names(lookup) <- 1:10
ovol <- map_values(vol, lookup)

## Lookup table is a matrix. The first column is the key, and the second column is the value
names(lookup) <- 1:length(lookup)
lookup.mat <- cbind(as.numeric(names(lookup)), unlist(lookup))
ovol2 <- map_values(vol, lookup.mat)
all.equal(as.vector(ovol2), as.vector(ovol))




cleanEx()
nameEx("mask-methods")
### * mask-methods

flush(stderr()); flush(stdout())

### Name: mask
### Title: Extract Mask from Neuroimaging Object
### Aliases: mask mask,ClusteredNeuroVol-method
###   mask,FileBackedNeuroVec-method mask,MappedNeuroVec-method
###   mask,NeuroHyperVec-method mask,NeuroSlice-method
###   mask,DenseNeuroVec-method mask,DenseNeuroVol-method
###   mask,LogicalNeuroVol-method mask,AbstractSparseNeuroVec-method
###   mask,SparseNeuroVecSource-method

### ** Examples

# Create a dense volume
vol <- NeuroVol(array(rnorm(64^3), c(64,64,64)), NeuroSpace(c(64,64,64)))
m <- mask(vol)  # Returns all TRUE mask

# Create a sparse vector with explicit mask
mask_array <- array(runif(64^3) > 0.5, c(64,64,64))
mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(c(64,64,64)))
# Data must be a matrix (time x masked voxels)
sparse_data <- matrix(rnorm(sum(mask_array) * 10), nrow = 10, ncol = sum(mask_array))
svec <- SparseNeuroVec(sparse_data, NeuroSpace(c(64,64,64,10)), mask_vol)
m2 <- mask(svec)  # Returns the stored mask




cleanEx()
nameEx("matricized_access-methods")
### * matricized_access-methods

flush(stderr()); flush(stdout())

### Name: matricized_access
### Title: Extract values from a 4D tensor using a matrix of time-space
###   indices.
### Aliases: matricized_access
###   matricized_access,SparseNeuroVec,matrix-method
###   matricized_access,SparseNeuroVec,integer-method
###   matricized_access,SparseNeuroVec,numeric-method
###   matricized_access,BigNeuroVec,matrix-method
###   matricized_access,BigNeuroVec,integer-method
###   matricized_access,BigNeuroVec,numeric-method

### ** Examples

# Create a sparse 4D neuroimaging vector
bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)

# Extract specific timepoint-voxel pairs
# Get value at timepoint 1, voxel 1 and timepoint 2, voxel 2
idx_mat <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
vals <- matricized_access(svec, idx_mat)

# Get full time series for voxels 1 and 2
ts_mat <- matricized_access(svec, c(1,2))
# Each column in ts_mat contains the full time series for that voxel



cleanEx()
nameEx("meta_info")
### * meta_info

flush(stderr()); flush(stdout())

### Name: meta_info
### Title: Lightweight metadata for neuroimaging files
### Aliases: meta_info meta_info,FileMetaInfo-method
###   meta_info,character-method

### ** Examples





cleanEx()
nameEx("ndim-methods")
### * ndim-methods

flush(stderr()); flush(stdout())

### Name: ndim
### Title: Extract the number of dimensions of an object
### Aliases: ndim ndim,AxisSet-method ndim,ClusteredNeuroVec-method
###   ndim,NeuroObj-method ndim,NeuroHyperVec-method ndim,NeuroSpace-method

### ** Examples


x = NeuroSpace(c(10,10,10), spacing=c(1,1,1))
ndim(x) == 3
x = NeuroSpace(c(10,10,10,3), spacing=c(1,1,1))
ndim(x) == 4




cleanEx()
nameEx("not-methods")
### * not-methods

flush(stderr()); flush(stdout())

### Name: not-methods
### Title: Logical Negation for Neuroimaging Volumes
### Aliases: not-methods !,DenseNeuroVol-method !,SparseNeuroVol-method

### ** Examples

sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
mask <- LogicalNeuroVol(array(sample(c(TRUE, FALSE), 125, replace = TRUE),
                              c(5, 5, 5)), sp)
inv <- !mask




cleanEx()
nameEx("num_clusters-methods")
### * num_clusters-methods

flush(stderr()); flush(stdout())

### Name: num_clusters
### Title: Number of Clusters
### Aliases: num_clusters num_clusters,ClusteredNeuroVec-method
###   num_clusters,ClusteredNeuroVol-method

### ** Examples

# Create a simple 3D volume and mask
space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
mask_vol <- LogicalNeuroVol(vol_data > 0, space)

# Get coordinates of masked voxels for clustering
mask_idx <- which(mask_vol)
coords <- index_to_coord(mask_vol, mask_idx)

# Cluster the coordinates into 10 groups using k-means
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 10)

# Create a clustered volume
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Get the number of clusters
n_clusters <- num_clusters(clustered_vol)
n_clusters == 10



cleanEx()
nameEx("orientation_utils")
### * orientation_utils

flush(stderr()); flush(stdout())

### Name: orientation_utils
### Title: Orientation utility functions
### Aliases: orientation_utils affine_to_orientation orientation_transform
###   apply_orientation orientation_inverse_affine orientation_to_axcodes
###   axcodes_to_orientation affine_to_axcodes

### ** Examples

aff <- diag(4)
ornt <- affine_to_orientation(aff)
orientation_to_axcodes(ornt)

arr <- array(1:24, dim = c(2, 3, 4))
tx <- orientation_transform(
  axcodes_to_orientation(c("R", "A", "S")),
  axcodes_to_orientation(c("A", "R", "S"))
)
out <- apply_orientation(arr, tx)

inv_aff <- orientation_inverse_affine(tx, dim(arr))
inv_aff




cleanEx()
nameEx("origin-methods")
### * origin-methods

flush(stderr()); flush(stdout())

### Name: origin
### Title: Extract Image Origin
### Aliases: origin origin,NeuroHyperVec-method origin,NeuroSpace-method
###   origin,NeuroVol-method origin,NeuroVec-method

### ** Examples

bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
stopifnot(origin(bspace) == c(0,0,0))




cleanEx()
nameEx("parse_afni_extension")
### * parse_afni_extension

flush(stderr()); flush(stdout())

### Name: parse_afni_extension
### Title: Parse AFNI Extension
### Aliases: parse_afni_extension

### ** Examples

## Not run: 
##D # Read a NIfTI file with AFNI extension
##D hdr <- read_nifti_header("afni_file.nii")
##D afni_ext <- hdr$extensions[[1]]
##D parsed <- parse_afni_extension(afni_ext)
## End(Not run)




cleanEx()
nameEx("parse_extension")
### * parse_extension

flush(stderr()); flush(stdout())

### Name: parse_extension
### Title: Parse NIfTI Extension Data
### Aliases: parse_extension

### ** Examples

# Parse a comment extension
ext <- NiftiExtension(ecode = 6L, data = "Test comment")
parse_extension(ext)  # Returns "Test comment"




cleanEx()
nameEx("partition-methods")
### * partition-methods

flush(stderr()); flush(stdout())

### Name: partition
### Title: Partition an image into a set of disjoint clusters
### Aliases: partition partition,LogicalNeuroVol,integer-method
###   partition,LogicalNeuroVol,numeric-method
###   partition,DenseNeuroVol,numeric-method

### ** Examples

# Load an example 3D image
library(neuroim2)
img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Partition the image into 5 clusters using default options
clusters <- partition(img, 5)





cleanEx()
nameEx("patch_set")
### * patch_set

flush(stderr()); flush(stdout())

### Name: patch_set
### Title: Generate a set of coordinate "patches" of fixed size from an
###   image object.
### Aliases: patch_set

### ** Examples

# Create a simple 3D volume
space <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
vol <- NeuroVol(array(rnorm(1000), c(10,10,10)), space)

# Create a mask with some active voxels
mask <- LogicalNeuroVol(vol > 0, space)

# Extract 3x3x3 patches centered at each active voxel
patches <- patch_set(vol, dims=c(3,3,3), mask=mask)

# Access the first patch
patch1 <- patches[[1]]
dim(patch1)  # Should be c(27) (flattened 3x3x3 patch)



cleanEx()
nameEx("perm_mat-methods")
### * perm_mat-methods

flush(stderr()); flush(stdout())

### Name: perm_mat
### Title: Extract permutation matrix associated with an image
### Aliases: perm_mat perm_mat,AxisSet2D-method perm_mat,AxisSet3D-method
###   perm_mat,NeuroSpace-method

### ** Examples


fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vol <- read_vol(fname)
pmat <- perm_mat(space(vol))

vox <- c(12,12,8)
pvox <- vox %*% perm_mat(space(vol))

stopifnot(all(pvox == c(-12,12,8)))



cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

### Name: plot,NeuroSlice-method
### Title: Plot a NeuroSlice
### Aliases: plot,NeuroSlice-method plot,NeuroSlice,ANY-method
###   plot,NeuroVol-method plot,NeuroVol,missing-method
###   plot,NeuroVol,NeuroVol-method

### ** Examples

# Create example slice
slice_space <- NeuroSpace(c(100, 100))
slice_data <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(slice_data, slice_space)


dat <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(dat, NeuroSpace(c(100,100)))



cleanEx()
nameEx("random_searchlight")
### * random_searchlight

flush(stderr()); flush(stdout())

### Name: random_searchlight
### Title: Create a spherical random searchlight iterator
### Aliases: random_searchlight

### ** Examples

# Create a simple brain mask
mask_data <- array(TRUE, c(10, 10, 10))
mask_data[1, 1, 1] <- FALSE
mask <- LogicalNeuroVol(mask_data, NeuroSpace(c(10,10,10)))

# Generate random searchlight iterator with a radius of 2 voxels




cleanEx()
nameEx("read_columns-methods")
### * read_columns-methods

flush(stderr()); flush(stdout())

### Name: read_columns
### Title: Read a set of column vector from an input source (e.g.
###   'ColumnReader')
### Aliases: read_columns read_columns,ColumnReader,integer-method
### Keywords: internal

### ** Examples

# Create a reader function that returns random data
reader_func <- function(cols) {
  matrix(rnorm(100 * length(cols)), 100, length(cols))
}

# Create a ColumnReader with 100 rows and 10 columns
col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)

# Read columns 1, 3, and 5
cols <- read_columns(col_reader, c(1L, 3L, 5L))
dim(cols) == c(100, 3)




cleanEx()
nameEx("read_elements-methods")
### * read_elements-methods

flush(stderr()); flush(stdout())

### Name: read_elements
### Title: Read a sequence of elements from an input source
### Aliases: read_elements read_elements,BinaryReader,numeric-method
### Keywords: internal

### ** Examples

# Create a temporary binary file with test data
tmp <- tempfile()
con <- file(tmp, "wb")
test_data <- rnorm(100)
writeBin(test_data, con, size = 8)
close(con)

# Create a BinaryReader and read the data
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "double", bytes_per_element = 8L)
data <- read_elements(reader, 100)
close(reader)

# Clean up
unlink(tmp)



cleanEx()
nameEx("read_header")
### * read_header

flush(stderr()); flush(stdout())

### Name: read_header
### Title: read header information of an image file
### Aliases: read_header

### ** Examples

hdr <- read_header(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
dim(hdr)                  # image dimensions
hdr@header$pixdim[5]      # TR in seconds



cleanEx()
nameEx("read_image")
### * read_image

flush(stderr()); flush(stdout())

### Name: read_image
### Title: read_image
### Aliases: read_image

### ** Examples

vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))




cleanEx()
nameEx("read_meta_info-methods")
### * read_meta_info-methods

flush(stderr()); flush(stdout())

### Name: read_meta_info
### Title: Generic function to read image meta info given a file
### Aliases: read_meta_info read_meta_info,NIFTIFormat-method
###   read_meta_info,AFNIFormat-method

### ** Examples

# Create a NIFTI format descriptor
fmt <- new("NIFTIFormat",
           file_format = "NIFTI",
           header_encoding = "raw",
           header_extension = "nii",
           data_encoding = "raw",
           data_extension = "nii")

# Read metadata from a NIFTI file




cleanEx()
nameEx("read_vec")
### * read_vec

flush(stderr()); flush(stdout())

### Name: read_vec
### Title: read_vec
### Aliases: read_vec

### ** Examples


# Load a single NIfTI file
img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))


# Memory-mapped loading for large files
big_img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"), mode="mmap")

# Load masked data for memory efficiency
mask <- as.logical(big_img[[1]])
masked_data <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"),
               mask=mask, mode="bigvec")





cleanEx()
nameEx("read_vol")
### * read_vol

flush(stderr()); flush(stdout())

### Name: read_vol
### Title: Load an image volume from a file
### Aliases: read_vol

### ** Examples

fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
x <- read_vol(fname)
print(dim(x))
space(x)




cleanEx()
nameEx("reorient-methods")
### * reorient-methods

flush(stderr()); flush(stdout())

### Name: reorient
### Title: Remap the grid-to-world coordinates mapping of an image.
### Aliases: reorient reorient,NeuroSpace,character-method

### ** Examples

# Create a NeuroSpace object in LPI (Left-Posterior-Inferior) orientation
space <- NeuroSpace(c(64, 64, 40), c(2, 2, 2))

# Reorient to RAS (Right-Anterior-Superior) orientation
# Use individual axis codes: "R" for Right, "A" for Anterior, "S" for Superior
space_ras <- reorient(space, c("R", "A", "S"))

# The transformation matrix will be updated to reflect the new orientation
# Original and reoriented spaces will have different coordinate mappings
coords <- c(32, 32, 20)
orig_world <- grid_to_coord(space, coords)
new_world <- grid_to_coord(space_ras, coords)



cleanEx()
nameEx("resample-methods")
### * resample-methods

flush(stderr()); flush(stdout())

### Name: resample
### Title: Resample an Image to Match the Space of Another Image
### Aliases: resample resample,NeuroVol,NeuroVol-method
###   resample,NeuroVol,NeuroSpace-method
###   resample,ClusteredNeuroVol,NeuroSpace-method
###   resample,ClusteredNeuroVol,NeuroVol-method

### ** Examples


img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
rspace <- space(img)


newtrans4X3 <- trans(img)[1:4, 1:3]
newtrans4X3 <- newtrans4X3 * c(.5,.5,.5,1)
newtrans <- cbind(newtrans4X3, c(space(img)@origin,1))

rspace <- NeuroSpace(rspace@dim*2, rspace@spacing/2, origin=rspace@origin, trans=trans(img))


# Create source and target volumes
src_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
targ_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Resample source to match target
resampled <- resample(src_vol, targ_vol, interpolation=1)





cleanEx()
nameEx("resample_to")
### * resample_to

flush(stderr()); flush(stdout())

### Name: resample_to
### Title: Resample an image with readable method names
### Aliases: resample_to

### ** Examples




cleanEx()
nameEx("resampled_searchlight")
### * resampled_searchlight

flush(stderr()); flush(stdout())

### Name: resampled_searchlight
### Title: Create a resampled searchlight iterator
### Aliases: resampled_searchlight bootstrap_searchlight

### ** Examples

# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate a resampled searchlight iterator with radii drawn from {4,6,8}
searchlights <- resampled_searchlight(mask, radius = c(4, 6, 8))

# Use a custom shape: random ellipsoid scaled along each axis
ellipsoid_fun <- function(mask, center, radius, iter, nonzero) {
  scales <- runif(3, 0.5, 1.5)        # axis-wise stretch/compress
  vox <- spherical_roi(mask, center, radius, nonzero = FALSE)@coords
  ctr_mat <- matrix(center, nrow(vox), 3, byrow = TRUE)
  keep <- rowSums(((vox - ctr_mat) * scales)^2) <= radius^2
  vox[keep, , drop = FALSE]
}
ellip_searchlights <- resampled_searchlight(mask, radius = c(4, 6),
                                            iter = 50, shape_fun = ellipsoid_fun)

# Or use built-in named shapes
ellip_builtin <- resampled_searchlight(mask, radius = 6, shape_fun = "ellipsoid")
cube_builtin  <- resampled_searchlight(mask, radius = 6, shape_fun = "cube")





cleanEx()
nameEx("scale_series-methods")
### * scale_series-methods

flush(stderr()); flush(stdout())

### Name: scale_series
### Title: Generic functions to scale (center and/or normalize by standard
###   deviation) each series of a 4D image That is, if the 4th dimension is
###   'time' each series is a 1D time series.
### Aliases: scale_series scale_series,NeuroVec,logical,missing-method
###   scale_series,DenseNeuroVec,logical,logical-method
###   scale_series,SparseNeuroVec,logical,logical-method
###   scale_series,NeuroVec,logical,logical-method
###   scale_series,NeuroVec,missing,logical-method
###   scale_series,NeuroVec,missing,missing-method

### ** Examples

bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
res <- scale_series(bvec, TRUE, TRUE)



cleanEx()
nameEx("searchlight")
### * searchlight

flush(stderr()); flush(stdout())

### Name: searchlight
### Title: Create an exhaustive searchlight iterator
### Aliases: searchlight

### ** Examples

# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate an exhaustive searchlight iterator with a radius of 6 mm

searchlights <- searchlight(mask, radius = 6, eager = FALSE)





cleanEx()
nameEx("searchlight_coords")
### * searchlight_coords

flush(stderr()); flush(stdout())

### Name: searchlight_coords
### Title: Create an exhaustive searchlight iterator for voxel coordinates
###   using spherical_roi
### Aliases: searchlight_coords

### ** Examples

# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate an exhaustive searchlight iterator with a radius of 6 mm

searchlights <- searchlight_coords(mask, radius = 6)





cleanEx()
nameEx("searchlight_shape_functions")
### * searchlight_shape_functions

flush(stderr()); flush(stdout())

### Name: searchlight_shape_functions
### Title: Convenience shape generators for 'resampled_searchlight()'
### Aliases: searchlight_shape_functions ellipsoid_shape cube_shape
###   blobby_shape

### ** Examples

mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Ellipsoid stretched along z with modest per-iteration jitter
sl_ellip <- resampled_searchlight(mask, radius = 6,
                                   shape_fun = ellipsoid_shape(scales = c(1, 1, 1.4),
                                                              jitter = 0.1))

# Simple axis-aligned cube (Chebyshev ball)
sl_cube <- resampled_searchlight(mask, radius = 5, shape_fun = "cube")

# Blobby sphere with 40% dropout on boundary voxels
sl_blob <- resampled_searchlight(mask, radius = 6,
                                 shape_fun = blobby_shape(drop = 0.4, edge_fraction = 0.6))




cleanEx()
nameEx("series-methods")
### * series-methods

flush(stderr()); flush(stdout())

### Name: series
### Title: Extract one or more series from object
### Aliases: series series,ClusteredNeuroVec,numeric-method
###   series,NeuroHyperVec,ANY-method series,NeuroVec,matrix-method
###   series_roi,NeuroVec,matrix-method series,NeuroVec,ROICoords-method
###   series_roi,NeuroVec,ROICoords-method
###   series,NeuroVec,LogicalNeuroVol-method
###   series,NeuroVec,NeuroVol-method
###   series_roi,NeuroVec,LogicalNeuroVol-method
###   series,NeuroVec,integer-method series,DenseNeuroVec,integer-method
###   series,NeuroVec,numeric-method series_roi,NeuroVec,numeric-method
###   series,NeuroVecSeq,integer-method series,NeuroVecSeq,numeric-method
###   series,NeuroVecSeq,matrix-method series_roi,NeuroVecSeq,matrix-method
###   series,AbstractSparseNeuroVec,ROICoords-method
###   series,AbstractSparseNeuroVec,matrix-method
###   series,AbstractSparseNeuroVec,numeric-method
###   series,AbstractSparseNeuroVec,integer-method

### ** Examples

# Create a simple 4D neuroimaging vector (10x10x10 volume with 20 timepoints)
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Extract time series using linear indices
ts1 <- series(vec, 1:10)  # Get time series for first 10 voxels

# Extract time series using 3D coordinates
coords <- matrix(c(1,1,1, 2,2,2, 3,3,3), ncol=3, byrow=TRUE)
ts2 <- series(vec, coords)  # Get time series for 3 specific voxel locations

# Extract single time series using x,y,z coordinates
ts3 <- series(vec, 5, 5, 5)  # Get time series from middle voxel




cleanEx()
nameEx("series_roi")
### * series_roi

flush(stderr()); flush(stdout())

### Name: series_roi
### Title: Extract time series from specific voxel coordinates and return
###   as ROI object
### Aliases: series_roi

### ** Examples

# Create a simple 4D neuroimaging vector
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Extract time series for first 100 voxels as ROI
roi1 <- series_roi(vec, 1:100)

# Extract time series using 3D coordinates
coords <- matrix(c(1,1,1, 2,2,2, 3,3,3), ncol=3, byrow=TRUE)
roi2 <- series_roi(vec, coords)



cleanEx()
nameEx("simulate_fmri")
### * simulate_fmri

flush(stderr()); flush(stdout())

### Name: simulate_fmri
### Title: Simulate fMRI Data
### Aliases: simulate_fmri

### ** Examples

# Create a simple spherical mask
dims <- c(32, 32, 20)
mask_array <- array(FALSE, dims)
center <- dims / 2
for (i in 1:dims[1]) {
  for (j in 1:dims[2]) {
    for (k in 1:dims[3]) {
      if (sum(((c(i,j,k) - center) / (dims/3))^2) <= 1) {
        mask_array[i,j,k] <- TRUE
      }
    }
  }
}

mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))

# Simulate 100 time points
sim_data <- simulate_fmri(mask, n_time = 100, seed = 42)

# Check dimensions
dim(sim_data)  # Should be c(32, 32, 20, 100)




cleanEx()
nameEx("slices-methods")
### * slices-methods

flush(stderr()); flush(stdout())

### Name: slices
### Title: Extract an ordered series of 2D slices from a 3D or 4D object
### Aliases: slices slices,NeuroVol-method

### ** Examples

# Create a simple 3D volume
space <- NeuroSpace(c(10,10,10), c(1,1,1))
vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), space)

# Get all slices along the z-axis
slc <- slices(vol)

# Number of slices equals the z dimension
length(slc) == dim(vol)[3]

# Each slice is a 2D matrix
dim(slc[[1]]) == c(10,10)



cleanEx()
nameEx("space-methods")
### * space-methods

flush(stderr()); flush(stdout())

### Name: space
### Title: Extract Geometric Properties of an Image
### Aliases: space space,ClusteredNeuroVec-method
###   space,IndexLookupVol-method space,ROICoords-method
###   space,NeuroObj-method space,NeuroHyperVec-method
###   space,NeuroSpace-method

### ** Examples

# Create a NeuroSpace object with dimensions (10, 10, 10) and voxel size (1, 1, 1)
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Create a NeuroVol object with random data and the specified NeuroSpace
vol <- NeuroVol(rnorm(10 * 10 * 10), x)

# Retrieve the geometric properties of the NeuroVol object
identical(x, space(vol))

space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
space(ilv)  # Get the associated NeuroSpace object





cleanEx()
nameEx("space_utils")
### * space_utils

flush(stderr()); flush(stdout())

### Name: space_utils
### Title: Space utility functions
### Aliases: space_utils output_aligned_space vox2out_vox
###   slice_to_volume_affine slice2volume

### ** Examples

sp <- NeuroSpace(c(10L, 8L, 6L), spacing = c(2, 2, 2))
out <- output_aligned_space(sp)
out$shape
out$affine

slice_aff <- slice_to_volume_affine(index = 3, axis = 3, shape = c(10, 8, 6))
slice_aff



cleanEx()
nameEx("spacing-methods")
### * spacing-methods

flush(stderr()); flush(stdout())

### Name: spacing
### Title: Extract Voxel Dimensions of an Image
### Aliases: spacing spacing,ROICoords-method spacing,NeuroObj-method
###   spacing,NeuroHyperVec-method spacing,NeuroSpace-method

### ** Examples

bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
all.equal(spacing(bspace), c(2, 2, 2))




cleanEx()
nameEx("spherical_roi")
### * spherical_roi

flush(stderr()); flush(stdout())

### Name: spherical_roi
### Title: Create a Spherical Region of Interest
### Aliases: spherical_roi

### ** Examples

 sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))
 # create an ROI centered around the integer-valued positive voxel coordinate: i=5, j=5, k=5
 cube <- spherical_roi(sp1, c(5,5,5), 3.5)
 vox <- coords(cube)
 cds <- coords(cube, real=TRUE)
 ## fill in ROI with value of 6
 cube1 <- spherical_roi(sp1, c(5,5,5), 3.5, fill=6)
 all(cube1 == 6)

 ## Create multiple spherical ROIs at once (preferred):
 centers <- rbind(c(5,5,5), c(3,3,3), c(7,7,7))
 vols <- spherical_roi_set(bvol = sp1,
                          centroids = centers, radius = 3.5, fill = 1)
 length(vols)  # 3

 ## Equivalent, less efficient lapply variant:
 vols2 <- lapply(seq_len(nrow(centers)), function(i) {
   spherical_roi(sp1, centers[i,], radius = 3.5, fill = 1)
 })

 # create an ROI centered around the real-valued coordinates: x=5, y=5, z=5
 vox <- coord_to_grid(sp1, c(5, 5, 5))
 cube <- spherical_roi(sp1, vox, 3.5)



cleanEx()
nameEx("spherical_roi_set")
### * spherical_roi_set

flush(stderr()); flush(stdout())

### Name: spherical_roi_set
### Title: Create Multiple Spherical Regions of Interest
### Aliases: spherical_roi_set

### ** Examples

# Create a NeuroSpace object
sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))

# Create multiple ROIs centered at different voxel coordinates
centroids <- matrix(c(5,5,5, 3,3,3, 7,7,7), ncol=3, byrow=TRUE)
rois <- spherical_roi_set(sp1, centroids, 3.5)

# Create ROIs with specific fill values
rois <- spherical_roi_set(sp1, centroids, 3.5, fill=c(1,2,3))




cleanEx()
nameEx("split_blocks-methods")
### * split_blocks-methods

flush(stderr()); flush(stdout())

### Name: split_blocks
### Title: Cut a vector-valued object into a list of sub-blocks
### Aliases: split_blocks split_blocks,NeuroVec,integer-method
###   split_blocks,NeuroVec,factor-method

### ** Examples

# Create a 4D neuroimaging vector with 20 timepoints
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Split into 4 blocks by assigning timepoints to blocks 1-4 repeatedly
block_indices <- rep(1:4, length.out=20)
blocks <- split_blocks(vec, block_indices)




cleanEx()
nameEx("split_clusters-methods")
### * split_clusters-methods

flush(stderr()); flush(stdout())

### Name: split_clusters
### Title: Cut an object into a list of spatial or spatiotemporal clusters
### Aliases: split_clusters
###   split_clusters,NeuroVec,ClusteredNeuroVol-method
###   split_clusters,NeuroVec,integer-method
###   split_clusters,NeuroVol,ClusteredNeuroVol-method
###   split_clusters,NeuroVol,integer-method
###   split_clusters,NeuroVol,numeric-method
###   split_clusters,ClusteredNeuroVol,missing-method
###   split_clusters,NeuroVec,numeric-method

### ** Examples



# Create a simple example space and data
space <- NeuroSpace(c(10, 10, 10,4))
data <- array(rnorm(1000*4), dim = c(10, 10, 10,4))
vec <- NeuroVec(data, space)

# Create a mask for clustering (e.g., values > 0)
mask <- vec[,,,1] > 0
mask_vol <- LogicalNeuroVol(as.array(mask), NeuroSpace(c(10, 10, 10)))

# Get coordinates of masked voxels for clustering
mask_idx <- which(mask)
coords <- index_to_coord(mask_vol, mask_idx)

# Perform clustering on the coordinates (3 clusters for example)
set.seed(123) # for reproducibility
kmeans_result <- kmeans(coords, centers = 3)

# Create a ClusteredNeuroVol
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Split the NeuroVec by clusters
split_result <- split_clusters(vec, clustered_vol)

# Calculate mean value for each cluster
cluster_means <- sapply(split_result, function(x) mean(values(x)))
print(cluster_means)

# Alternative: using integer cluster assignments
cluster_indices <- numeric(prod(dim(space)[1:3]))
cluster_indices[mask_idx] <- kmeans_result$cluster
split_result2 <- split_clusters(vec, as.integer(cluster_indices))

# Verify both methods give same results
cluster_means2 <- sapply(split_result2, function(x) mean(values(x)))
print(all.equal(sort(cluster_means), sort(cluster_means2)))



cleanEx()
nameEx("split_fill-methods")
### * split_fill-methods

flush(stderr()); flush(stdout())

### Name: split_fill
### Title: Fill Disjoint Sets of Values with the Output of a Function
### Aliases: split_fill split_fill,NeuroVol,factor,function-method

### ** Examples

## Summarize with mean -- FUN returns a scalar
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
vol <- NeuroVol(rnorm(10 * 10 * 10), x)
fac <- factor(rep(1:10, length.out=1000))
ovol.mean <- split_fill(vol, fac, mean)
identical(dim(ovol.mean), dim(vol))
length(unique(as.vector(ovol.mean))) == 10

## Transform by reversing vector -- FUN returns a vector
ovol2 <- split_fill(vol, fac, rev)




cleanEx()
nameEx("split_reduce-methods")
### * split_reduce-methods

flush(stderr()); flush(stdout())

### Name: split_reduce
### Title: Summarize Subsets of an Object by Splitting by Row and Applying
###   a Summary Function
### Aliases: split_reduce split_reduce,matrix,integer,function-method
###   split_reduce,matrix,factor,missing-method
###   split_reduce,matrix,factor,function-method
###   split_reduce,NeuroVec,factor,function-method
###   split_reduce,NeuroVec,factor,missing-method

### ** Examples

mat = matrix(rnorm(100*100), 100, 100)
fac = factor(sample(1:3, nrow(mat), replace=TRUE))
## Compute column means of each sub-matrix
ms <- split_reduce(mat, fac)
all.equal(row.names(ms), levels(fac))

## Compute column medians of each sub-matrix
ms <- split_reduce(mat, fac, median)

## Compute time-series means grouped over voxels.
## Here, 'length(fac)' must equal the number of voxels: 'prod(dim(bvec)[1:3])'
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
ms <- split_reduce(bvec, fac)
ms2 <- split_reduce(bvec, fac, mean)
all.equal(row.names(ms), levels(fac))
all.equal(ms, ms2)




cleanEx()
nameEx("split_scale-methods")
### * split_scale-methods

flush(stderr()); flush(stdout())

### Name: split_scale
### Title: Center and/or Scale Row-subsets of a Matrix or Matrix-like
###   Object
### Aliases: split_scale split_scale,matrix,factor,logical,logical-method
###   split_scale,matrix,factor,missing,missing-method
###   split_scale,DenseNeuroVec,factor,missing,missing-method
###   split_scale,DenseNeuroVec,factor,logical,missing-method
###   split_scale,DenseNeuroVec,factor,logical,logical-method

### ** Examples


M <- matrix(rnorm(1000), 10, 100)
fac <- factor(rep(1:2, each=5))
Ms <- split_scale(M, fac)

## Correctly centered
all(abs(apply(Ms[fac == 1,], 2, mean)) < .000001)
all(abs(apply(Ms[fac == 2,], 2, mean)) < .000001)

## Correctly scaled
all.equal(apply(Ms[fac == 1,], 2, sd), rep(1, ncol(Ms)))
all.equal(apply(Ms[fac == 2,], 2, sd), rep(1, ncol(Ms)))



cleanEx()
nameEx("square_roi")
### * square_roi

flush(stderr()); flush(stdout())

### Name: square_roi
### Title: Create a square region of interest
### Aliases: square_roi

### ** Examples

sp1 <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
square <- square_roi(sp1, c(5, 5, 5), 1)
vox <- coords(square)
## a 3 X 3 X 1 grid
nrow(vox) == 9



cleanEx()
nameEx("strip_extension-methods")
### * strip_extension-methods

flush(stderr()); flush(stdout())

### Name: strip_extension
### Title: Generic function to strip extension from file name, given a
###   'FileFormat' instance.
### Aliases: strip_extension strip_extension,FileFormat,character-method

### ** Examples

# Create a FileFormat for NIFTI files
fmt <- new("FileFormat",
           header_extension = "nii",
           data_extension = "nii")

# Strip extension from a NIFTI file
strip_extension(fmt, "brain_scan.nii")  # Returns "brain_scan"




cleanEx()
nameEx("sub_clusters")
### * sub_clusters

flush(stderr()); flush(stdout())

### Name: sub_clusters
### Title: Select a Subset of Clusters
### Aliases: sub_clusters sub_clusters,ClusteredNeuroVec,integer-method
###   sub_clusters,ClusteredNeuroVec,numeric-method
###   sub_clusters,ClusteredNeuroVec,character-method
###   sub_clusters,ClusteredNeuroVol,integer-method
###   sub_clusters,ClusteredNeuroVol,numeric-method
###   sub_clusters,ClusteredNeuroVol,character-method

### ** Examples

sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
mask <- LogicalNeuroVol(array(c(rep(TRUE, 500), rep(FALSE, 500)),
                              c(10, 10, 10)), sp)
clusters <- rep(1:5, length.out = 500)
cvol <- ClusteredNeuroVol(mask, clusters,
          label_map = list(A = 1, B = 2, C = 3, D = 4, E = 5))
# By integer ID
sub <- sub_clusters(cvol, c(1L, 3L))
# By name
sub2 <- sub_clusters(cvol, c("A", "C"))




cleanEx()
nameEx("sub_vector-methods")
### * sub_vector-methods

flush(stderr()); flush(stdout())

### Name: sub_vector
### Title: Generic function to extract a sub-vector from a 'NeuroVec'
###   object.
### Aliases: sub_vector sub_vector,FileBackedNeuroVec,numeric-method
###   sub_vector,NeuroVec,numeric-method
###   sub_vector,NeuroVecSeq,numeric-method
###   sub_vector,AbstractSparseNeuroVec,numeric-method

### ** Examples

bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
vec <- sub_vector(bvec,1:2)
all.equal(2, dim(vec)[4])

vec <- sub_vector(bvec, c(1,3,5,7))
all.equal(4, dim(vec)[4])

mask <- LogicalNeuroVol(rep(TRUE, 24*24*24), NeuroSpace(c(24,24,24), c(1,1,1)))
svec <- SparseNeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)),
NeuroSpace(c(24,24,24,24), c(1,1,1)), mask)
vec <- sub_vector(svec, c(1,3,5))
all.equal(3, dim(vec)[4])



cleanEx()
nameEx("trans-methods")
### * trans-methods

flush(stderr()); flush(stdout())

### Name: trans
### Title: Extract image coordinate transformation
### Aliases: trans trans,MetaInfo-method trans,NIFTIMetaInfo-method
###   trans,NeuroObj-method trans,NeuroHyperVec-method
###   trans,NeuroSpace-method

### ** Examples

bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
trans(bspace)
all.equal(dim(trans(bspace)), c(4,4))



cleanEx()
nameEx("values-methods")
### * values-methods

flush(stderr()); flush(stdout())

### Name: values
### Title: Extract Data Values of an Object
### Aliases: values values,ClusteredNeuroVec-method
###   values,DenseNeuroVol-method values,SparseNeuroVol-method
###   values,ROIVol-method values,ROIVec-method

### ** Examples

x <- NeuroSpace(c(10,10,10), c(1,1,1))
vol <- NeuroVol(rnorm(10 * 10 * 10), x)
values(vol)



cleanEx()
nameEx("vec_from_vols")
### * vec_from_vols

flush(stderr()); flush(stdout())

### Name: vec_from_vols
### Title: Create NeuroVec from list of NeuroVol objects
### Aliases: vec_from_vols

### ** Examples

# Create a simple NeuroVec from list of volumes
spc <- NeuroSpace(c(10, 10, 10))
vol1 <- NeuroVol(rnorm(10*10*10), spc)
vol2 <- NeuroVol(rnorm(10*10*10), spc)
vec <- vec_from_vols(list(vol1, vol2))
print(dim(vec))  # Should be c(10, 10, 10, 2)




cleanEx()
nameEx("vectors-methods")
### * vectors-methods

flush(stderr()); flush(stdout())

### Name: vectors
### Title: Extract an ordered list of 1D vectors.
### Aliases: vectors vectors,NeuroVec,missing-method
###   vectors,DenseNeuroVec,missing-method vectors,NeuroVec,numeric-method
###   vectors,NeuroVec,logical-method vectors,NeuroVecSeq,missing-method
###   vectors,NeuroVecSeq,numeric-method vectors,NeuroVecSeq,logical-method
###   vectors,ROIVec,missing-method vectors,matrix,missing-method
###   vectors,ROIVec,integer-method vectors,matrix,integer-method
###   vectors,matrix,numeric-method vectors,ROIVec,numeric-method
###   vectors,ROIVec,logical-method vectors,SparseNeuroVec,missing-method

### ** Examples


file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vec <- read_vec(file_name)
v <- vectors(vec)
mean(v[[1]])



cleanEx()
nameEx("vols-methods")
### * vols-methods

flush(stderr()); flush(stdout())

### Name: vols
### Title: Extract an ordered series of 3D volumes.
### Aliases: vols vols,NeuroVec,numeric-method vols,NeuroVec,missing-method

### ** Examples

vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
vs <- vols(vec)
length(vs) == dim(vec)[4]

vs <- vols(vec, indices=1:3)
length(vs) == 3



cleanEx()
nameEx("voxels-methods")
### * voxels-methods

flush(stderr()); flush(stdout())

### Name: voxels
### Title: extract voxel coordinates
### Aliases: voxels voxels,Kernel-method

### ** Examples

# Create a 3D kernel with dimensions 3x3x3 and voxel size 1x1x1
kern <- Kernel(kerndim = c(3,3,3), vdim = c(1,1,1))

# Get voxel coordinates centered at origin (0,0,0)
vox <- voxels(kern)
# Returns a matrix where each row is a voxel coordinate
# relative to the kernel center

# Get voxel coordinates centered at specific point (5,5,5)
vox_centered <- voxels(kern, center_voxel = c(5,5,5))
# Returns coordinates shifted to be centered at (5,5,5)




cleanEx()
nameEx("which_dim-methods")
### * which_dim-methods

flush(stderr()); flush(stdout())

### Name: which_dim
### Title: Find Dimensions of a Given Axis
### Aliases: which_dim which_dim,NeuroSpace,NamedAxis-method

### ** Examples


x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
which_dim(x, x@axes@j) == 2



cleanEx()
nameEx("write_elements-methods")
### * write_elements-methods

flush(stderr()); flush(stdout())

### Name: write_elements
### Title: Write a sequence of elements from an input source
### Aliases: write_elements write_elements,BinaryWriter,numeric-method

### ** Examples

# Create a temporary binary file for writing
tmp <- tempfile()
writer <- BinaryWriter(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)

# Write some random data
data <- rnorm(100)
write_elements(writer, data)
close(writer)

# Read back the data to verify
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "double", bytes_per_element = 8L)
read_data <- read_elements(reader, 100)
close(reader)

# Verify data was written correctly
all.equal(data, read_data)

# Clean up
unlink(tmp)



cleanEx()
nameEx("write_vec-methods")
### * write_vec-methods

flush(stderr()); flush(stdout())

### Name: write_vec
### Title: Write a 4d image vector to disk
### Aliases: write_vec
###   write_vec,NeuroHyperVec,character,missing,missing-method
###   write_vec,NeuroHyperVec,character,character,missing-method
###   write_vec,NeuroHyperVec,character,missing,character-method
###   write_vec,NeuroHyperVec,character,missing,character,ANY-method
###   write_vec,ROIVec,character,missing,missing-method
###   write_vec,NeuroVec,character,missing,missing-method
###   write_vec,NeuroVec,character,character,missing-method
###   write_vec,NeuroVec,character,missing,character-method
###   write_vec,NeuroVec,character,missing,character,ANY-method

### ** Examples


bvec <- NeuroVec(array(0, c(10,10,10,10)), NeuroSpace(c(10,10,10,10), c(1,1,1)))



cleanEx()
nameEx("write_vol-methods")
### * write_vol-methods

flush(stderr()); flush(stdout())

### Name: write_vol
### Title: Write a 3d image volume to disk
### Aliases: write_vol write_vol,NeuroVol,character,missing,missing-method
###   write_vol,ClusteredNeuroVol,character,missing,missing-method
###   write_vol,NeuroVol,character,character,missing-method
###   write_vol,ROIVol,character,character,missing-method
###   write_vol,NeuroVol,character,missing,character-method

### ** Examples


bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

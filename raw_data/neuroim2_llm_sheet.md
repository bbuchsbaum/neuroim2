<metainfo>
Package: neuroim2
Description: R package for neuroimaging data analysis and manipulation
Core Classes: NeuroSpace, NeuroVol, NeuroVec (Dense/Sparse), NeuroHyperVec
File Formats: NIfTI (.nii, .nii.gz)
Memory Modes: normal, mmap, bigvec, filebacked
Version: 2.0.0
</metainfo>

<prompt>
You are teaching another llm how to use the neuroim2 package. Focus on the main points: reading volumes/vectors, creating NeuroSpace objects, creating NeuroVol objects, creating NeuroVec (Dense/Sparse), accessing data via [, series, also describe other functions that seem improtant due to frequency of occurence or your intuition.
</prompt>

<core_concepts>
1. NeuroSpace: Fundamental spatial geometry container
```r
# Constructor
NeuroSpace(dim, origin = NULL, spacing = NULL, axes = NULL, trans = NULL)

# Key Slots
- dim: integer vector of dimensions
- origin: numeric vector of world coordinates for first voxel
- spacing: numeric vector of voxel sizes
- axes: AxisSet defining orientation
- trans: transformation matrix

# Examples
space_3d <- NeuroSpace(
  dim = c(64L, 64L, 40L),
  spacing = c(2, 2, 2),
  origin = c(-90, -126, -72)
)

# Coordinate transformations
world_coords <- c(0, 0, 0)
vox_idx <- coord_to_index(space_3d, world_coords)
back_to_world <- index_to_coord(space_3d, vox_idx)
```

2. NeuroVol: 3D volumetric data container
```r
# Constructor
NeuroVol(data, space, label = "", indices = NULL)

# Types
- DenseNeuroVol: Full volume storage
- SparseNeuroVol: Masked storage for efficiency

# Examples
# Create from array
vol_data <- array(rnorm(64*64*32), dim=c(64,64,32))
vol <- NeuroVol(vol_data, space, label="example_vol")

# Read from file
vol <- read_vol("image.nii", index=1)

# Write to file
write_vol(vol, "output.nii")
```

3. NeuroVec: 4D data container (3D + time/series)
```r
# Constructor
NeuroVec(data, space = NULL, mask = NULL, label = "")

# Input Types
- 4D array
- Matrix (voxels × time)
- List of NeuroVol objects

# Examples
# Create from 4D array
vec_data <- array(rnorm(64*64*32*100), dim=c(64,64,32,100))
vec <- NeuroVec(vec_data, space)

# Create from NeuroVol list
vol_list <- list(vol1, vol2, vol3)
vec <- NeuroVec(vol_list)

# Create sparse vector with mask
mask <- array(runif(64*64*32) > 0.5, dim=c(64,64,32))
sparse_vec <- NeuroVec(vec_data, space, mask=mask)
```
</core_concepts>

<data_access>
1. Array-like indexing:
```r
# Volume access
vol_value <- vol[x, y, z]
vol_slice <- vol[, , 1]  # First slice

# Vector access
vec_value <- vec[x, y, z, t]
time_point <- vec[, , , 1]  # First timepoint
```

2. Series extraction:
```r
# Get time series at voxel
ts <- series(vec, index)  # Linear index
ts <- series(vec, x, y, z)  # 3D coordinates

# Get multiple time series
mask <- vol > threshold
ts_list <- vectors(vec, mask)

# Get volume sequence
vol_list <- vols(vec)  # All volumes
vol_subset <- vols(vec, 1:10)  # First 10 volumes
```

3. ROI operations:
```r
# Create spherical ROI
roi <- spherical_roi(mask, centroid=c(10,10,10), radius=8)

# Extract ROI time series
roi_ts <- series_roi(vec, coords(roi))

# Convert ROI to sparse vector
sparse_vec <- as(roi, "SparseNeuroVec")
```
</data_access>

<key_operations>
1. Data Loading:
```r
# Basic loading
vol <- read_vol("anatomy.nii")
vec <- read_vec("functional.nii")

# Memory-efficient loading
vec_mmap <- read_vec("large_timeseries.nii", mode="mmap")
vec_masked <- read_vec("functional.nii", mask=mask, mode="bigvec")

# Load multiple files
vec_list <- read_vol_list(c("run1.nii", "run2.nii"))
```

2. Data Manipulation:
```r
# Concatenate vectors
vec_combined <- concat(vec1, vec2)

# Extract subset
sub_vec <- sub_vector(vec, 1:10)

# Resample volume
vol_resampled <- resample(vol, new_space, interpolation=0)

# Map values
vol_mapped <- map_values(vol, lookup=list("0"=2, "1"=5))
```

3. Advanced Features:
```r
# Create 5D data (NeuroHyperVec)
hvec <- NeuroHyperVec(data, space5d, mask)

# Access 5D data
feature_data <- hvec[,,,2,1]  # specific trial and feature
voxel_series <- series(hvec, 5, 5, 5)  # all trials/features for voxel
```
</key_operations>

<best_practices>
1. Memory Management:
- Use "bigvec" mode with mask for large datasets
- Use "mmap" mode for memory-efficient access
- Consider "filebacked" for flexible caching
- Use SparseNeuroVec when working with masked data

2. Data Types:
- DenseNeuroVol/Vec: Full data storage
- SparseNeuroVol/Vec: Masked storage
- ROIVol/Vec: Region analysis
- NeuroHyperVec: 5D data (e.g., multiple trials/features)

3. File Operations:
- Support for .nii and .nii.gz
- Always check file existence
- Verify dimension compatibility
- Use appropriate data types for saving

4. Performance Tips:
- Pre-compute and reuse masks
- Use linear indexing when possible
- Leverage sparse representations
- Consider memory-mapped files for large datasets
</best_practices>
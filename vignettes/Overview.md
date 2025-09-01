---
title: "Overview: Getting Started with neuroim2"
date: "2025-08-31"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Overview: Getting Started with neuroim2}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Introduction

**neuroim2** is a comprehensive R package for neuroimaging data analysis, providing efficient data structures and methods for handling 3D brain volumes and 4D time-series data. Whether you're working with structural MRI, functional MRI, or other volumetric brain imaging data, neuroim2 offers the tools you need.

## Key Features

- **Efficient Data Structures**: Optimized representations for both dense and sparse neuroimaging data
- **Flexible I/O**: Read and write common neuroimaging formats (NIfTI, AFNI)
- **ROI Analysis**: Create and analyze regions of interest with various shapes
- **Searchlight Methods**: Implement searchlight analyses for pattern detection
- **Memory Management**: Handle large datasets with file-backed and memory-mapped arrays
- **Spatial Operations**: Resample, reorient, filter, and transform brain images
- **Parcellation Support**: Work with brain atlases and parcellated data

## Quick Start


``` r
library(neuroim2)

# Load a 3D brain volume
img <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))

# Inspect the spatial properties
space(img)      # Complete spatial information
#> 
#>  NeuroSpace Object 
#> 
#>  >> Dimensions 
#>   Grid Size: 64 x 64 x 25
#>   Memory:   5.9 KB
#> 
#>  >> Spatial Properties 
#>   Spacing:   3.50 x 3.50 x 3.70 mm
#>   Origin:    112.00 x -108.50 x -46.25 mm
#> 
#>  >> Anatomical Orientation 
#>   X: Right-to-Left  |  Y: Posterior-to-Anterior  |  Z: Inferior-to-Superior 
#> 
#>  >> World Transformation 
#>   Forward (Voxel to World): 
#>     -3.500  -0.000  -0.000   112.000
#>  0.000   3.500  -0.000  -108.500
#> -0.000   0.000   3.700   -46.250
#>  0.000   0.000   0.000     1.000 
#>   Inverse (World to Voxel): 
#>     -0.286  -0.000  -0.000  32.000
#>  0.000   0.286   0.000  31.000
#>  0.000   0.000   0.270  12.500
#>  0.000   0.000   0.000   1.000 
#> 
#>  >> Bounding Box 
#>   Min Corner: -108.5, -108.5, -46.2 mm
#>   Max Corner: 112.0, 112.0, 42.6 mm
#> 
#> ==================================================
dim(img)        # Dimensions
#> [1] 64 64 25
spacing(img)    # Voxel sizes
#> [1] 3.5 3.5 3.7
origin(img)     # Origin coordinates
#> [1]  112.00 -108.50  -46.25
```

# Core Data Structures

## 3D Volumes: NeuroVol

The `NeuroVol` class represents 3D brain volumes (structural images, masks, single time points):


``` r
# Create a synthetic 3D volume
dims <- c(64, 64, 40)
dat <- array(rnorm(prod(dims)), dims)
vol <- NeuroVol(dat, NeuroSpace(dims))

# Basic operations
vol_mean <- mean(vol)
vol_sd <- sd(vol)
cat("Volume mean:", vol_mean, "SD:", vol_sd, "\n")
#> Volume mean: 0.001293415 SD: 1.001037

# Logical volumes for masks
mask <- vol > 0
cat("Voxels above zero:", sum(mask), "\n")
#> Voxels above zero: 81901
```

## 4D Time-Series: NeuroVec

The `NeuroVec` class handles 4D data (e.g., fMRI time-series):


``` r
# Create a 4D time-series (small example)
dims_4d <- c(10, 10, 10, 20)  # 10x10x10 volume, 20 time points
dat_4d <- array(rnorm(prod(dims_4d)), dims_4d)
vec_4d <- NeuroVec(dat_4d, NeuroSpace(dims_4d))

# Extract time-series at a specific voxel
ts <- series(vec_4d, 5, 5, 5)
cat("Time-series length at voxel (5,5,5):", length(ts), "\n")
#> Time-series length at voxel (5,5,5): 20

# Extract a single volume at time point 10
vol_t10 <- vec_4d[,,,10]
cat("Volume at t=10 dimensions:", dim(vol_t10), "\n")
#> Volume at t=10 dimensions: 10 10 10
```

# Region of Interest (ROI) Analysis

## Creating ROIs

neuroim2 provides multiple ways to define regions of interest:


``` r
# Spherical ROI - most common for searchlight analyses
sphere <- spherical_roi(img, c(30, 30, 20), radius = 5)
cat("Spherical ROI contains", length(sphere), "voxels\n")
#> Spherical ROI contains 11 voxels

# Cuboid ROI - rectangular box
cube <- cuboid_roi(space(img), c(30, 30, 20), surround = 3)
cat("Cuboid ROI contains", length(cube), "voxels\n")
#> Cuboid ROI contains 343 voxels

# Extract values from the original image using ROI
roi_values <- img[coords(sphere)]
cat("Mean value in spherical ROI:", mean(roi_values), "\n")
#> Mean value in spherical ROI: 1
```

## Searchlight Analysis

Searchlight is a powerful technique for local pattern analysis:


``` r
# Create searchlights with 6mm radius
lights <- searchlight(img, radius = 6, eager = FALSE)

# Process first few searchlights (normally you'd process all)
first_5 <- lights[1:5]
means <- sapply(first_5, function(roi) mean(img[coords(roi)]))
cat("First 5 searchlight means:", round(means, 2), "\n")
#> First 5 searchlight means: 0.95 0.95 0.84 0.74 0.79
```

# Coordinate Systems and Transformations

## Coordinate Conversions

neuroim2 handles conversions between different coordinate systems:


``` r
# Voxel coordinates to world coordinates (mm)
voxel_coords <- c(30, 30, 20)
world_coords <- grid_to_coord(img, matrix(voxel_coords, nrow = 1))
cat("Voxel", voxel_coords, "-> World", round(world_coords, 2), "mm\n")
#> Voxel 30 30 20 -> World 10.5 -7 24.05 mm

# World coordinates back to voxel
voxel_back <- coord_to_grid(img, world_coords)
cat("World", round(world_coords, 2), "-> Voxel", voxel_back, "\n")
#> World 10.5 -7 24.05 -> Voxel 30 30 19.99999

# Linear indices
idx <- grid_to_index(img, matrix(voxel_coords, nrow = 1))
cat("Voxel", voxel_coords, "-> Linear index", idx, "\n")
#> Voxel 30 30 20 -> Linear index 79710
```

# Memory-Efficient Operations

## Sparse Representations

For data with many zero values (e.g., masks, ROIs):


``` r
# Create a sparse representation directly from an ROI
roi <- spherical_roi(img, c(30, 30, 20), radius = 8, fill = 1)

# Convert ROI to sparse volume
sparse_roi <- as.sparse(roi)

# Compare memory usage
dense_vol <- as.dense(sparse_roi)  # Convert back for comparison
#> Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'as.dense' for signature '"SparseNeuroVol"'
cat("Dense size:", format(object.size(dense_vol), units = "auto"), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'format': object 'dense_vol' not found
cat("Sparse size:", format(object.size(sparse_roi), units = "auto"), "\n")
#> Sparse size: 8.3 Kb
cat("Non-zero voxels:", length(roi), "out of", prod(dim(img)), "total\n")
#> Non-zero voxels: 49 out of 102400 total
cat("Space savings:", round((1 - as.numeric(object.size(sparse_roi)) / 
                              as.numeric(object.size(dense_vol))) * 100, 1), "%\n")
#> Error: object 'dense_vol' not found
```

## File-Backed Arrays

For datasets too large to fit in memory:


``` r
# Create a file-backed 4D dataset (not run in vignette)
big_vec <- FileBackedNeuroVec(
  dims = c(91, 109, 91, 1000),  # Standard MNI space, 1000 volumes
  dtype = "float32",
  file_name = "big_data.dat"
)

# Access works like regular arrays but data stays on disk
subset <- big_vec[45, 55, 45, 1:10]  # Load only what you need
```

# Spatial Filtering and Processing

## Smoothing Operations


``` r
# Gaussian smoothing with single sigma value
img_smooth <- gaussian_blur(img, sigma = 2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'space': argument "mask" is missing, with no default

# Compare original vs smoothed
orig_vals <- img[30:32, 30:32, 20]
smooth_vals <- img_smooth[30:32, 30:32, 20]
#> Error: object 'img_smooth' not found
cat("Original variance:", var(as.vector(orig_vals)), "\n")
#> Original variance: 0
cat("Smoothed variance:", var(as.vector(smooth_vals)), "\n")
#> Error: object 'smooth_vals' not found
```

## Resampling


``` r
# Downsample by factor of 2
img_down <- downsample(img, c(2, 2, 2))
cat("Original dimensions:", dim(img), "\n")
#> Original dimensions: 64 64 25
cat("Downsampled dimensions:", dim(img_down), "\n")
#> Downsampled dimensions: 112 112 46
```

# Working with Parcellations

## ClusteredNeuroVol for Atlas-Based Analysis


``` r
# Create a simple parcellation
coords <- index_to_coord(img, which(as.logical(img)))
set.seed(123)
k <- 10  # 10 parcels
if (nrow(coords) > k) {
  km <- kmeans(coords, centers = k, iter.max = 100)
  
  # Create clustered volume
  cvol <- ClusteredNeuroVol(img, km$cluster)
  cat("Created", num_clusters(cvol), "parcels\n")
  
  # Get centroids of each parcel
  centers <- centroids(cvol)
  cat("First parcel centroid:", round(centers[1,], 1), "mm\n")
}
#> Created 10 parcels
#> First parcel centroid: 41.3 26.2 18.8 mm
```

# Input/Output Operations

## Reading and Writing


``` r
# Write a volume to a temporary file
tmp_file <- tempfile(fileext = ".nii.gz")
write_vol(img, tmp_file)
cat("Wrote volume to:", tmp_file, "\n")
#> Wrote volume to: /var/folders/9h/nkjq6vss7mqdl4ck7q1hd8ph0000gp/T//RtmpKXhOQ3/filebe04d43a54d.nii.gz

# Read it back
img_read <- read_vol(tmp_file)
cat("Read volume with dimensions:", dim(img_read), "\n")
#> Read volume with dimensions: 64 64 25

# Clean up
file.remove(tmp_file)
#> [1] TRUE
```

## Working with Multiple Files


``` r
# Read multiple volumes (not run)
files <- c("scan1.nii", "scan2.nii", "scan3.nii")
vols <- read_vec(files)  # Creates a NeuroVecSeq

# Or read as a single concatenated 4D volume
vec_concat <- read_vec(files, concatenate = TRUE)
```

# Practical Example: ROI-Based Time-Series Analysis

Here's a complete workflow combining multiple features:


``` r
# Create synthetic fMRI-like data
dims_fmri <- c(20, 20, 15, 50)  # Small for example
fmri_data <- array(rnorm(prod(dims_fmri), mean = 1000, sd = 50), dims_fmri)
fmri <- NeuroVec(fmri_data, NeuroSpace(dims_fmri))

# Define an ROI
roi <- spherical_roi(space(fmri), c(10, 10, 8), radius = 3)
cat("ROI size:", length(roi), "voxels\n")
#> ROI size: 123 voxels

# Extract time-series from ROI
roi_ts <- series_roi(fmri, roi)
roi_mat <- values(roi_ts)  # T x N matrix
cat("ROI time-series matrix:", dim(roi_mat), "\n")
#> ROI time-series matrix: 50 123

# Compute mean time-series
mean_ts <- rowMeans(roi_mat)

# Z-score the mean time-series
z_ts <- scale(mean_ts)
#> Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'scale' for signature '"numeric"'
cat("Mean time-series - Mean:", round(mean(z_ts), 4), 
    "SD:", round(sd(z_ts), 4), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'mean': object 'z_ts' not found

# Find peak activation time
peak_time <- which.max(z_ts)
#> Error: object 'z_ts' not found
cat("Peak activation at time point:", peak_time, 
    "with z-score:", round(z_ts[peak_time], 2), "\n")
#> Error: object 'peak_time' not found
```

# Performance Tips

1. **Use appropriate data structures**: 
   - `SparseNeuroVol` for masks and ROIs
   - `FileBackedNeuroVec` for very large datasets
   - Regular arrays for small-to-medium data

2. **Vectorize operations**: 
   - Use `spherical_roi_set()` instead of loops for multiple ROIs
   - Process searchlights in parallel when possible

3. **Preallocate memory**: 
   - Know your output dimensions and allocate arrays upfront

4. **Choose the right format**: 
   - Use compressed NIfTI (.nii.gz) for storage
   - Keep working data uncompressed for speed

# Getting Help

## Package Documentation

``` r
# List all functions
help(package = "neuroim2")

# Search for specific topics
help.search("roi", package = "neuroim2")
```

## Vignettes for Deep Dives

- **3D Volumes**: `vignette("ImageVolumes", package = "neuroim2")`
- **4D Time-Series**: `vignette("NeuroVector", package = "neuroim2")`  
- **ROI Analysis**: `vignette("regionOfInterest", package = "neuroim2")`
- **Parcellations**: `vignette("clustered-neurovec", package = "neuroim2")`
- **Pipelines**: `vignette("pipelines", package = "neuroim2")`

## Quick Reference

Common operations at a glance:

| Task | Function |
|------|----------|
| Read NIfTI file | `read_vol()`, `read_vec()` |
| Write NIfTI file | `write_vol()`, `write_vec()` |
| Create spherical ROI | `spherical_roi()` |
| Extract time-series | `series()`, `series_roi()` |
| Smooth image | `gaussian_blur()` |
| Resample image | `resample()` |
| Coordinate conversion | `coord_to_grid()`, `grid_to_coord()` |
| Create mask | `LogicalNeuroVol()` |
| Sparse representation | `as.sparse()` |
| Concatenate volumes | `concat()` |

# Summary

neuroim2 provides a comprehensive toolkit for neuroimaging analysis in R. Its efficient data structures, flexible ROI tools, and memory-conscious design make it suitable for both interactive exploration and large-scale processing pipelines. Start with the basic NeuroVol and NeuroVec classes, explore ROI creation for your specific needs, and leverage sparse or file-backed arrays when working with large datasets.

For more advanced usage and specific workflows, consult the topic-specific vignettes listed above.

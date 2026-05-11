# Advanced ROI Construction

This article is now the detailed ROI reference to read after
[`vignette("AnalysisWorkflows")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md).
The core workflow vignette shows when to use ROIs; this article goes
deeper into the different ROI constructors and ROI data structures
available in `neuroim2`.

### ROI types in neuroim2

The package supports several ROI styles, depending on whether you want
one simple neighborhood, a grid-based box, a parcel map, or a weighted
kernel:

| ROI Type | Function | Description | Use Case |
|----|----|----|----|
| **Spherical** | [`spherical_roi()`](https://bbuchsbaum.github.io/neuroim2/reference/spherical_roi.md) | Sphere around a center point | Searchlight analysis, seed-based connectivity |
| **Cuboid** | [`cuboid_roi()`](https://bbuchsbaum.github.io/neuroim2/reference/cuboid_roi.md) | Rectangular box region | Grid-based parcellation |
| **Square** | [`square_roi()`](https://bbuchsbaum.github.io/neuroim2/reference/square_roi.md) | 2D square in one plane | Slice-based analysis |
| **Clustered** | `ClusteredNeuroVol` | Parcellation-based regions | Atlas-based analysis |
| **Custom** | `Kernel` | Weighted regions | Gaussian-weighted, gradient-based |

### Quick start

``` r

file_name <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
vol <- read_vol(file_name)

roi <- spherical_roi(vol, c(20, 20, 20), radius = 5)
cat("ROI contains", length(roi), "voxels\n")
#> ROI contains 11 voxels
```

------------------------------------------------------------------------

## Basic ROI operations

The introductory ROI extraction and searchlight workflow now lives in
[`vignette("AnalysisWorkflows")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md).
This article stays focused on ROI construction, representation, and less
common ROI variants.

### Creating a Spherical ROI

Spherical ROIs are still the base constructor most other workflows build
on, so this article starts there and then moves into less common
patterns.

``` r

# Load an example brain volume
file_name <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
vol <- read_vol(file_name)

# Create a spherical ROI centered at voxel coordinates [20,20,20]
# with a 5mm radius, filling all values with 100
sphere <- spherical_roi(vol, c(20, 20, 20), radius = 5, fill = 100)

# Examine the ROI
cat("Number of voxels in ROI:", length(sphere), "\n")
#> Number of voxels in ROI: 11
cat("ROI dimensions:", dim(sphere), "\n")
#> ROI dimensions: 11 3
cat("All values are 100:", all(sphere == 100), "\n")
#> All values are 100: TRUE
```

#### Performance Note: C++ vs Pure R Implementation

``` r

# The spherical_roi function can use either C++ (fast) or pure R (slower)
# C++ is the default and recommended for performance
sphere_cpp <- spherical_roi(vol, c(20, 20, 20), radius = 5, use_cpp = TRUE)
sphere_r   <- spherical_roi(vol, c(20, 20, 20), radius = 5, use_cpp = FALSE)

# They produce the same voxel set (up to ordering and type)
all.equal(
  coords(sphere_cpp)[order(coords(sphere_cpp)[,1], coords(sphere_cpp)[,2], coords(sphere_cpp)[,3]),],
  coords(sphere_r)[order(coords(sphere_r)[,1], coords(sphere_r)[,2], coords(sphere_r)[,3]),],
  check.attributes = FALSE
)
#> [1] "Mean relative difference: 0.05"
```

### Creating a Spherical ROI Around Real-World Coordinates

Often, you’ll have coordinates in millimeter space (e.g., from published
studies) rather than voxel indices.

``` r

# Define a real-world coordinate in mm
rpoint <- c(-34, -28, 10)

# Convert real-world coordinates to voxel coordinates
vox <- coord_to_grid(vol, rpoint)
cat("Real coordinate:", rpoint, "-> Voxel coordinate:", vox, "\n")
#> Real coordinate: -34 -28 10 -> Voxel coordinate: 42.71429 24 16.2027

# Create spherical ROI with 10mm radius
sphere <- spherical_roi(vol, vox, radius = 10, fill = 1)
cat("ROI contains", length(sphere), "voxels\n")
#> ROI contains 85 voxels

# Verify the center of mass is close to our target
coords <- index_to_coord(vol, indices(sphere))
center_of_mass <- colMeans(coords)
cat("Center of mass:", center_of_mass, "\n")
#> Center of mass: -33.25 -29.75 11.1
cat("Original point:", rpoint, "\n")
#> Original point: -34 -28 10
```

### Creating Multiple Spherical ROIs Efficiently

When creating many ROIs, use the vectorized
[`spherical_roi_set()`](https://bbuchsbaum.github.io/neuroim2/reference/spherical_roi_set.md)
function for better performance:

``` r

library(neuroim2)
vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
d <- dim(vol)

# Define multiple ROI centers
centers <- rbind(
  c(floor(d[1]/3), floor(d[2]/3), floor(d[3]/3)),
  c(floor(d[1]/2), floor(d[2]/2), floor(d[3]/2)),
  c(floor(2*d[1]/3), floor(2*d[2]/3), floor(2*d[3]/3))
)

# Efficient vectorized creation
rois <- spherical_roi_set(bvol = vol, centroids = centers, radius = 5, fill = 1)
cat("Created", length(rois), "ROIs\n")
#> Created 3 ROIs

# Compare with loop approach (slower but equivalent)
rois2 <- lapply(seq_len(nrow(centers)), function(i) {
  spherical_roi(vol, centers[i,], 5, fill = 1)
})
cat("Loop method also created", length(rois2), "ROIs\n")
#> Loop method also created 3 ROIs
```

### Creating Cuboid and Square ROIs

#### Cuboid ROIs (3D boxes)

``` r

# Create a cuboid ROI - a 3D rectangular box
sp1 <- NeuroSpace(c(20, 20, 20), c(1, 1, 1))

# Create a 7x7x7 cube (surround=3 means 3 voxels on each side of center)
cube <- cuboid_roi(sp1, c(10, 10, 10), surround = 3, fill = 5)
cat("Cuboid ROI contains", length(cube), "voxels\n")
#> Cuboid ROI contains 343 voxels
cat("All values are 5:", all(cube == 5), "\n")
#> All values are 5: TRUE

# Get the coordinates
vox_coords <- coords(cube)
cat("Coordinate ranges:\n")
#> Coordinate ranges:
cat("  X:", range(vox_coords[,1]), "\n")
#>   X: 7 13
cat("  Y:", range(vox_coords[,2]), "\n")
#>   Y: 7 13
cat("  Z:", range(vox_coords[,3]), "\n")
#>   Z: 7 13
```

#### Square ROIs (2D in 3D space)

``` r

# Create a square ROI - a 2D square fixed in one dimension
sp1 <- NeuroSpace(c(20, 20, 20), c(1, 1, 1))

# Create a 5x5 square in the z=10 plane (fixdim=3 fixes the z dimension)
square <- square_roi(sp1, c(10, 10, 10), surround = 2, fill = 3, fixdim = 3)
cat("Square ROI contains", length(square), "voxels (should be 25)\n")
#> Square ROI contains 25 voxels (should be 25)

# Verify it's planar
vox_coords <- coords(square)
cat("All z-coordinates are the same:", 
    length(unique(vox_coords[,3])) == 1, "\n")
#> All z-coordinates are the same: TRUE
cat("Z-plane:", unique(vox_coords[,3]), "\n")
#> Z-plane: 10
```

------------------------------------------------------------------------

## Working with ROI Data

### Converting ROIs to Sparse Volumes

Sparse volumes are memory-efficient representations that only store
non-zero values:

``` r

# Create a spherical ROI
sphere <- spherical_roi(vol, c(50, 10, 10), radius = 10, fill = 1)
cat("ROI contains", length(sphere), "voxels\n")
#> ROI contains 85 voxels

# Convert to SparseNeuroVol - memory efficient storage
# Prefer the provided coercion helper
sparsevol <- as.sparse(sphere)

# Verify properties
cat("Sum of values match:", sum(sparsevol) == sum(sphere), "\n")
#> Sum of values match: TRUE
cat("Dimensions match:", all(dim(sparsevol) == dim(vol)), "\n")
#> Dimensions match: TRUE

# Memory comparison
roi_size <- object.size(sphere)
sparse_size <- object.size(sparsevol)
cat("Memory usage - ROI:", format(roi_size, units = "auto"), "\n")
#> Memory usage - ROI: 9.8 Kb
cat("Memory usage - Sparse:", format(sparse_size, units = "auto"), "\n")
#> Memory usage - Sparse: 8.9 Kb
```

The basic 4D ROI extraction workflow with
[`series_roi()`](https://bbuchsbaum.github.io/neuroim2/reference/series_roi.md)
now lives in
[`vignette("AnalysisWorkflows")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md).
The current article keeps the lower-level ROI data structures and
construction details.

### Using ROIVec for 4D ROI Data

The `ROIVec` class is designed for efficient storage of 4D ROI data:

``` r

# Create a 4D NeuroSpace
vspace <- NeuroSpace(dim = c(10, 10, 10, 20), spacing = c(1, 1, 1))

# Define ROI coordinates
roi_coords <- matrix(c(5,5,5, 6,5,5, 5,6,5, 5,5,6), ncol = 3, byrow = TRUE)

# Create synthetic time-series data for each voxel
n_timepoints <- 20
n_voxels <- nrow(roi_coords)
roi_data <- matrix(rnorm(n_timepoints * n_voxels), 
                   nrow = n_timepoints, ncol = n_voxels)

# Create ROIVec object
roi_vec <- ROIVec(vspace, roi_coords, roi_data)
cat("ROIVec created with", nrow(roi_coords), "voxels and", 
    nrow(roi_data), "timepoints\n")
#> ROIVec created with 4 voxels and 20 timepoints

# Access as matrix of values (T x N)
roi_matrix <- values(roi_vec)
cat("Matrix dimensions (T x N):", dim(roi_matrix), "\n")
#> Matrix dimensions (T x N): 20 4
```

Searchlight workflows, including
[`searchlight()`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight.md),
[`random_searchlight()`](https://bbuchsbaum.github.io/neuroim2/reference/random_searchlight.md),
and
[`clustered_searchlight()`](https://bbuchsbaum.github.io/neuroim2/reference/clustered_searchlight.md),
are now covered in
[`vignette("AnalysisWorkflows")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md).
That keeps the current article focused on constructing ROIs and ROI-like
data objects rather than on the first-pass analyses built on top of
them.

------------------------------------------------------------------------

## Image Patches

Patches are regular, fixed-size regions that tile the image space,
useful for convolutional approaches:

``` r

# Create 3x3x1 patches covering the volume
pset <- patch_set(vol, dims = c(3, 3, 1))
cat("Number of patches:", length(pset), "\n")
#> Number of patches: 102400

# Compute statistics for each patch
patch_means <- pset %>% purrr::map_dbl(~ mean(.))
cat("Patch mean range:", range(patch_means, na.rm = TRUE), "\n")
#> Patch mean range: 0 1

# Restrict patches to masked regions only
pset_masked <- patch_set(vol, dims = c(3, 3, 1), mask = as.logical(vol))
cat("Masked patches:", length(pset_masked), "\n")
#> Masked patches: 29532

# Patches are guaranteed to be equal size (padded at edges if needed)
patch_sizes <- pset_masked %>% purrr::map_int(~ length(.))
cat("All patches same size:", length(unique(patch_sizes)) == 1, "\n")
#> All patches same size: TRUE
cat("Patch size:", unique(patch_sizes), "voxels\n")
#> Patch size: 18 voxels
```

------------------------------------------------------------------------

## Advanced ROI Techniques

### Custom Weighted ROIs with Kernels

Create Gaussian-weighted or other custom-shaped ROIs:

``` r

# Create a Gaussian kernel for weighted ROI
kern_dim <- c(5, 5, 5)  # 5x5x5 kernel
voxel_dim <- c(1, 1, 1)  # 1mm isotropic voxels

# Create Gaussian-weighted kernel
gauss_kernel <- Kernel(kerndim = kern_dim, vdim = voxel_dim, 
                       FUN = dnorm, sd = 1.5)

# Embed kernel at a specific location
embedded <- embed_kernel(gauss_kernel, space(vol), 
                         center_voxel = c(20, 20, 20))
cat("Kernel weights sum to:", sum(embedded), "\n")
#> Kernel weights sum to: 1
```

### Working with ClusteredNeuroVec for Parcellated Data

Efficiently handle parcellated 4D data where voxels are grouped into
regions:

``` r

# Create synthetic 4D data
sp4 <- NeuroSpace(c(10, 10, 10, 20), c(1, 1, 1))
arr <- array(rnorm(10*10*10*20), dim = c(10, 10, 10, 20))
vec <- NeuroVec(arr, sp4)

# Create mask for central region
sp3 <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
mask_arr <- array(FALSE, dim = c(10, 10, 10))
mask_arr[3:8, 3:8, 3:8] <- TRUE
mask <- LogicalNeuroVol(mask_arr, sp3)

# Assign voxels to 5 random clusters
n_voxels <- sum(mask_arr)
clusters <- sample(1:5, n_voxels, replace = TRUE)
cvol <- ClusteredNeuroVol(mask, clusters)

# Create clustered representation - one time-series per cluster
cv <- ClusteredNeuroVec(vec, cvol)

# Access cluster time-series efficiently
cluster_matrix <- as.matrix(cv)  # T x K matrix
cat("Cluster matrix dimensions:", dim(cluster_matrix), "\n")
#> Cluster matrix dimensions: 20 5
cat("(", dim(cluster_matrix)[1], "timepoints x", 
    dim(cluster_matrix)[2], "clusters)\n")
#> ( 20 timepoints x 5 clusters)
```

### ROI Set Operations

Combine and manipulate multiple ROIs:

``` r

# Create two overlapping spherical ROIs
roi1 <- spherical_roi(vol, c(20, 20, 20), radius = 6, fill = 1)
roi2 <- spherical_roi(vol, c(23, 20, 20), radius = 6, fill = 1)

# Get indices for set operations
idx1 <- indices(roi1)
idx2 <- indices(roi2)

# Intersection - voxels in both ROIs
intersection_idx <- intersect(idx1, idx2)
cat("Intersection:", length(intersection_idx), "voxels\n")
#> Intersection: 0 voxels

# Union - voxels in either ROI
union_idx <- union(idx1, idx2)
cat("Union:", length(union_idx), "voxels\n")
#> Union: 38 voxels

# Difference - voxels in roi1 but not roi2
diff_idx <- setdiff(idx1, idx2)
cat("Difference (roi1 - roi2):", length(diff_idx), "voxels\n")
#> Difference (roi1 - roi2): 19 voxels

# Calculate overlap percentage
overlap_pct <- length(intersection_idx) / length(union_idx) * 100
cat("Overlap:", round(overlap_pct, 1), "%\n")
#> Overlap: 0 %
```

------------------------------------------------------------------------

## Best Practices and Performance

### Memory Management

When working with many ROIs, consider memory usage:

``` r

# Compare memory usage of different ROI storage methods
n_rois <- 100
centers <- matrix(runif(n_rois * 3, 5, 15), ncol = 3)

# Method 1: List of ROI objects (flexible but more memory)
  roi_list <- lapply(1:nrow(centers), function(i) {
  spherical_roi(vol, centers[i,], radius = 6, fill = 1)
  })

# Method 2: Single sparse volume with labeled regions (ensure increasing indices)
all_indices <- list()
all_labels <- list()
for (i in 1:nrow(centers)) {
  roi <- spherical_roi(vol, centers[i,], radius = 6)
  idx <- indices(roi)
  idx <- sort(unique(idx))
  all_indices[[i]] <- idx
  all_labels[[i]] <- rep(i, length(idx))
}
idx_all <- unlist(all_indices)
lab_all <- unlist(all_labels)
ord <- order(idx_all)
idx_all <- idx_all[ord]
lab_all <- lab_all[ord]
# Ensure strictly increasing indices by removing duplicates
keep <- !duplicated(idx_all)
idx_all <- idx_all[keep]
lab_all <- lab_all[keep]
combined_sparse <- SparseNeuroVol(lab_all, space(vol), indices = idx_all)

cat("Memory usage:\n")
#> Memory usage:
cat("  ROI list:", format(object.size(roi_list), units = "auto"), "\n")
#>   ROI list: 771.1 Kb
cat("  Sparse combined:", format(object.size(combined_sparse), units = "auto"), "\n")
#>   Sparse combined: 23.2 Kb
```

### Choosing ROI Sizes

Guidelines for selecting appropriate ROI sizes:

``` r

# Demonstrate effect of ROI size on coverage and overlap
radii <- c(6, 8, 10, 12)
coverage_stats <- data.frame(
  radius = radii,
  n_voxels = numeric(length(radii)),
  pct_overlap = numeric(length(radii))
)

center1 <- c(20, 20, 20)
center2 <- c(25, 20, 20)  # 5 voxels apart

for (i in seq_along(radii)) {
  roi1 <- spherical_roi(vol, center1, radius = radii[i])
  roi2 <- spherical_roi(vol, center2, radius = radii[i])
  
  coverage_stats$n_voxels[i] <- length(roi1)
  
  overlap <- length(intersect(indices(roi1), indices(roi2)))
  total <- length(union(indices(roi1), indices(roi2)))
  coverage_stats$pct_overlap[i] <- overlap / total * 100
}

print(coverage_stats)
#>   radius n_voxels pct_overlap
#> 1      6       19    0.000000
#> 2      8       49    0.000000
#> 3     10       85    0.000000
#> 4     12      163    5.844156
```

### Parallel Processing Tips

For large-scale ROI analyses, consider parallel processing:

``` r

# Example of parallel searchlight analysis (not run in vignette)
library(parallel)
library(foreach)
library(doParallel)

# Setup parallel backend
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parallel searchlight computation
searchlight_results <- foreach(roi = searchlight_list, 
                               .packages = c("neuroim2")) %dopar% {
  # Your analysis function here
  mean(data[coords(roi)])
}

# Clean up
stopCluster(cl)
```

------------------------------------------------------------------------

## Troubleshooting and Tips

### Common Issues and Solutions

#### Issue: ROI extends outside brain mask

``` r

# Problem: ROI at edge of brain
edge_vox <- c(2, 2, 2)  # Near edge of volume

# Solution 1: Use nonzero=TRUE to keep only mask voxels
roi_filtered <- spherical_roi(vol, edge_vox, radius = 5, nonzero = TRUE)
cat("Filtered ROI size:", length(roi_filtered), "voxels\n")
#> Filtered ROI size: 0 voxels

# Solution 2: Check if ROI is valid before analysis
if (length(roi_filtered) < 10) {
  cat("Warning: ROI too small for reliable analysis\n")
}
#> Warning: ROI too small for reliable analysis
```

#### Issue: Different coordinate systems

``` r

# Always verify your coordinate system
test_coord_mm <- c(-34, -28, 10)  # MNI coordinates
test_coord_vox <- coord_to_grid(vol, test_coord_mm)

# Verify round-trip conversion
back_to_mm <- grid_to_coord(vol, matrix(test_coord_vox, nrow = 1))
cat("Original (mm):", test_coord_mm, "\n")
#> Original (mm): -34 -28 10
cat("Voxel:", test_coord_vox, "\n")  
#> Voxel: 42.71429 23.85714 16.18919
cat("Back to mm:", back_to_mm, "\n")
#> Back to mm: -34 -28.00001 10.00001
```

------------------------------------------------------------------------

## See Also

For related functionality, see these other vignettes and help pages:

- [`vignette("NeuroVector")`](https://bbuchsbaum.github.io/neuroim2/articles/NeuroVector.md) -
  Working with 4D neuroimaging data
- [`?read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
  and
  [`?read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md) -
  Reading neuroimaging files
- [`?NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace.md) -
  Understanding coordinate systems and spaces
- [`?ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md) -
  Parcellation-based analyses
- [`?SparseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md) -
  Memory-efficient sparse representations
- [`?series`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md)
  and
  [`?series_roi`](https://bbuchsbaum.github.io/neuroim2/reference/series_roi.md) -
  Time-series extraction

For the complete list of ROI-related functions:

``` r

help(package = "neuroim2", topic = "roi")
```

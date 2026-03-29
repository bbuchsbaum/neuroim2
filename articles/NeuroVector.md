# Working with 4D NeuroVectors (NeuroVec)

## Working with neuroimaging time-series data

The `neuroim2` package contains data structures and functions for
reading, accessing, and processing 4-dimensional neuroimaging data.

In this vignette, we’ll introduce `NeuroVec` (4D images) and related
helpers you’ll use most often:

- Read/write on-disk images (`read_vec`, `write_vec`)
- Spatial metadata via `NeuroSpace` (dimensions, spacing, origin)
- Voxel- and ROI-based access (`series`, `series_roi`, `coords`)
- Common transforms (z-scoring with `scale_series`, concatenation with
  `concat`)
- Dense vs. sparse representations (`DenseNeuroVec`, `SparseNeuroVec`,
  `as.sparse`)

### Reading a four-dimensional NifTI image with read_vec

Here we read in an example image. This file is 4D in the package data
(64×64×25×4), so
[`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
returns a `NeuroVec` with four volumes (timepoints).

``` r
      library(purrr)
      library(ggplot2)
      file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
      vec <- read_vec(file_name)
      dim(vec)
#> [1] 64 64 25  4
      vec
#> <DenseNeuroVec> [3.1 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25 (4 timepoints)
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.288 +/- 0.453 (t=1)
#>   Label         : /home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii
```

Now imagine we have a set of images. We can read multiple files with
`read_vec`. Passing multiple paths returns a `NeuroVecSeq` (a sequence
of vectors) rather than a single concatenated 4D vector.

``` r
    
      file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
      vec_multi <- read_vec(c(file_name, file_name, file_name))
      dim(vec_multi)
#> [1] 64 64 25 12
      
      vec2 <- read_vec(rep(file_name, 10))
      vec2
#> <NeuroVecSeq> [10 vectors] 
#> ── Sequence ──────────────────────────────────────────────────────────────────── 
#>     [1]         : DenseNeuroVec 64x64x25 x 4t
#>     [2]         : DenseNeuroVec 64x64x25 x 4t
#>     [3]         : DenseNeuroVec 64x64x25 x 4t
#>     [4]         : DenseNeuroVec 64x64x25 x 4t
#>     [5]         : DenseNeuroVec 64x64x25 x 4t
#>   ... and 5 more
```

To extract a subset of volumes from a 4D vector, use `sub_vector`
directly:

``` r
      # Extract a subset of volumes (first 3 timepoints)
      vec_1_3 <- sub_vector(vec, 1:3)
      dim(vec_1_3)
#> [1] 64 64 25  3
      vec_1_3
#> <DenseNeuroVec> [2.4 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25 (3 timepoints)
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.288 +/- 0.453 (t=1)
#>   Label         : none
```

### Extracting time-series data using the `series` and `series_roi` functions

To get the time-series at voxel (1,1,1) we can use the `series`
function:

``` r
      
      series(vec_1_3, 1,1,1)
#> [1] 0 0 0
```

We can extract a 4d region of interest with the `series_roi` as follows:

``` r
      file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
      vol <- read_vol(file_name)
      roi <- spherical_roi(vol, c(12,12,12), radius=8)
      rvec1 <- series_roi(vec, roi)
      
      ## or alternatively as a pipeline
      rvec2 <- read_vol(file_name) %>% spherical_roi(c(12,12,12), radius=8) %>% series_roi(vec,.)
      rvec2
#> <ROIVec> [49 voxels x 4 features] 
#> ── Structure ─────────────────────────────────────────────────────────────────── 
#>   Voxels        : 49
#>   Features      : 4
#>   Size          : 9.6 Kb
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 0.000]
      
      ## we can extract the ROI values with the `values` method.
      stopifnot(all(values(rvec1) == values(rvec2)))
      stopifnot(all(coords(rvec1) == coords(rvec2)))
      
```

We can also extract an ROI using 1d indices:

``` r

r1 <- series_roi(vec, 1:100)
r1
#> <ROIVec> [100 voxels x 4 features] 
#> ── Structure ─────────────────────────────────────────────────────────────────── 
#>   Voxels        : 100
#>   Features      : 4
#>   Size          : 11.2 Kb
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 0.000]
```

Or we can extract a plain matrix using the `series` function:

``` r
r2 <- series(vec, 1:100)
dim(r2)
#> [1]   4 100
```

We can also use coordinate indexing using voxel coordinates. First we
load a binary mask with the same spatial dimensions as our NeuroVec:

``` r
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
```

Now we convert indices to voxels and extract a matrix of values at the
specified locations:

``` r
vox <- index_to_grid(mask, 1:100)

r3 <- series(vec, vox)
dim(r3)
#> [1]   4 100
```

And the same using `series_roi`:

``` r
r4 <- series_roi(vec,vox)
r4
#> <ROIVec> [100 voxels x 4 features] 
#> ── Structure ─────────────────────────────────────────────────────────────────── 
#>   Voxels        : 100
#>   Features      : 4
#>   Size          : 12.4 Kb
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 0.000]
```

### Inspecting spatial metadata with NeuroSpace

Every `NeuroVec` carries a `NeuroSpace` describing its geometry.

``` r
sp <- space(vec)
sp                   # dimensions, spacing, origin, axes
#> <NeuroSpace> [4D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25 x 4
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#>   Voxels        : 409,600
dim(vec)             # 4D dims: X × Y × Z × T
#> [1] 64 64 25  4
spacing(vec)         # voxel spacing (mm)
#> [1] 3.5 3.5 3.7
origin(vec)          # image origin
#> [1]  112.0 -108.0  -46.2
ndim(vec)            # == 4 for time series
#> [1] 4
```

The default mask for dense vectors is “all voxels are valid”:

``` r
m <- mask(vec)
m
#> <DenseNeuroVol> [406.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [1.000, 1.000]
```

### Creating a NeuroVec from an in-memory array

You can build a `NeuroVec` directly from arrays or matrices:

``` r
set.seed(1)
dims <- c(24, 24, 24, 5)
arr  <- array(rnorm(prod(dims)), dims)
sp4  <- NeuroSpace(dims, spacing = c(2,2,2))
dvec <- NeuroVec(arr, sp4)
dim(dvec)
#> [1] 24 24 24  5
```

You can also start from a matrix (voxels × time or time × voxels) using
`DenseNeuroVec`:

``` r
mat  <- matrix(rnorm(prod(dims)), nrow = prod(dims[1:3]), ncol = dims[4])
dvec2 <- DenseNeuroVec(mat, sp4)
all.equal(dim(dvec), dim(dvec2))
#> [1] TRUE
```

### Time-series transforms: z-scoring and summary volumes

Z-score each voxel’s time-series (center and scale across time):

``` r
vec_z <- scale_series(dvec, center = TRUE, scale = TRUE)
dim(vec_z)
#> [1] 24 24 24  5
```

Compute a mean volume across time and return a 3D `NeuroVol`:

``` r
M      <- as.matrix(dvec)              # voxels × time
vmean  <- rowMeans(M)                  # per-voxel mean
mean3d <- NeuroVol(vmean, drop_dim(space(dvec)))
mean3d
#> <DenseNeuroVol> [114.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 24 x 24 x 24
#>   Spacing       : 2 x 2 x 2 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [-1.643, 1.756]
```

### Concatenating along time

Append time points by concatenating vectors or volumes:

``` r
dvec_more <- concat(dvec, dvec)        # doubles the time dimension
dim(dvec_more)
#> [1] 24 24 24 10
```

### Dense ↔︎ sparse workflows

Sparse representations store only voxels within a mask. This is handy
for large ROIs or brain masks.

``` r
# Build a random mask and convert a dense vec to sparse
mask_arr <- array(runif(prod(dims[1:3])) > 0.7, dims[1:3])
mask_vol <- LogicalNeuroVol(mask_arr, drop_dim(sp4))

svec <- as.sparse(dvec, mask_vol)     # SparseNeuroVec with explicit mask
svec                                 # note the stored mask and cardinality
#> <SparseNeuroVec> [306 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 24 x 24 x 24
#>   Time Points   : 5
#>   Spacing       : 2 x 2 x 2
#>   Origin        : 0, 0, 0
#> ── Sparse ────────────────────────────────────────────────────────────────────── 
#>   Cardinality   : 4128

# Convert back to dense if needed
dense_again <- as.dense(svec)
all.equal(dim(dense_again), dim(dvec))
#> [1] TRUE
```

Tip: For file-backed or memory-mapped vectors, convert to
`DenseNeuroVec` via a matrix if you need dense-only operations:

``` r
dv_dense <- DenseNeuroVec(as.matrix(vec), space(vec))
```

### Writing vectors to disk

You can write `NeuroVec` and `NeuroVol` objects as NIfTI files:

``` r
tmp_vec <- tempfile(fileext = ".nii.gz")
write_vec(vec_1_3, tmp_vec)
file.exists(tmp_vec)
#> [1] TRUE
unlink(tmp_vec)
```

### Putting it together with an ROI

Combine ROI extraction with time-series transforms:

``` r
vol3d <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
roi   <- spherical_roi(vol3d, c(12,12,12), radius = 6)
rts   <- series_roi(vec, roi)          # ROIVec (T × N with coords)

# z-score each column (voxel) across time, then average within ROI
mat_roi  <- values(rts)                # T × N
mat_z    <- base::scale(mat_roi, center=TRUE, scale=TRUE)
roi_mean <- rowMeans(mat_z)
length(roi_mean)                       # matches time dimension
#> [1] 4
```

That’s the core workflow for 4D data in `neuroim2`: load or create a
`NeuroVec`, inspect its `NeuroSpace`, access time-series via voxels or
ROIs, apply simple transforms, and optionally move between dense and
sparse representations.

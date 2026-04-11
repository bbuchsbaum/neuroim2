# Advanced NeuroVec Patterns

This article is now the deeper companion to
[`vignette("VolumesAndVectors")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md).
Read that first if you want the compact introduction to `NeuroVec`. Stay
here when you want more detail on multi-file reads, manual construction,
matrix views, and dense/sparse conversions for 4D data.

## Working with neuroimaging time-series data

This vignette focuses on the parts of the 4D workflow that remain useful
after you already understand the basic container model:

- reading many 4D files into a `NeuroVecSeq`
- creating `NeuroVec` objects from arrays or matrices
- time-series transforms and concatenation
- dense/sparse conversion patterns
- writing vectors back to disk

### Reading multiple four-dimensional images

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vec <- read_vec(file_name)
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
write_vec(sub_vector(vec, 1:3), tmp_vec)
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

This article now complements the core path rather than replacing it. Use
[`vignette("VolumesAndVectors")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md)
for the basic 4D container model and come back here when you need denser
transformation and conversion patterns.

# Working with 3D Image Volumes

## Reading a NIFTI formatted image volume

The way to read an volumetric image file is to use `read_vol`:

``` r
    suppressPackageStartupMessages(library(neuroim2))
    file_name <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
    vol <- read_vol(file_name)
```

## Working with image volumes

Information about the geometry of the image volume is shown here:

``` r
    print(vol)
#> <DenseNeuroVol> [806.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108.5, -46.25
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 1.000]
```

`read_vol` returns an object of class `NeuroVol` object which extends an
R `array` and has 3 dimensions (x,y,z).

``` r
    class(vol)
#> [1] "DenseNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
    
    is.array(vol)
#> [1] TRUE
    
    dim(vol)
#> [1] 64 64 25
    
    vol[1,1,1]
#> [1] 0
    
    vol[64,64,24]
#> [1] 0
```

Arithmetic can be performed on images as if they were ordinary `array`s:

``` r
    
    vol2 <- vol + vol
    sum(vol2) == 2 * sum(vol)
#> [1] TRUE
    
    vol3 <- vol2 - 2*vol
    all(vol3 == 0)
#> [1] TRUE
```

## Inspecting geometry and spatial metadata

Each `NeuroVol` has an associated `NeuroSpace` describing its geometry
(dimensions, spacing, origin, axes/orientation).

``` r
    sp <- space(vol)
    sp                 # human-readable summary
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108.5, -46.25
#>   Orientation   : LAS
#>   Voxels        : 102,400
    dim(vol)           # spatial dimensions (x, y, z)
#> [1] 64 64 25
    spacing(vol)       # voxel size in mm
#> [1] 3.5 3.5 3.7
    origin(vol)        # image origin
#> [1]  112.00 -108.50  -46.25
```

You can convert between indices, voxel grid coordinates, and real-world
coordinates:

``` r
    idx <- 1:5
    g   <- index_to_grid(vol, idx)     # 1D index -> (i,j,k)
    w   <- index_to_coord(vol, idx)    # 1D index -> world coords
    idx2 <- coord_to_index(vol, w)     # back to indices
    all.equal(idx, idx2)
#> [1] "Mean relative difference: 0.3333333"
```

A numeric image volume can be converted to a binary image as follows:

``` r
    
    vol2 <- as.logical(vol)
    class(vol2)
#> [1] "LogicalNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
    print(vol2[1,1,1])
#> [1] FALSE
```

## Masks and LogicalNeuroVol

Create a mask from a threshold or an explicit set of indices. Masks are
`LogicalNeuroVol` and align with the 3D space.

``` r
    # Threshold-based mask
    mask1 <- as.mask(vol > 0.5)
    mask1
#> <DenseNeuroVol> [406.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108.5, -46.25
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 1.000]

    # Index-based mask
    idx_hi <- which(vol > 0.8)
    mask2 <- as.mask(vol, idx_hi)
    sum(mask2) == length(idx_hi)
#> [1] TRUE

    # Use a mask to compute a summary
    mean_in_mask <- mean(vol[mask1@.Data])
    mean_in_mask
#> [1] 1
```

We can also create a `NeuroVol` instance from an `array` or `numeric`
vector. First we consruct a standard R `array`:

``` r
    x <- array(0, c(64,64,64))
```

Now we reate a `NeuroSpace` instance that describes the geometry of the
image including, at minimum, its dimensions and voxel spacing.

``` r
    bspace <- NeuroSpace(dim=c(64,64,64), spacing=c(1,1,1))
    vol <- NeuroVol(x, bspace)
    vol
#> <DenseNeuroVol> [2 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 64
#>   Spacing       : 1 x 1 x 1 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 0.000]
```

We do not usually have to create `NeuroSpace` objects, because geometric
information about an image is automatically determined from information
stored in the image file header. Thus, `NeuroSpace` objects are usually
copied from existing images using the `space` extractor function when
needed:

``` r
    vol2 <- NeuroVol((vol+1)*25, space(vol))
    max(vol2)
#> [1] 25
    space(vol2)
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 64
#>   Spacing       : 1 x 1 x 1 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#>   Voxels        : 262,144
```

## Slicing and quick visualization

The easiest way to view a volume is with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), which shows a
3 x 3 montage of evenly-spaced axial slices:

``` r
plot(vol)
```

![3x3 montage of axial
slices.](ImageVolumes_files/figure-html/plot_montage-1.png)

Default plot() montage

You can also extract a single 2D slice for display using standard array
indexing:

``` r
    z <- ceiling(dim(vol)[3] / 2)
    image(vol[,,z], main = paste("Slice z=", z))
```

![Mid-slice of example volume (grayscale
image).](ImageVolumes_files/figure-html/slice_mid-1.png)

Mid-slice of example volume

## Reorienting and resampling

You can change an image’s orientation and voxel spacing. Use
[`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)
to remap axes (e.g., to RAS) and
[`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
to match a target space.

``` r
    # Reorient the space (LPI -> RAS) and compare coordinate mappings
    sp_lpi <- space(vol)
    sp_ras <- reorient(sp_lpi, c("R","A","S"))
    g     <- t(matrix(c(10, 10, 10)))
    world_lpi <- grid_to_coord(sp_lpi, g)
    world_ras <- grid_to_coord(sp_ras, g)
    # world_lpi and world_ras differ due to axis remapping
```

Resample to a new spacing or match a target `NeuroSpace`:

``` r
    # Create a target space with 2x finer resolution
    sp  <- space(vol)
    sp2 <- NeuroSpace(sp@dim * c(2,2,2), sp@spacing/2, origin=sp@origin, trans=trans(vol))

    # Resample (trilinear)
    vol_resamp <- resample_to(vol, sp2, method = "linear")
    dim(vol_resamp)
```

## Downsampling

Reduce spatial resolution to speed up downstream operations.

``` r
    # Downsample by target spacing
    vol_ds1 <- downsample(vol, spacing = spacing(vol)[1:3] * 2)
    dim(vol_ds1)
#> [1] 32 32 32

    # Or by factor
    vol_ds2 <- downsample(vol, factor = 0.5)
    dim(vol_ds2)
#> [1] 32 32 32
```

## Writing a NIFTI formatted image volume

When we’re ready to write an image volume to disk, we use `write_vol`

``` r
    write_vol(vol2, "output.nii")
    
    ## adding a '.gz' extension results ina gzipped file.
    write_vol(vol2, "output.nii.gz")
```

You can also write to a temporary file during workflows:

``` r
    tmp <- tempfile(fileext = ".nii.gz")
    write_vol(vol2, tmp)
    file.exists(tmp)
#> [1] TRUE
    unlink(tmp)
```

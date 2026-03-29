# neuroim2 ![neuroim2 hex logo](https://raw.githubusercontent.com/bbuchsbaum/neuroim2/master/docs/hex-neuroim2.png)

Data structures and I/O for volumetric brain imaging with a focus on
fMRI. This is the successor to
[`neuroim`](https://github.com/bbuchsbaum/neuroim) with improved S4
classes, sparse/dense 3D–4D representations, and fast
resampling/filtering.

**Website:** <https://bbuchsbaum.github.io/neuroim2/> **CRAN:**
<https://cran.r-project.org/package=neuroim2>

## Key Features

- **NIfTI / AFNI support** — read and write `.nii`, `.nii.gz`, and AFNI
  `.BRIK`/`.HEAD` files
- **S4 class hierarchy** — `NeuroVol` (3D volumes) and `NeuroVec` (4D
  time-series) with dense, sparse, memory-mapped, and clustered variants
- **Coordinate systems** — full affine-aware transforms between voxel,
  grid, and world coordinates
- **Searchlights & ROIs** — spherical, cuboid, and random searchlight
  iterators plus flexible ROI construction
- **Spatial filtering** — Gaussian blur, guided filter, bilateral
  filter, and connected-component labelling
- **Resampling** — volume-to-volume resampling with configurable
  interpolation
- **Visualization** — slice montages, orthographic views, and overlay
  plotting via ggplot2

## Installation

### CRAN

``` r
install.packages("neuroim2")
```

### R-universe (daily builds)

``` r
install.packages("neuroim2",
  repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))
```

### Development version (GitHub)

``` r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/neuroim2")
```

## Usage

Read a NIFTI image and perform simple operations:

``` r
library(neuroim2)
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vol <- read_vol(fname)

v1 <- vol[1,1,1]
vol2 <- vol + vol
all(vol == (vol2 - vol))
```

``` R
## [1] TRUE
```

Create a 4D image from volumes:

``` r
vec <- vec_from_vols(list(vol, vol, vol2))
series1 <- vec[1,1,1,]
length(series1)
```

``` R
## [1] 3
```

## Vignettes

Browse the full set of
[articles](https://bbuchsbaum.github.io/neuroim2/articles/index.html) on
the pkgdown site:

| Getting Started                                                                              | Analysis Workflows                                                                             | Advanced                                                                                     |
|----------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
| [Overview](https://bbuchsbaum.github.io/neuroim2/articles/Overview.html)                     | [Analysis Workflows](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.html)    | [Image Volumes](https://bbuchsbaum.github.io/neuroim2/articles/ImageVolumes.html)            |
| [Choosing Backends](https://bbuchsbaum.github.io/neuroim2/articles/ChoosingBackends.html)    | [Slice Visualization](https://bbuchsbaum.github.io/neuroim2/articles/slice-visualization.html) | [NeuroVector](https://bbuchsbaum.github.io/neuroim2/articles/NeuroVector.html)               |
| [Coordinate Systems](https://bbuchsbaum.github.io/neuroim2/articles/coordinate-systems.html) | [Cookbook](https://bbuchsbaum.github.io/neuroim2/articles/Cookbook.html)                       | [Regions of Interest](https://bbuchsbaum.github.io/neuroim2/articles/regionOfInterest.html)  |
| [Volumes & Vectors](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.html)   | [Smoothing](https://bbuchsbaum.github.io/neuroim2/articles/Smoothing.html)                     | [Clustered NeuroVec](https://bbuchsbaum.github.io/neuroim2/articles/clustered-neurovec.html) |
| [Resampling](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.html)                 |                                                                                                | [Pipelines](https://bbuchsbaum.github.io/neuroim2/articles/pipelines.html)                   |

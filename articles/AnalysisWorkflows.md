# ROIs, Searchlights, and Pipelines

Most analysis work in `neuroim2` reduces to the same pattern: define a
spatial support, extract values from it, then summarize or map those
values into a new result. This article shows the three main ways to do
that.

## Extract a time series from one ROI

Start with a 4D image and define one spherical region of interest.

``` r

vec_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vol <- read_vol(vec_file)
vec <- read_vec(vec_file)
```

``` r

roi <- spherical_roi(vol, c(12, 12, 12), radius = 6)
roi_ts <- series_roi(vec, roi)

dim(values(roi_ts))
#> [1]  4 19
```

``` r

stopifnot(length(roi) > 0L)
stopifnot(nrow(values(roi_ts)) == dim(vec)[4])
```

That gives you one compact object containing the time series from every
voxel inside the ROI.

## Split a masked volume into parcels

When you already have a parcellation or cluster assignment,
[`split_clusters()`](https://bbuchsbaum.github.io/neuroim2/reference/split_clusters-methods.md)
turns one `NeuroVec` into a list of region-wise objects.

``` r

set.seed(1)

mask_vol <- vol > 0
cluster_ids <- sample(1:4, sum(mask_vol), replace = TRUE)
clustered <- ClusteredNeuroVol(mask_vol, cluster_ids)
parts <- split_clusters(vec, clustered)

length(parts)
#> [1] 4
```

``` r

part_means <- map_dbl(parts, ~ mean(values(.)))
part_means
#> 1 2 3 4 
#> 1 1 1 1
```

``` r

stopifnot(length(parts) == 4L)
stopifnot(all(is.finite(part_means)))
```

This is the right pattern when the support is fixed in advance by an
atlas, parcel map, or clustering step.

## Build searchlights lazily

Searchlights define many overlapping ROIs, one centered at each voxel in
a mask. The lazy form is useful because you only realize the
neighborhoods you actually touch.

``` r

sl <- searchlight(mask_vol, radius = 4, eager = FALSE, nonzero = FALSE)
first_sl <- sl[[1]]

nrow(coords(first_sl))
#> [1] 6
```

``` r

stopifnot(nrow(coords(first_sl)) > 0L)
```

The eager form is better when you want to iterate repeatedly over the
full set and can afford the up-front construction cost.

## Map a simple statistic over the first few searchlights

You do not need a specialized pipeline framework to use split-map-reduce
style workflows. A small list of ROIs plus `purrr::map_*()` is often
enough.

``` r

first_five <- lapply(seq_len(5), function(i) sl[[i]])
first_five_means <- map_dbl(first_five, ~ mean(values(.)))

first_five_means
#> [1] 0.5000000 0.8333333 0.8333333 0.5000000 0.5000000
```

``` r

stopifnot(length(first_five_means) == 5L)
stopifnot(all(is.finite(first_five_means)))
```

The pattern is the same whether the pieces come from parcels, ROIs,
connected components, or searchlights:

1.  define the spatial pieces
2.  extract the values you need
3.  reduce them into a scalar, vector, or model result

## Where to go next

- [`vignette("VolumesAndVectors")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md)
  for the core data containers
- [`vignette("regionOfInterest")`](https://bbuchsbaum.github.io/neuroim2/articles/regionOfInterest.md)
  for the full ROI API surface
- [`vignette("pipelines")`](https://bbuchsbaum.github.io/neuroim2/articles/pipelines.md)
  for older split-map-reduce examples
- [`vignette("clustered-neurovec")`](https://bbuchsbaum.github.io/neuroim2/articles/clustered-neurovec.md)
  for parcel-based data structures

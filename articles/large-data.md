# Large data

``` r

library(neuroim2)
```

Dense is the default because it is the simplest thing that works, not
because it is the right thing at scale. An R array holds `double`s
whatever the file said, so a 56 MB `int16` run is 236 MB the moment you
read it, and a session of runs will not fit. The two out-of-core routes
below are not a consolation prize for running out of memory: for the
scattered, voxel-at-a-time access that searchlight, ROI and connectivity
work actually perform, they are also *faster* than loading everything
first.

On an ordinary run (64 × 64 × 36 × 200), pulling 5,000 scattered voxel
time courses measured **44 ms** through `read_vec(path, mask = brain)`
and **144 ms** through `mode = "mmap"` including the open, against **4.8
s** for reading the whole image and then indexing it. So:

- **Pass a `mask`** whenever the analysis is confined to one — which is
  nearly always. It is the cheapest option on both axes.
- **Use `mode = "mmap"`** when the access order is scattered and unknown
  ahead of time.
- **Stay dense** when you genuinely need every voxel of a small image in
  memory at once.

What follows measures each backend so the choice is yours rather than a
recommendation you have to take on faith. The backends differ enormously
in *which* access pattern they make cheap, and picking the wrong one is
worse than any of the defaults.

Everything below measures the same dataset, written to an uncompressed
NIfTI so that every backend can read it:

``` r

mask <- demo_mask()
brain <- as.mask(mask)
bold <- demo_bold(n_time = 60)

path <- tempfile(fileext = ".nii")
write_vec(bold, path)

dim(bold)
#> [1] 64 64 25 60
round(file.size(path) / 1e6, 1)
#> [1] 24.6
```

The file is uncompressed. Every backend can read `.nii.gz` too, but two
of them pay for it: memory mapping silently decompresses to a temporary
uncompressed copy, and file-backed access constructs happily and then
fails at the first read. Write uncompressed if you intend to map.

## What each backend costs

``` r

dense <- read_vec(path)
sparse <- read_vec(path, mask = brain)
mapped <- read_vec(path, mode = "mmap")
filebacked <- read_vec(path, mode = "filebacked")

mb <- function(x) round(as.numeric(object.size(x)) / 1e6, 1)

data.frame(
  backend = c("DenseNeuroVec", "SparseNeuroVec", "MappedNeuroVec", "FileBackedNeuroVec"),
  memory_MB = c(mb(dense), mb(sparse), mb(mapped), mb(filebacked))
)
#>              backend memory_MB
#> 1      DenseNeuroVec      49.2
#> 2     SparseNeuroVec      15.1
#> 3     MappedNeuroVec       0.0
#> 4 FileBackedNeuroVec       0.0
```

They all hold the same data:

``` r

ref <- as.numeric(series(dense, 32, 32, 12))

vapply(
  list(sparse = sparse, mapped = mapped, filebacked = filebacked),
  function(x) isTRUE(all.equal(as.numeric(series(x, 32, 32, 12)), ref)),
  logical(1)
)
#>     sparse     mapped filebacked 
#>       TRUE       TRUE       TRUE
```

## The number that actually decides it

Memory is the easy axis. The one that catches people is access
*pattern*. Here is the same work — twenty whole volumes, then twenty
voxel time courses — on each backend:

``` r

per_call <- function(f, n) unname(system.time(for (i in seq_len(n)) f(i))["elapsed"] / n)

bench <- function(x, n = 20) {
  c(
    volumes = per_call(function(i) invisible(x[[i]]), n),
    series = per_call(function(i) series(x, 32, 32, 12), n)
  )
}

timings <- rbind(
  dense = bench(dense),
  sparse = bench(sparse),
  mapped = bench(mapped),
  filebacked = bench(filebacked, n = 5)
)

signif(timings, 2)
#>            volumes series
#> dense       0.0014  5e-05
#> sparse      0.0015  1e-04
#> mapped      0.0043  1e-04
#> filebacked  0.0320  1e-02
```

``` r

round(timings["filebacked", "series"] / timings[c("dense", "sparse", "mapped"), "series"])
#>  dense sparse mapped 
#>    200    100    100
```

Seconds per call, and the ratio is computed rather than quoted, so it
stays honest on your machine. Fetching one voxel’s time course from a
`FileBackedNeuroVec` costs three to four orders of magnitude more than
it costs anywhere else, because each query materialises the whole 4D
array as doubles before discarding all but one row — the disk read is a
small part of it. **Use it only when you genuinely work volume at a
time.** `MappedNeuroVec` is cheap on both patterns and is the better
default.

## The ladder

**Dense** is what you get with no arguments. Predictable arithmetic,
indexing and plotting, and nothing to think about — at the cost of
holding every voxel as a `double`. Right for a single small image, wrong
for a session.

**Sparse** stores only the voxels a mask admits. Here the brain is 29%
of the grid, and the object is about a third the size:

``` r

c(in_mask = sum(brain), total = prod(dim(brain)))
#> in_mask   total 
#>   29532  102400
c(dense_MB = mb(dense), sparse_MB = mb(sparse))
#>  dense_MB sparse_MB 
#>      49.2      15.1
```

Dimensions are unchanged — sparsity is storage, not shape. Reach for it
when a mask defines your analysis anyway.

**Memory-mapped** leaves the data on disk and lets the operating system
page it in. Pointed at a `.nii.gz` it decompresses to a temporary
uncompressed copy first, so you pay the disk either way — better to
store uncompressed deliberately than to have a full copy appear in
[`tempdir()`](https://rdrr.io/r/base/tempfile.html).

``` r

class(mapped)[1]
#> [1] "MappedNeuroVec"
dim(mapped)
#> [1] 64 64 25 60
```

[`as_mmap()`](https://bbuchsbaum.github.io/neuroim2/reference/as_mmap.md)
converts an in-memory object, writing an uncompressed file if it has to:

``` r

mmap_path <- tempfile(fileext = ".nii")
converted <- as_mmap(dense, file = mmap_path, overwrite = TRUE)

class(converted)[1]
#> [1] "MappedNeuroVec"
isTRUE(all.equal(as.numeric(series(converted, 32, 32, 12)), ref))
#> [1] TRUE
```

Give it an explicit `file` when you can. Left to itself it writes into
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) and you have a
second full copy of the data with no obvious owner.

**BigNeuroVec** backs a masked series with an on-disk `bigstatsr`
matrix, which is the option for data larger than RAM that still needs
voxelwise access:

``` r

big <- read_vec(path, mode = "bigvec", mask = brain)

class(big)[1]
#> [1] "BigNeuroVec"
round(mb(big), 2)
#> [1] 1
```

### What object.size() does not see

The two on-disk backends report near-zero above, and that is
[`object.size()`](https://rdrr.io/r/utils/object.size.html) telling you
the *R object* is a handle. It is not telling you the data is free:
mapped pages are charged to the process as you touch them, and
`BigNeuroVec` keeps a backing file — about 14 MB here — beside the
original. Measure resident memory and disk when the difference matters.

## Dropping to parcels

Every backend so far keeps one time course per voxel. The largest
reduction available comes from deciding you do not need them: a
`ClusteredNeuroVec` stores one series per *parcel*.

Parcels have to be spatially contiguous for any of this to mean
anything. Assigning labels at random would put all 100 parcels
everywhere at once and collapse their centroids onto the same point, so
cluster them by position:

``` r

set.seed(1)
coords_in <- index_to_grid(mask, which(as.vector(brain)))
km <- kmeans(coords_in, centers = 100, iter.max = 30)

parcels <- ClusteredNeuroVol(brain, km$cluster)
cvec <- ClusteredNeuroVec(bold, parcels)

num_clusters(cvec)
#> [1] 100
dim(cvec)
#> [1] 64 64 25 60
round(apply(centroids(cvec), 2, sd), 1)
#> [1]  9.3 10.8  6.2
```

Those centroid standard deviations are voxels: the parcels are spread
across the brain, which is what makes centroid distances meaningful
below.

It still behaves like a 4D image —
[`dim()`](https://rdrr.io/r/base/dim.html) is unchanged and indexing
works — but only 100 time courses exist behind it:

``` r

dim(as.matrix(cvec, by = "cluster"))
#> [1]  60 100

c(dense_MB = mb(dense), clustered_MB = mb(cvec))
#>     dense_MB clustered_MB 
#>         49.2          1.3
```

``` r

length(series(cvec, 32, 32, 12))
#> [1] 60
dim(cvec[, , , 1])
#> [1] 64 64 25

c(voxel_series = sum(brain), stored_series = num_clusters(cvec),
  fold = round(sum(brain) / num_clusters(cvec)))
#>  voxel_series stored_series          fold 
#>         29532           100           295
```

Every voxel in a parcel returns that parcel’s series: roughly a 300-fold
reduction in stored series here, and no within-parcel detail at all. For
atlas-based connectivity that is exactly right; for anything needing
voxelwise variation it is not.

One sharp edge: voxels outside the mask have no parcel, and
[`series()`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md)
returns `NA` for them rather than zero or an error. Restrict to the mask
before doing arithmetic over the whole grid.

Cluster-level searchlights come with it, using centroid distances rather
than voxel neighbourhoods:

``` r

windows <- cluster_searchlight_series(cvec, k = 3)

length(windows)
#> [1] 100
dim(values(windows[[1]]))
#> [1] 60  3
```

## Choosing

| Situation | Backend | Why |
|:---|:---|:---|
| It fits in memory | `DenseNeuroVec` | simplest, fastest, no caveats |
| A mask defines the analysis | `SparseNeuroVec` | stores only admitted voxels |
| Large uncompressed file, mixed access | `MappedNeuroVec` | cheap on both access patterns |
| Large file, strictly volume at a time | `FileBackedNeuroVec` | never for voxel series |
| Larger than RAM, needs voxel access | `BigNeuroVec` | on-disk matrix |
| Parcel-level analysis is enough | `ClusteredNeuroVec` | orders of magnitude fewer series |

If you are unsure: read with
[`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md),
and change only when a measurement tells you to.
[`object.size()`](https://rdrr.io/r/utils/object.size.html) and
[`system.time()`](https://rdrr.io/r/base/system.time.html) will get you
most of the way, as long as you remember that the first does not see
mapped pages or backing files.

## Where to go next

- [`vignette("volumes-and-vectors")`](https://bbuchsbaum.github.io/neuroim2/articles/volumes-and-vectors.md)
  — the containers these all implement
- [`vignette("regions-and-searchlights")`](https://bbuchsbaum.github.io/neuroim2/articles/regions-and-searchlights.md)
  — parcels and searchlights in depth
- [`?read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md),
  [`?as_mmap`](https://bbuchsbaum.github.io/neuroim2/reference/as_mmap.md),
  [`?ClusteredNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVec.md),
  [`?BigNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/BigNeuroVec-methods.md)

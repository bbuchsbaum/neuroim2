# Where neuroim2 stands against nilearn

Measured against `nilearn` 0.14.0 / scipy 1.17.1 / scikit-learn 1.8.0, with
neuroim2 built from source under R 4.3.3 with a real `RNiftyReg` 2.8.5. The
harness is in `dev/bench/nilearn/` and re-runs from scratch. Absolute times are
specific to this machine; the ratios are what transfer.

**Status.** All seven defects below have been fixed; see `NEWS.md` for what
changed and which changes alter results. Each section keeps the original finding
and adds the measurement after the fix, so the before/after is visible. The
seventh was found while closing the last two performance gaps, not by the
harness: `conn_comp()`'s local maxima were never pruned at all.

This is the companion to `dev/nibabel-gap-analysis.md`. nibabel is the
reference for file I/O and does no image processing; nilearn is the reference
for the processing neuroim2 also does.

---

## Scope

nilearn is a machine-learning library for neuroimaging. Most of it — decoders,
connectivity estimators, parcellation learners, GLMs, its plotting gallery — is
not what neuroim2 is for, and comparing there would produce a list of things
neuroim2 has chosen not to build. The overlap that *is* worth measuring is the
image-processing layer underneath: smoothing, resampling, reorientation,
connected components, masking, and getting signals out of ROIs, parcels and
searchlights.

Nine such areas were compared. Reference volume: 60x72x56 at 2 mm with an
83,563-voxel brain; run: the same grid, 60 timepoints.

**Summary as first measured: neuroim2 was correct almost everywhere it
overlaps, and faster than nilearn on the signal-extraction paths it was built
for, but had three real correctness defects and one very large performance
hole. All are now fixed, along with a fourth defect -- inert local maxima --
that only surfaced once the code around it was made fast enough to profile.**

---

## Correctness defects

### 1. `gaussian_blur()` does not deliver the smoothing it is asked for

`window` truncates the kernel at a fixed number of voxels, independent of
`sigma`. The default `window = 1` is a 3x3x3 kernel — one voxel each side —
so on 2 mm data a `sigma = 2 mm` request is truncated at 1 SD and renormalised.
Measured by smoothing an impulse and reading the FWHM back out:

| call | delivered FWHM | requested | shortfall |
|---|---:|---:|---:|
| `gaussian_blur(sigma = 2, window = 1)` | 3.49 mm | 4.71 mm | **26%** |
| `gaussian_blur(sigma = 3, window = 1)` | 3.70 mm | 7.06 mm | **48%** |
| `gaussian_blur(sigma = 4, window = 1)` | 3.76 mm | 9.42 mm | **60%** |
| `gaussian_blur(sigma = 2, window = 2)` | 4.53 mm | 4.71 mm | 4% |
| `gaussian_blur(sigma = 4, window = 4)` | 8.72 mm | 9.42 mm | 7% |
| `smooth_img(fwhm = 4.71)` | 4.71 mm | 4.71 mm | 0% |
| `smooth_img(fwhm = 9.42)` | 9.42 mm | 9.42 mm | 0% |

Three things make this worse than a documentation problem:

- **It is silent.** Nothing warns that the kernel has been clipped.
- **It is voxel-size dependent.** `sigma = 2, window = 1` gives 3.49 mm FWHM on
  2 mm voxels and 1.88 mm on 1 mm voxels — the same call smooths two
  acquisitions of the same brain differently.
- **The default is the worst case.** At `window = 1` the delivered smoothing
  saturates around 3.5–3.8 mm no matter how large `sigma` is.

scipy sizes the kernel from sigma (`truncate = 4.0`), so the requested
smoothing is what is applied. neuroim2 also had no FWHM interface, which is the
unit the field actually specifies smoothing in.

**Fixed.** `window` defaults to `NULL`, meaning "derive the half-width from
sigma and the voxel size", using scipy's own rounding. `fwhm` and `truncate` are
accepted. Requests are now honoured to better than 1% of FWHM on isotropic and
anisotropic grids alike, and matched-kernel timing is at parity (0.007 s against
nilearn's 0.008 s in 3-D). Passing `window` explicitly still truncates exactly
as before, so old analyses reproduce. This changes results for anyone who relied
on the old default.

### 2. `as_canonical()` silently discards data

`reorient(NeuroSpace, target)` rotates the affine but keeps the source `dim`.
When the reorientation permutes axes — which is the whole point — the output
grid no longer covers the data:

| | in | out | nonzero kept | pure permutation |
|---|---|---|---:|---|
| `reorder_img` | AIR 60x72x56 | RAS **56x60x72** | **100%** | yes |
| `as_canonical` | AIR 60x72x56 | RAS **60x72x56** | **76.4%** | no |

Nearly a quarter of the image is thrown away, with no warning. Two further
problems come from the same design decision — `as_canonical()` implements a
pure axis permutation by calling `resample()`, i.e. by going through an image
*registration* library:

- values are interpolated rather than permuted, so a reorientation is lossy
  where it should be exact;
- RNiftyReg refuses images smaller than 4 voxels on any axis, so
  `as_canonical()` errors on a 3x4x5 volume with *"Source image should have
  width at least 4 in all dimensions"*.

nibabel and nilearn do this with an integer axis permutation and flip: exact,
and 0.001 s.

**Fixed.** `as_canonical()` permutes and flips the array and rebuilds the affine
using the orientation primitives the package already exported, which reproduce
`nibabel::as_closest_canonical` voxel for voxel. 100% of the data is retained,
the output is a permutation of the input (checked with `identical()`), the
four-voxel minimum is gone, and it runs in 0.009 s instead of 0.376 s.
`reorient()` on a `NeuroSpace` permutes `dim` and `spacing` with the affine,
which was the bug underneath.

### 3. `searchlight()` and `searchlight_coords()` disagree about the centres

Same mask, same radius, on the 83,563-voxel brain:

```
searchlight(mask, radius = 6)          83,563 elements   (in-mask voxels)
searchlight_coords(mask, radius = 6)  241,920 elements   (every voxel in the grid)
```

Both are documented as *"visits every non-zero voxel in the mask as a potential
center voxel"*. `searchlight_coords()` defaults to `nonzero = FALSE`, which in
that function means "centre on every voxel", so the documented behaviour needs
`nonzero = TRUE`. Worse, `nonzero` is overloaded there: it selects the centre
set *and* filters the returned coordinates, where in `searchlight()` it only
does the latter. A user who writes the same arguments for both gets 2.9x the
work and most of the searchlights centred outside the brain.

**Fixed.** `searchlight_coords()` centres on nonzero mask voxels like every
other iterator in the package, and `nonzero` decides only what each searchlight
contains. Both now return 83,563 elements with identical contents. This changes
results for callers who did not pass `nonzero = TRUE`.

### 4. `NeuroSpace(dim, spacing, origin)` is a trap as a resample target

It builds a positive-diagonal affine. A typical NIfTI is LAS-to-RAS with a
negative x, so the "obvious" way to write a target grid mirrors the image:

```
resample(vol, NeuroSpace(dim, spacing = c(3,3,3), origin = origin(vol)))
  source mean 34.51  ->  resampled mean 0.0024
```

The target barely overlaps the source and the result is essentially empty, with
no warning. This is not a bug in `resample()` — given a faithful target affine
it is correct (below) — but a resample whose output has almost no overlap with
its input should say so.

---

### 5. `conn_comp()`'s local maxima were never pruned, and the rule was wrong

Not from the harness -- nilearn has no counterpart -- but from reading
`.pruneCoords()` while making it fast. Two defects, one hiding the other.

`.pruneCoords()` read the neighbour distances out of `dbscan::kNN()`'s result as
`ret$distances`. That component is called `dist`, and `$` does not partially
match a name *longer* than the one it is looking for, so the lookup returned
`NULL`. `NULL[, 2]` is `NULL`, `NULL < mindist & ...` is `logical(0)`, and
`any(logical(0))` is `FALSE` -- so the pruning loop exited on its first pass
every time and nothing was ever dropped. `local_maxima_dist` had no effect at
any value:

```
conn_comp(vol, threshold = 110, connect = "26-connect")
  12 components, 25,538 in-mask voxels
  local_maxima: 25,538 rows      <- every in-mask voxel, 25,527 from one component
```

Underneath that, the rule itself did not deliver the minimum distance the
argument documents. It compared each voxel only against its single nearest
neighbour, so with `local_maxima_dist = 15` three points at 0, 5 and 12 mm
valued 5, 1 and 9 keep the 5 *and* the 9 -- twelve millimetres apart -- because
the 5's nearest neighbour is the 1. "The nearest neighbour" was not even well
defined: on a voxel grid more than half of all points have several neighbours at
exactly the same distance, and which one came back depended on the kd-tree's
traversal order.

The rule is now the one the argument names: keep a voxel when no other voxel of
its component within `mindist` has a larger value. That is deterministic, needs
a single pass, and holds the minimum distance by construction.

```
                        local_maxima rows    min distance between maxima
26-connect   was              25,538                1 voxel (3 mm)
             now                  73                15.1 mm
6-connect    was              25,538                1 voxel (3 mm)
             now               5,039                15.1 mm
```

Pinned in `tests/testthat/test-conncomp-equivalence.R` against the rule stated
directly as an all-pairs scan, over randomised voxel lattices with ties in both
distance and value, plus the 0/5/12 mm case above.

---

## Correct, and competitive or better

### Signal extraction — neuroim2 wins

This is the layer neuroim2 was built for, and it shows. Values agree with
nilearn to float32 storage precision throughout.

| 60-timepoint run | nilearn | neuroim2 | |
|---|---:|---:|---|
| 200 spheres, r = 6 mm | 1.78 s | **0.064 s** | **28x faster** |
| 200 spheres, r = 8 mm | 1.71 s | **0.086 s** | **20x faster** |
| masked extraction, 83,563 voxels | 0.39 s | **0.071 s** | **5x faster** |
| 95 parcel means | 0.21 s | 0.50 s | 2.4x slower |
| max abs difference vs nilearn | — | 1e-4 (spheres), 4e-6 (parcels) | — |

### Searchlight neighbourhoods — identical, and faster

Against the `sklearn.neighbors` radius graph that `nilearn.decoding.SearchLight`
builds:

| r | centres | mean voxels | median | sklearn | neuroim2 |
|---|---:|---:|---:|---:|---:|
| 6 mm | 83,563 | 115.0 / 115.0 | 123 / 123 | 1.77 s | **1.06 s** |
| 8 mm | 83,563 | 235.8 / 235.8 | 257 / 257 | 2.71 s | **1.25 s** |

Neighbourhood statistics match exactly, neuroim2 is ~1.7x faster, and its
iterator is lazy where sklearn materialises the whole sparse graph.

### Connected components — exactly right

Component counts and largest-component sizes match `scipy.ndimage.label()` on
every connectivity pattern and at every mask density tested (5% to 80% of the
grid):

```
6-connect    4989 components, largest 1150    both
18-connect     87 components, largest 25423   both
26-connect     12 components, largest 25527   both
```

### Resampling — correct

Identity resample is bit-exact. Against nilearn at matched interpolation order:

```
2mm -> 3mm trilinear   cor 0.999993   RMS 0.18   (values ~0-200)
2mm -> 3mm cubic       cor 0.999194   RMS 1.94   (different spline basis)
nearest on 96 labels   96 -> 93 labels, same as nilearn
```

The cubic residual is RNiftyReg's cubic B-spline against scipy's cubic spline
interpolation — a basis difference, not an error. Timing is within 2.2x.

### Brain masking — better

`automask()` recovers the ground-truth mask exactly (Dice 1.0000) where
`compute_epi_mask` gets 0.9999. See the performance section for its cost.

### Utilities

Thresholding, binarising and arithmetic are 2-4x *faster* than nilearn's
equivalents. `mean()` on a 4-D image and `concat()` are slower (below).

---

## Performance: one very large hole (fixed)

### `conn_comp()` was 1000x slower than scipy, and super-linear

`conn_comp_3D()` (`R/conncomp.R`) was a pure-R two-pass union-find: a loop over
in-mask voxels, each iteration building a 26x3 neighbour matrix with `t()`,
indexing the label array with a coordinate matrix, and running an R-level
`find()`. `src/conncomp.cpp` existed but contained only a commented-out sketch
of that R code.

| volume | in-mask voxels | components | scipy | was | compiled labeller | + wrapper |
|---|---:|---:|---:|---:|---:|---:|
| 40x48x38 | 7,312 | 1,171 | 0.001 s | 1.4 s | 0.35 s | **0.024 s** |
| 91x109x91 | 90,357 | 12,396 | 0.016 s | 38.4 s | 2.10 s | **0.31 s** |
| 182x218x182 | 720,916 | 94,328 | 0.152 s | **916 s** | 16.1 s | **3.37 s** |
| *labelling alone*, 91x109x91 | 90,357 | 12,396 | 0.016 s | -- | 0.027 s | **0.027 s** |

Counts match exactly at every size, throughout. The compiled labeller is at
parity with scipy; what is left in `conn_comp()` is the per-cluster voxel lists,
the cluster table and the local maxima it returns and `scipy.ndimage.label()`
does not. Three other things came out of the same work: `ClusteredNeuroVol()`
filled its cluster map with a name-by-name list lookup that was quadratic in the
cluster count (2.2 s at 12,396 clusters); the wrapper re-derived a coordinate
matrix per component with `as.matrix()`, and then subset a data frame per
component to build the voxel lists (0.39 s at 4,989 components, against 0.03 s
for the labelling); and the local-maxima search called `dbscan::kNN()` once per
component per round, where the kd-tree build and dbscan's own argument checking
cost more than the search. That last one is now a single compiled pass over a
uniform grid -- see defect 5 above, which is what made it worth looking at.

The remaining ratio no longer grows with size: 18-22x across a 100x range in
voxel count, against 15-1000x before, and flat because both sides are now
compiled.

Output is unchanged except for the local maxima, which were wrong; the rest is
pinned against the previous implementation in
`tests/testthat/test-conncomp-equivalence.R`.

### `automask()` inherited it

| | voxels | Dice | time |
|---|---:|---:|---:|
| `compute_epi_mask` | 83,547 | 0.9999 | 0.031 s |
| `automask` (was) | 83,563 | **1.0000** | 14.0 s (448x) |
| `automask` (now) | 83,563 | **1.0000** | **0.39 s** (12x) |

92% of its runtime was in the labeller, which it calls up to four times.

### Smaller gaps

| operation | was | now | nilearn |
|---|---:|---:|---:|
| `as.matrix()` on a 4-D image | 0.180 s | **~0 s** | -- |
| `mean()` on a 4-D image | 0.197 s | **0.126 s** | 0.012 s |
| `resample()` 2mm -> 3mm | 0.066 s | **0.018 s** | 0.006 s |
| Gaussian smooth, matched kernel, 3-D | 0.009 s | **0.007 s** | 0.007 s |
| Gaussian smooth, matched kernel, 4-D, 60 volumes | 1.66 s | **0.20 s** | 0.51-0.71 s |

`as.matrix()` on a `DenseNeuroVec` duplicated the payload to reshape it to
voxels-by-time; a reshape is a change of `dim`, and R copies on the first write
anyway. `mean()` now goes through `.rowMeans()` with no reshape at all -- the
residual 10x is structural, since an R image holds `double` where numpy reads
the file's `float32`. `resample()` was scanning the image to infer an output
datatype for headers whose geometry is all it uses.

## What is left

The seven defects are fixed and both remaining performance gaps are closed.
What the harness still shows, in the order it would be worth looking at:

**1. `mean()` and other 4-D reductions are ~10x off numpy**, for a structural
reason: an R image holds `double` where numpy reads the file's `float32`, so
every reduction moves twice the bytes. Not fixable without changing what a
`NeuroVec` stores.

**2. `conn_comp()` is 14-37x off `scipy.ndimage.label()`** on the whole call,
but the labelling alone is 1.6-2.0x, and that is the like-for-like number --
`conn_comp()` also returns per-cluster voxel lists, a cluster table and the
local maxima. Turning those off with `cluster_table = FALSE, local_maxima =
FALSE` brings the 6-connect case from 0.092 s to 0.044 s. Whether they should
stay on by default is worth revisiting; the R-level assembly of the voxel lists
is what is left in the middle.

**3. `automask()` is 8x off `compute_epi_mask`** (0.26 s against 0.031 s), down
from 448x. It calls the labeller up to four times and is now dominated by the
rest of its own R code rather than by labelling.

**4. `concat()` is ~4x off `concat_imgs`**, and was not investigated.

**5. `bilateral_filter`, `guided_filter`, `enhance_stat_map` and the `cgb_*`
graph-smoothing functions** were never profiled or probed -- none has a close
nilearn counterpart, so the harness says nothing about them.

### Where neuroim2 is now ahead

| operation | neuroim2 | nilearn | |
|---|---:|---:|---|
| Gaussian smooth, 4-D, 60 volumes | 0.20 s | 0.51-0.71 s | 2.6-3.6x faster |
| Spherical ROI signal, r = 6 mm | 0.044 s | 1.78 s | 40x faster |
| Spherical ROI signal, r = 8 mm | 0.063 s | 1.71 s | 27x faster |
| Masked extraction | 0.030 s | 0.388 s | 13x faster |
| Searchlight neighbourhoods, r = 8 mm | 1.18 s | 2.71 s | 2.3x faster |
| `concat_imgs` equivalent | 0.118 s | 0.243 s | 2.1x faster |
| Threshold / binarise | 0.002-0.003 s | 0.003-0.004 s | 1.3x faster |

## Reproducing

```sh
cd dev/bench/nilearn
export NEUROIM2_NL_DIR=/tmp/nl
python3 make_data.py && python3 ref_py.py && Rscript probe_r.R && python3 compare.py
# add NEUROIM2_NL_SCALING=1 for the connected-component scaling probe, which
# now takes seconds rather than the quarter of an hour it first did
```

Caveats. One machine, medians of 1–3 repetitions, with `gc()` outside the timed
region so that one repetition does not pay for the previous one's garbage --
without it, smoothing a 60-volume run alternates 0.19 s and 0.86 s depending on
which call trips the collector. `conn_comp()` also returns cluster sizes,
per-cluster voxel lists, a cluster table and local maxima where
`scipy.ndimage.label()` returns only the labelling, which is most of what is
left of that gap and none of what it was. The smoothing timing is reported twice because the default
neuroim2 kernel is narrower than nilearn's, so the naive comparison flatters
neuroim2. Nothing here examines `bilateral_filter`, `guided_filter`,
`enhance_stat_map`, the `cgb_*` graph-smoothing functions or the plotting
layer, none of which have a close nilearn counterpart to compare against.

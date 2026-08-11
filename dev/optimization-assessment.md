# neuroim2 optimization assessment

Ranked by impact on the paths users actually run. Every number below was measured,
not estimated; see **Method** at the end for how, and for what the measurements do
and do not cover.

Reference workload throughout: a 91×109×91 volume at 2 mm with a 290,137-voxel brain
mask, radius-8 mm searchlights (257-voxel spheres), 250 timepoints. That is an
ordinary whole-brain MVPA setup.

## Summary

| # | Area | Now | Achievable | Gain |
|---|------|-----|-----------|------|
| 1 | S4 construction in ROI loops | 454 µs/object | 41 µs | **11×**, bit-identical objects |
| 2 | `spherical_roi()` per-centre work | 660 µs/ROI | 65 µs | **10×** combined with #1 |
| 3 | Searchlight enumeration, whole brain | 191 s | 19 s (0.4 s index-only) | **10×** / **495×** |
| 4 | `series()` on a dense vec | 0.47–0.71 ms/ROI | 0.05 ms | **9–14×** |
| 5 | `series()` on a sparse vec | 645 µs/ROI | 342 µs | **1.9×** (memory-bound) |
| 6a | `gaussian_blur`, `normalize = TRUE`, window 2–3 | 0.07–0.16 s | 0.054–0.055 s | **1.3–2.9×**, exact |
| 6b | `gaussian_blur`, `normalize = FALSE`, window 2–4 | 0.15–0.78 s | 0.023–0.051 s | **6.6–15.2×**, exact |
| 7 | Scalar R loops (clustered/hyper/kernel) | — | — | **16×+** on the ones measured |

The headline: **enumerating a whole-brain searchlight currently costs ~191 seconds
before any analysis runs.** Roughly 70% of that is S4 bookkeeping, not geometry.

---

## 1. S4 object construction dominates every ROI loop

`new("ROIVolWindow", ...)` costs **454 µs** and the cost is essentially independent
of ROI size:

```
new() with     1 voxels     437 us
new() with    33 voxels     460 us
new() with   257 voxels     497 us
new() with  2000 voxels     467 us
```

That is fixed overhead, and it is 69% of `spherical_roi()`'s 660 µs.

**Root cause.** Decomposing the class hierarchy by adding one feature at a time:

```
coords slot only, no .Data, no space, no validity :   37 us
+ numeric .Data part                              :  143 us
+ NeuroSpace slot                                 :  156 us
+ validity function                               :  168 us
+ one more inheritance level (= ROIVolWindow)     :  474 us   <-- 2.8x jump
```

The jump is the second `contains=` level. `ROIVolWindow` extends `ROIVol`
(`R/all_class.R:1550`), which itself extends `c("ROICoords", "numeric")`
(`R/all_class.R:1516`). R's default `initialize` for a class that inherits a basic
type through a multi-level chain walks that chain with `callNextMethod` and runs
`validObject` at each level — `validObject` alone measures 167 µs on this object,
because it recurses into the `NeuroSpace` slot and its nested `AxisSet3D`/`NamedAxis`
tree.

**Fix.** An internal constructor that clones a cached prototype and assigns slots:

```r
.roi_window_proto <- new("ROIVolWindow", numeric(0), space = <sp>,
                         coords = matrix(0, 0, 3),
                         center_index = NA_integer_, parent_index = NA_integer_)

.new_roi_window <- function(vals, space, coords, center_index, parent_index) {
  o <- .roi_window_proto
  o@.Data <- vals; o@space <- space; o@coords <- coords
  o@center_index <- center_index; o@parent_index <- parent_index
  o
}
```

**454 µs → 41 µs (11×)**, and `identical(new(...), .new_roi_window(...))` is `TRUE` —
same class, same slots, same S4 bit. Validity is not skipped so much as not
*re-run per object*: the invariants (`ncol(coords) == 3`, `length(.Data) ==
nrow(coords)`) are cheap to assert directly in the one constructor that builds these,
and every caller in the package already satisfies them by construction.

This applies wherever ROI objects are minted in a loop — `spherical_roi()`,
`spherical_roi_set()`, `random_searchlight()`, `resampled_searchlight()`,
`clustered_searchlight()`.

*Longer-term alternative:* flattening `ROIVolWindow` to inherit `numeric` directly
(one `contains=` level, keeping `coords`/`space` as plain slots) would remove the
2.8× jump at the source. That is an API-visible change to the class graph, so it is a
separate decision from the constructor fix, which is invisible.

## 2. `spherical_roi()` redoes work that is constant across centres

`R/roi.R:657`. Per-ROI cost breakdown at radius 8 mm:

```
local_sphere()             15.0 us
order() sort               27.5 us
bvol[grid] values           6.0 us
rowSums centre match       16.0 us
gridToIndex3DCpp            3.5 us
new("ROIVolWindow")       454.0 us
```

Four separate problems:

- **The sphere is recomputed for every centre.** `make_spherical_grid()` →
  `local_sphere()` (`src/indexFuns.cpp:275`) evaluates `sqrt(pow(...) + pow(...) +
  pow(...))` over the whole candidate cube on every call. The offset set is identical
  for all centres — only the boundary clipping differs. Compute `(dx, dy, dz)` once,
  then per centre it is an add-and-clip.

- **`order()` at `R/roi.R:698` is a no-op.** `local_sphere` emits `i` outer, `j`
  middle, `k` inner, so rows already arrive in lexicographic `(x, y, z)` order.
  Verified across 300 random centres: `order(...)` is always `seq_len(nrow)`. It costs
  27.5 µs/ROI — the second most expensive line in the function — and buys nothing.

- **The centre match at `R/roi.R:725`** materializes an `n × 3` comparison matrix via
  `rowSums(grid == matrix(centroid, nrow(grid), 3, byrow = TRUE))`. With a precomputed
  offset template the centre's row is known before any clipping, or is a three-way
  `&` comparison at worst. 16 µs → ~1 µs.

- **Anisotropic voxels scan a cube 3× too large.** `src/indexFuns.cpp:277` uses a
  single `window = ceil(radius / min(spacing))` for all three axes:

  ```
  spacing 2/2/2      r=8mm:  729 candidates scanned,  729 needed, 257 kept  [1.0x]
  spacing 1/1/4      r=8mm: 4913 candidates scanned, 1445 needed, 489 kept  [3.4x]
  spacing 0.8/0.8/3  r=8mm: 9261 candidates scanned, 3087 needed, 1137 kept [3.0x]
  ```

  Per-axis windows fix this. Anisotropic acquisition is common enough that this is
  not a corner case.

## 3. What the searchlight actually costs

Combining #1 and #2, measured end to end, with output objects verified
`identical()` to what the current code returns:

```
current spherical_roi()        660.0 us/ROI   191.5 s whole brain
offset template + prototype     65.3 us/ROI    19.0 s   [10.1x]  bit-identical
index-only C++ iterator          1.3 us/ROI     0.4 s   [495x]
```

The third row is the interesting one. Most searchlight consumers do
`series(vec, coords(sl))` and never touch the ROI object as an object — they want the
voxel set. A C++ iterator that emits linear indices (plus the centre's row) for all
centres in one pass makes enumeration free: 0.4 s instead of 191 s. `searchlight_coords()`
already has the right shape for this; it just calls `spherical_roi()` underneath
(`R/searchlight.R:595`) and pays the full object-construction price to then throw the
object away.

Recommended shape: a fused C++ core that returns index sets, with the existing
`ROIVolWindow`-returning API layered on top for callers that need it (now using the
cheap constructor). RcppParallel is already in `LinkingTo` and the loop over centres
is embarrassingly parallel.

**Related:** `spherical_roi_set()` (`R/roi.R:1108`) documents itself as "more efficient
than calling `spherical_roi` multiple times", but its body is a `for` loop calling
`spherical_roi` once per centroid. It is exactly as fast, and `searchlight(eager = TRUE)`
routes through it (`R/searchlight.R:687`) expecting a batched path that does not exist.
This is the natural place to put the fused implementation, at which point the docstring
becomes true.

## 4. `series()` on dense vectors

Two methods, both improvable:

```
series(NeuroVec, matrix)         0.707 ms   generic path
series(DenseNeuroVec, matrix)    0.467 ms   per-timepoint R loop
C++ gather                       0.050 ms   [14.1x / 9.3x]
```

- `R/neurovec.R:670` — the generic method builds an `(n·T) × 4` index matrix
  (`i[rep(1:nrow(i), each = d4), ]` then `cbind(..., 1:d4)`) and does one big matrix
  index. For a 257-voxel ROI over 250 timepoints that is a 64,250×4 allocation per
  call.
- `R/neurovec.R:785` — the dense method is better but loops over timepoints, and
  `lin + (t - 1L) * nels` promotes to double on every iteration because `nels` is
  `prod()`'s double. So it allocates a fresh double index vector per timepoint and
  indexes with doubles.

A single-pass C++ gather over `(voxel, timepoint)` handles both and is ~10× faster.
This is the most-called function in any ROI or searchlight analysis.

## 5. `series()` on sparse vectors

```
current  matrix(0) + out[,nz]<-    645 us/ROI   187 s whole brain
all-in-mask shortcut               436 us/ROI   127 s   [1.5x]
single-pass C++ gather             342 us/ROI    99 s   [1.9x]
```

`R/sparse_neurovec.R:265` allocates a zero matrix and then scatters into it with
`out[, nz] <- ...`. Inside a brain mask, 94% of searchlight voxels are in-mask on
average, so the zero-fill is almost entirely overwritten and the scatter goes through
R's general matrix subassign path.

Be honest about the ceiling here: 257 voxels × 250 timepoints × 8 bytes is 514 KB
gathered from scattered columns per ROI. This path is close to memory-bandwidth-bound,
so **~2× is the realistic gain, not 10×**. Still worth it — 187 s → 99 s on a
whole-brain sweep — but it should not be the first thing fixed.

## 6. `gaussian_blur` is not exploiting separability

`src/indexFuns.cpp` applies the full `(2w+1)³` kernel. On the unnormalized dense
path, a Gaussian is a tensor product, so three 1-D passes give the same answer in
`3(2w+1)` taps:

```
window=1  taps  27 vs  9   dense 0.049s  separable 0.026s  [ 1.9x]  max|diff| 6.7e-16
window=2  taps 125 vs 15   dense 0.152s  separable 0.023s  [ 6.6x]  max|diff| 1.0e-15
window=3  taps 343 vs 21   dense 0.381s  separable 0.025s  [15.2x]  max|diff| 2.4e-15
window=4  taps 729 vs 27   dense 0.777s  separable 0.051s  [15.2x]  max|diff| 2.1e-15
```

The mask-insulated default (`normalize = TRUE`) is *also* separable, as a ratio of two
separable convolutions — `blur(y) / blur(m)` — since the numerator and denominator
the current code accumulates are exactly those two convolutions. Here `y` must be
constructed as a zero-filled volume with `x` assigned only at in-mask voxels; it must
not be formed as `x * m`, because IEEE arithmetic makes `0 * NaN` and `0 * Inf`
non-finite and would leak those values into nearby in-mask outputs:

```
normalize=TRUE window=1: dense 0.020s  separable-ratio 0.048s  [0.4x]  max|diff| 1.1e-15
normalize=TRUE window=2: dense 0.070s  separable-ratio 0.054s  [1.3x]  max|diff| 1.0e-15
normalize=TRUE window=3: dense 0.158s  separable-ratio 0.055s  [2.9x]  max|diff| 1.5e-15
```

Agreement is at machine epsilon in both cases. The masked variant loses at `window = 1`
because the separable passes sweep the whole volume while the current code visits only
mask voxels — so gate on `window`: keep the current kernel at `window = 1`, use
separable above it. The default is `window = 1`, but anyone smoothing at a realistic
FWHM raises it, and that is where the 3–15× lives.

## 7. Scalar R loops over voxels

Four sites doing per-voxel work in interpreted R:

- **`R/clustered_neurovec.R:278`** — `series(ClusteredNeuroVec, i, j, k)` runs a nested
  loop over timepoints × voxels to copy values out of `x@ts`. The whole thing is one
  vectorized subset. Measured at 3,000 voxels × 250 timepoints: **0.124 s → 0.008 s
  (16×)**, output `identical()`.
- **`R/neurovol.R:906`** — `mapf(NeuroVol, Kernel)` originally convolved by looping
  over centres and doing `sum(x[cbind(ii, jj, kk)] * wts)` per voxel. The follow-up
  implementation uses each genuinely 3-D kernel tap as a fixed linear-index offset
  and vectorizes over centres. The old `kernel_filt_3d_cpp` was not usable here: it
  accepted two `NumericMatrix` objects and filtered only rows and columns, so wiring
  it into `mapf()` would have dropped the z dimension and changed boundary behavior.
- **`R/neurohypervec.R:448`** — `[` loops over voxels doing an `aperm` per voxel.
- **`R/neurohypervec.R:289`** — `as.dense` loops over features × trials rebuilding a
  volume each time.

## 8. Hygiene

- **Dead compiled code.** `kernel_filt_3d_cpp`, `radius_search_3d_nonisotropic`,
  `radius_search_3d_direct`, `radius_search_3d_precomputed`, and `local_spheres` were
  called from nowhere in `R/` or `tests/`. The follow-up implementation drops all five;
  `local_sphere` (singular) remains as the tested reference implementation.
- **`random_searchlight()`** (`R/searchlight.R:81`) rebuilds `remain_indices` with an
  O(n) filter on every iteration. With ~1,100 iterations over a 290k-voxel mask that is
  a few hundred million operations. A swap-with-last free list makes it O(1) amortized.
  Secondary to the above, but nearly free to fix.

---

## Suggested sequencing

1. **Cheap ROI constructor** (#1). Smallest diff, largest single win, bit-identical
   output, no API change. Roughly 3× on every searchlight path by itself.
2. **Offset template + drop the redundant `order()`/`rowSums`** (#2). Gets the combined
   gain to 10×, still bit-identical.
3. **C++ `series()` gather** (#4, #5). Touches the single most-called function.
4. **Fused index-only searchlight core** (#3) exposed through `spherical_roi_set()` /
   `searchlight_coords()`. This is where the 495× lives, and it is the right home for
   RcppParallel.
5. **Separable Gaussian** (#6) and the scalar loops (#7). Independent of the rest.

Steps 1–3 are the ones I would do first: they are localized, verifiable against the
existing output byte-for-byte, and together take a whole-brain searchlight from ~191 s
of enumeration plus ~187 s of extraction down to roughly 19 s plus 99 s.

## Method, and its limits

CRAN is unreachable from this container, so `RNifti`, `bigstatsr`, `RNiftyReg`, and
`deflist` could not be installed and **the package itself was never built or profiled**.
Instead I measured a standalone harness under R 4.3.3 + Rcpp:

- The C++ kernels (`local_sphere`, `indexToGridCpp`, `gridToIndex3DCpp`,
  `gaussian_weights`, `gaussian_blur_cpp`) are **verbatim copies** from `src/`.
- The S4 hierarchy (`NamedAxis` → `AxisSet3D` → `NeuroSpace`, `ROI` → `ROICoords` →
  `ROIVol` → `ROIVolWindow`) is **transcribed from `R/all_class.R`**, including both
  validity functions.
- The R-level hot paths (`spherical_roi` body, both `series` methods, the
  `ClusteredNeuroVec` loop) are **transcribed from the package sources**.
- Timings are medians of 5 repetitions after a warm-up call; proposed replacements are
  checked with `identical()` or `all.equal()` against the current implementation on
  every reported comparison.

What this means: the **ratios** should transfer, since they come from the same code
doing the same work. The **absolute** times are specific to this container and will
differ on other hardware. Two things the harness cannot see are S4 generic *dispatch*
overhead in the real package (measured separately at ~1.2 µs per `standardGeneric`
call, so a few percent at these ROI sizes, not a driver) and anything I/O-bound in the
NIfTI read path, which I did not examine.

Recommended next step before committing to the work: build the package normally and
confirm #1 and #2 with `profvis` on a real `searchlight()` call, which should take
about ten minutes once dependencies are available.

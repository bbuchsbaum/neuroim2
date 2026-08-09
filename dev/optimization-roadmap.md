# neuroim2 optimization roadmap

Companion to `dev/optimization-assessment.md`. That document says *where the time
goes*; this one says *what to do about it, in what order, and how each step gets
verified*. Phase 1 is done and on this branch; phases 2–4 are specified but not
implemented.

Every phase below names its verification up front, because the whole approach
here depends on being able to prove a change is invisible before claiming it is
fast.

---

## The verification harness (reusable, applies to every phase)

Three layers, all already in place:

1. **Golden sweep** — `dev/bench/golden.R` (scratch copy; see "Reproducing"
   below) captures `spherical_roi`, `spherical_roi_set`, `series`, and the four
   searchlight iterators across 8,613 configurations: 5 voxel spacings
   (isotropic and anisotropic), 3 volume shapes, 6 radii including exact
   boundary cases, centroids at every corner/edge/interior, `nonzero`
   TRUE/FALSE, `use_cpp` TRUE/FALSE, `NeuroSpace` and `NeuroVol` targets, fill
   variants, and the documented error paths. Capture before, compare after,
   require zero mismatches.
2. **Full test suite vs. a locked baseline** — 2,226 passing, 0 failing, 25
   errors (all attributable to offline dependency shims, see "Environment").
   Any new entry in the failing set is a regression.
3. **Permanent regression tests** — `tests/testthat/test-roi-series-fastpaths.R`
   (624 assertions), which pin equivalence rather than performance: the cheap
   constructor against `new()`, the offset template against `local_sphere()`,
   the compiled gathers against the R loops they replaced, the accessor-hook
   contract for lazy sparse subclasses, and a scaling guard that fails if a
   whole-volume copy is ever reintroduced into `series()`.

**Adversarial review is part of the process, not a courtesy.** Each phase gets
independent reviewers briefed to refute the equivalence claim, with distinct
attack surfaces (object/S4 semantics, geometry and floating-point boundaries,
C-level memory safety). Phase 1's review caught two things the golden sweep did
not, which is exactly the point:

- A **catastrophic performance inversion**: passing `x@.Data` into `.Call`
  duplicates the entire volume (165 ms per `series()` call on a 48×56×40×200
  image). The "optimization" was 100× *slower* than the code it replaced while
  producing bit-identical output — invisible to every correctness check.
- The **lazy-subclass accessor contract**: gating the sparse fast path on
  "is `@data` a double matrix" broke subclasses that keep a placeholder in
  `@data` and serve real values from an overridden `matricized_access()`. The
  gate is now on the concrete class.

Both are the same lesson: correctness tests do not catch performance bugs, and
type checks do not catch contract violations. Budget for review accordingly.

---

## Phase 1 — done (this branch)

| Change | File | Result |
|---|---|---|
| Prototype-cloning ROI constructor | `R/roi.R` | `new()` 454 µs → 41 µs, `identical()` output |
| Cached sphere-offset template, per-axis windows | `src/sphere_offsets.cpp`, `R/roi.R` | see below |
| Dropped the no-op `order()` on the compiled path | `R/roi.R` | ~27 µs/ROI |
| Three-column centre match instead of `rowSums` | `R/roi.R` | ~15 µs/ROI |
| Compiled `series()` gathers | `src/series_gather.cpp`, `R/neurovec.R`, `R/sparse_neurovec.R` | see below |

Measured end to end on the reference workload (91×109×91 @ 2 mm, 290,137-voxel
mask, radius 8 mm, 250 timepoints):

```
spherical_roi()                  660 us/ROI  ->  138 us/ROI      4.8x
series(DenseNeuroVec, matrix)   1.51 ms      -> 0.09 ms         16.8x
series(SparseNeuroVec, idx)      645 us      ->  342 us          1.9x
```

`spherical_roi` came in at 4.8× rather than the 10× the standalone harness
predicted. The harness modelled the arithmetic but not the S4 surface around it:
`space()`, `dim()`, `spacing()`, `grid_to_index()` dispatch and `[` on a
`LogicalNeuroVol` are now the majority of what remains. That residue is what
Phase 2 removes by not going through per-ROI R code at all.

---

## Phase 2 — done: thin the searchlight iterators

**What was actually built, and why it differs from the plan.** The plan projected
495x from emitting index sets for all centres in one pass. Measured first: a
whole-brain radius-8 mm searchlight over a 290k-voxel mask is ~0.9 GB of
coordinates, 1.5 GB as ROI objects. That is exactly why these iterators are lazy,
so bulk emission would have traded a memory hazard for a benchmark. The work went
instead into making the per-element call cheap:

- `src/searchlight_core.cpp` does the whole per-centre job in one pass: translate
  the cached template, clip, apply the mask, gather values, find the centre row.
- `.searchlight_plan()` hoists every invariant out of the iterator closures,
  which are built with `local()` so they retain only the plan and the centre grid.
- `spherical_roi_set()` expands all centres in one batched call;
  `searchlight(eager = TRUE)` inherits it.

```
searchlight_coords()        97.0 -> 7.7 us/elem    28.1 s -> 2.2 s   12.6x
searchlight(eager = FALSE)  92.3 -> 28.0 us/elem   26.8 s -> 8.1 s    3.3x
spherical_roi_set()                 31.5 us/ROI             9.1 s
```

ROI construction dropped again on the way: an S4 object extending a basic type
*is* that vector with the slots as attributes, so setting the attributes and
calling `asS4()` gives an `identical()` object. The mechanism costs ~3 us; with
the argument checks the constructor is ~11 us against `new()`'s 450-790 us.

**What the review found.** Four defects, all fixed:

- A **segfault** from `R_xlen_t` overflow in the length guards -- the same defect
  a previous review found in `series_gather.cpp`, reintroduced. Guards now use
  division-based checks.
- **`nonzero` disagreed about NA.** The compiled filter dropped missing voxels;
  the R filter (`vals != 0`) produced an `NA` logical subscript and emitted
  `NA NA NA` coordinate rows. Eager and lazy searchlights differed on 21 of 25
  windows on a volume with one `NA`. The rule now exists once, in C++.
- The lazy iterators had **stopped validating `radius < min(spacing)`** and
  silently returned one-voxel windows.
- `parent_index` was narrowed with a C cast past 2^31 voxels; now `NA`.

Also fixed while in the area: a pre-existing **SIGFPE** in `index_to_grid()` for
spaces whose running product of dimensions exceeds `int`.

**Still open, low value:** the `.sphere_offset_cache` is unbounded in entry *size*
(templates scale as radius cubed, and a pathological radius parks hundreds of MB
in the namespace). Bounding by total size rather than entry count would fix it.

---

## Phase 3 — separable Gaussian blur

**Why.** `gaussian_blur_cpp` applies the full `(2w+1)³` kernel. A Gaussian is a
tensor product, so three 1-D passes give the same answer in `3(2w+1)` taps.
Measured, agreeing to machine epsilon (max |diff| ~1e-15):

```
normalize = FALSE:  w=1  1.9x   w=2  6.6x   w=3 15.2x   w=4 15.2x
normalize = TRUE:   w=1  0.4x   w=2  1.3x   w=3  2.9x
```

The masked default (`normalize = TRUE`) is also separable, as the ratio of two
separable convolutions — `blur(x·m) / blur(m)` — because the numerator and
denominator the current code accumulates are exactly those two convolutions.

**What to build.** A separable kernel in `src/`, with `gaussian_blur()` choosing
by `window`: keep the current kernel at `window = 1` (where the separable version
loses, because it sweeps the whole volume while the current code visits only mask
voxels), switch above it. Gate on `window`, not on mask fraction — a simple,
predictable rule.

**Verification.** Golden capture of `gaussian_blur` over a grid of
(sigma, window, normalize, mask) and assert agreement to `1e-12` absolute. Add
the masked-NaN case from the docs: out-of-mask `NaN` must not leak in-mask.

**Risk.** Low, and the payoff is bounded by how often users raise `window` above
the default of 1. Worth doing, but it is not the searchlight.

---

## Phase 4 — scalar R loops and hygiene

Independent of everything above; each is small and self-contained.

| Item | File | Fix | Measured / expected |
|---|---|---|---|
| `series(ClusteredNeuroVec, i, j, k)` | `R/clustered_neurovec.R:278` | replace the nested timepoint × voxel loop with one vectorised subset (`ts[m, cids, drop=FALSE]`, NA for `cids == 0`) | **16×** measured, `identical()` output |
| `mapf(NeuroVol, Kernel)` | `R/neurovol.R:906` | per-voxel `sum(x[cbind(ii,jj,kk)] * wts)` in R; `kernel_filt_3d_cpp` already exists in `src/kernel_filter.cpp` and is never called — wire it up | large, unmeasured |
| `NeuroHyperVec` `[` | `R/neurohypervec.R:448` | per-voxel loop with an `aperm` per voxel; gather in one pass | large, unmeasured |
| `as.dense(NeuroHyperVec)` | `R/neurohypervec.R:289` | loops features × trials rebuilding a volume each time | moderate |
| `random_searchlight()` free list | `R/searchlight.R` | `remain_indices` is rebuilt with an O(n) filter every iteration; a swap-with-last free list makes it O(1) amortised | moderate |
| Dead compiled code | `src/radius_search_3d.cpp`, `src/kernel_filter.cpp` | `radius_search_3d_{nonisotropic,direct,precomputed}` and `local_spheres` are compiled and exported but called from nowhere in `R/` or `tests/` (~17 KB). Drop them, or wire `kernel_filt_3d_cpp` into `mapf()` | build time only |

---

## Resolved: the `use_cpp = FALSE` fallback (now deprecated)

Recorded here because the diagnosis is worth keeping. The legacy branch was
wrong in two independent ways:

1. **Off by +1 on every axis.** `make_spherical_grid()` documents a 0-based
   contract; the `dbscan` branch returned `cube`, built from
   `seq(centroid - w, centroid + w)` and therefore already 1-based, and
   `spherical_roi()` added 1 regardless. This also made `spherical_roi()` throw
   `subscript out of bounds` for centroids near the volume edge, because the
   shifted coordinates fell outside the volume.
2. **Missing boundary voxels.** It sized the candidate cube with `round()`
   rather than `ceiling()` and compared distances exclusively at the boundary.

Checked against brute-force enumeration — every in-bounds voxel whose centre
lies within `radius` mm of the centroid, computed directly — over 384 geometries
spanning isotropic and anisotropic spacing, corner/edge/interior centroids and
boundary radii: **the compiled path was correct in all 384; the fallback was
wrong in 35.**

Fixing it would have meant maintaining a second implementation of something the
compiled path already gets right, so `use_cpp` is now deprecated and ignored:
passing `FALSE` warns and returns the compiled result. The brute-force
comparison is now a permanent test, which is a stronger guarantee than
branch-to-branch equivalence ever was.

## Sequencing and effort

1. **Phase 2** — highest remaining value by a wide margin. ~1 day including the
   golden extension and adversarial review.
2. **Phase 4 scalar loops** — cheap, independent, each ~1 hour with a targeted
   equivalence test. Do the `ClusteredNeuroVec` one first; it is measured and
   trivially verifiable.
3. **Phase 3 separable blur** — ~half a day; the payoff depends on user smoothing
   parameters.
4. **`use_cpp = FALSE`** — decide delete vs fix. Deleting is ~10 minutes.
5. **Dead C++** — ~10 minutes, do it alongside anything else touching `src/`.

Phase 2 alone takes a whole-brain searchlight from ~40 s of enumeration to
roughly 1 s; combined with Phase 1's `series()` work, the end-to-end MVPA
pipeline on the reference workload goes from ~380 s of framework overhead to
roughly 100 s, of which almost all is the memory-bandwidth-bound sparse gather
that cannot be optimised further without changing the data layout.

---

## Environment note (applies to reproducing any of this)

CRAN is unreachable from the container this was developed in, so `RNifti`,
`RNiftyReg`, `bigstatsr`, `mmap`, and `deflist` were replaced with local shims to
make the package installable: `deflist` and `mmap` are functional
reimplementations (they sit on the code paths under test), the other three are
signature-only stubs that error when called. That accounts for all 25 errors in
the baseline — every one is `RNifti::writeNifti`, `RNifti::niftiHeader`, or
`bigstatsr::as_FBM`, and none touch `roi.R`, `searchlight.R`, or `series()`.

Before merging, re-run the golden sweep and the full suite against real
dependencies. Nothing in Phase 1 depends on the shimmed packages, but the
baseline's 25 errors should drop to 0 and that is worth confirming.

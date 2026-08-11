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
   (1,470 assertions), which pin equivalence rather than performance: the cheap
   constructor against `new()` (and against the class validity function, so the
   duplicated invariants cannot drift), the offset template against
   `local_sphere()`, `spherical_roi()` against brute-force enumeration, the
   compiled gathers against the R loops they replaced, the accessor-hook contract
   for lazy sparse subclasses, hostile-input handling for every compiled entry
   point, and a scaling guard that fails if a whole-volume copy is ever
   reintroduced into `series()`.
4. **GC torture** — not part of the suite (it is far too slow), but run by hand
   against every compiled entry point whenever the C++ changes. The last run
   covered `sphere_coords_cpp`, `sphere_roi_at_cpp`, `sphere_coords_batch_cpp`
   and the `.searchlight_plan()` path, with `NA`-bearing volumes and both
   `use_mask` settings, under `gctorture(TRUE)`: all results bit-identical, no
   crash. Re-run it after touching `src/` — a reviewer's earlier pass is not
   transferable once the signatures change.

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

## Phase 3 — done: separable Gaussian blur

**Why.** `gaussian_blur_cpp` applies the full `(2w+1)³` kernel. A Gaussian is a
tensor product, so three 1-D passes give the same answer in `3(2w+1)` taps. The
masked default (`normalize = TRUE`) is separable too, as the ratio of two
separable convolutions — `blur(x·m) / blur(m)` — because the numerator and
denominator the dense code accumulates are exactly those two convolutions.

**What was built.** `src/gaussian_blur_sep.cpp`, and a dispatcher
(`.gaussian_blur_engine()` / `.gaussian_blur_prefers_separable()` in
`R/spat_filter.R`) in front of both former `gaussian_blur_cpp()` call sites.
`x·m` is formed with an explicit branch, not a multiply: out-of-mask voxels are
documented to be allowed to hold `NaN`, and `0 * NaN` is `NaN`, which would leak
missingness into every in-mask neighbour.

Working memory is two vectors of `prod(dim)` doubles, or three plus a byte per
voxel when normalizing. The first 1-D pass reads the masked source through the
mask rather than materialising it, and the denominator pass is ordered to land
in its own buffer, which removed two full-volume sweeps and two of six vectors
from the first working version.

**The plan said to gate on `window` alone. That was wrong.** Keeping the dense
kernel at `window = 1` gives away 1.7-6.8x there, because the crossover depends
on the mask fraction as much as on the window. The first cost model was wrong in
the other direction: gating on tap counts alone
(`n_mask/n_vox > passes/(2w+1)²`) assumes the two kernels cost the same per tap,
and they do not — the separable passes walk contiguous memory with no per-tap
branch, while the dense kernel pays a bounds test and a mask lookup on every
tap. Measuring the crossover on a 91x109x91 volume over mask fractions
0.02-1.00, both contiguous and scattered, put it near
`0.15 * passes / (2w+1)^1.5`. Over that 224-cell grid:

```
rule                        worst case   >5% off   total time
tap counts, passes/(2w+1)^2     4.9x     59/224      7.02 s
0.15 * passes/(2w+1)^1.5        1.4x      7/224      5.51 s
always separable                4.9x     38/224      5.64 s
always dense (status quo)          -          -     48.64 s
per-cell oracle                    -          -      5.48 s
```

The constant is a calibration, not a law — it is a ratio of per-tap costs and
drifts with compiler and machine. It only has to be right near the crossover,
where by construction the two kernels cost nearly the same.

**Measured** (91x109x91 at 2 mm, contiguous mask, `normalize = TRUE`):

```
                     mask 20%  mask 25%  mask 32%  whole volume
window = 1              1.5x      1.7x      2.2x       6.8x
window = 2              2.3x      2.9x      3.7x       9.6x
window = 3              4.0x      5.0x      6.3x      17.3x
```

`normalize = FALSE` roughly doubles those (15x/19x/32x unmasked). The
whole-volume column is what `gaussian_blur(vol)` hits with no mask supplied.

**Verification.** 3,200 configurations of dimension (including degenerate axes),
spacing, sigma, window, mask shape (all / half / single voxel / `NaN` outside)
and `normalize`, compared against the dense kernel: 0 failures, max absolute
difference 2.6e-16. Where the relative difference is largest the output is a
near-zero value produced by cancellation, and against a brute-force reference
the separable result is the *closer* of the two. Golden harness 8,665/0; full
suite unchanged (same 10 pre-existing stub-package I/O errors);
`tests/testthat/test-gaussian-blur-separable.R` adds ~780 assertions covering
equivalence, the mask-insulation guarantee, dispatch monotonicity and argument
rejection.

**What the review found.** The kernel math survived a ~3,400-case differential
fuzz and 18 source mutations. What did not survive was the weaker claim that *the
dispatcher only changes which kernel runs, never the answer* — five reachable
counterexamples, every one a pre-existing defect that the dispatcher turned into
a data-dependent one, because the new kernel validates arguments the old one did
not:

- **A segfault, now mask-fraction dependent.** `gaussian_blur()` never checked
  that `mask` and `vol` had the same dimensions, and the dense kernel built its
  membership lookup with an unbounded `in_mask[mask_idx[i] - 1] = 1`. An
  oversized mask crashed the session — but only below the dispatch threshold,
  since above it the separable kernel's bounds check caught the same index.
  Fixed at both ends: `gaussian_blur()` and `enhance_stat_map()` require matching
  dimensions, and the dense kernel bounds-checks regardless.
- **`sigma` above ~1e203 returned zeros.** `gaussian_weights_impl` multiplied in
  a per-axis `sqrt(2*pi*sigma)` that cancelled exactly in its own normalisation
  but cubed to `Inf` first, making every weight `NaN`. Not formed any more.
- **`gaussian_weights(window >= 645)` killed the session** — `(2w+1)^3` wrapped
  negative in `int`. Guarded; `gaussian_blur()` also clamps the window to
  `max(dim)`, which is result-preserving because every tap that far out is
  already out of bounds.
- **`window = Inf` / `sigma = Inf` were accepted** and produced an all-zero
  volume; after the change they hit an internal error with a message from
  `as.integer()`. Both are now validated at the door with `cli_abort`.
- **`enhance_stat_map()` never validated `spatial_sigma`**, the one numeric
  argument it skipped, so `spatial_sigma = 0` silently returned zeros on one
  path and errored on the other. Validated, along with `intensity_sigma`.

Two divergences were left in place deliberately. A 4-D array is now kept on the
dense path by the dispatcher rather than rejected, so behaviour is unchanged
there. And at `window >= 23` with a non-finite value in the data, the dense
kernel's three-axis weight product underflows to exactly zero, so `0 * Inf`
poisoned the output to `NaN` while the separable form propagates the `Inf`;
neither answer means anything and reproducing an underflow artefact would be
wrong. It is in NEWS.

The review also found the test file's weak spots: the mask-insulation test used a
*constant* in-mask signal, which any normalised kernel reproduces regardless of
its shape, padding or axis order, and 404 of its assertions were algebraic
identities of the dispatch formula that would hold for any formula of that shape
— including one with a badly wrong constant. The signal is now structured and
swept over `NaN`/`NA`/`Inf`/`-Inf` outside the mask, and the dispatch test pins
the fitted thresholds numerically so re-calibration has to be deliberate. Cases
the mutation testing showed were unreachable — `window > 4`, windows larger than
the volume, integer `arr`, duplicated/unsorted/double `mask_idx`, the legacy
branch with non-finite input, and `.gaussian_blur_engine` itself — are now
covered.

**Risk.** The remaining exposure is the calibration constant, which costs at most
~1.4x on calls in the tens of milliseconds when it is wrong. The reviewer's
independent measurement put the `window = 1`, `normalize = TRUE` crossover nearer
0.20 than the 0.15 measured here, against a threshold of 0.173 — so a typical
brain mask can land in a band where the rule fires up to ~1.2x early. That spread
between two machines is the honest width of the calibration.

---

## Phase 4 — done: scalar R loops and hygiene

Five independent loops, each replaced by a vectorised equivalent and pinned
against a reference that reproduces the loop it replaced.

| Item | File | What it was | Measured |
|---|---|---|---|
| `ClusteredNeuroVec[i,j,k,m]` | `R/clustered_neurovec.R` | nested loop over timepoints × voxels | **39×** at 10,000 voxels × 100 timepoints |
| `mapf(NeuroVol, Kernel)` | `R/neurovol.R` | `sum(x[cbind(ii,jj,kk)] * wts)` per centre | **40×** (64³, 3×3×3), **14×** (5×5×5) |
| `random_searchlight()` | `R/searchlight.R` | O(remaining) free-list rebuild per iteration | **4.7×** on a 289,148-voxel mask |
| `NeuroHyperVec[i,j,k,l,m]` | `R/neurohypervec.R` | per-voxel subset + `aperm` | **1.8×** at 64,000 voxels |
| `as.dense(NeuroHyperVec)` | `R/neurohypervec.R` | fresh vector + array per (feature, trial) | **1.7×** |
| Dead compiled code | `src/` | 5 registered-but-uncalled entry points | ~29 KB of source gone |

**The roadmap's plan for `mapf` was wrong.** It said to wire up the existing
unused `kernel_filt_3d_cpp`. That function is a **2-D** filter — `NumericMatrix`
data, `NumericMatrix` kernel — despite the `3d` in its name, so it cannot express
a 3-D kernel convolution at all. It was deleted instead. The real fix needed no
C++: a kernel tap is a fixed coordinate offset, so in linear-index terms it is a
fixed *scalar* offset, and the convolution becomes one vectorised pass per tap
(a few dozen) instead of one R iteration per centre (up to a million) — the same
translate-a-template identity the ROI builders use. A first version tested the
per-centre bounds on every tap, which cost most of the win back at a 5×5×5 kernel
(2.9× rather than 14×); whether a tap can leave the volume is decided by the
extreme centres, so it is three scalar comparisons per tap.

**`mapf`'s boundary handling was broken and nobody had noticed.** With any mask
reaching within the kernel's radius of the boundary, a coordinate past the far
face raised `subscript out of bounds`, and one before the near face made
`x[cbind(ii,jj,kk)]` return a *shorter* vector — R drops zero subscripts silently
— which then recycled against the weights and produced a wrong number with only
a recycling warning. Verified on a `1 4 4` centre: 0.557 returned against 0.371
correct. Out-of-volume taps are now skipped. Its dimension check was also
`!all.equal(...)`, which raises `invalid argument type` on a mismatch because
`all.equal` returns a character description rather than `FALSE`.

I first claimed this affected only the masked path, "since unmasked centres are
held clear of the border". **That was wrong, and the reviewer disproved it.**
`hwidth[d]:(dims[d] - hwidth[d])` counts *backwards* when `dims[d] < 2*hwidth[d]`,
so for a volume thinner than the kernel the unmasked grid contains centres inside
the margin — on a 10 × 3 × 10 slab with a 3×3×3 kernel the j range is `2:1` — and
86 voxels change value there too. The new values are the correct clipped ones (the
old ones were recycled garbage), but the claim was false and is now stated
properly in NEWS, with a slab case in the test file. The descending range is a
pre-existing wart in `mapf`'s grid construction, independent of this change: the
honest reading is that such a volume has *no* valid interior centres. Left alone
because changing which voxels get written is a larger behaviour change than the
one being made here.

**`random_searchlight` could not be fixed without changing its output**, and that
is worth stating plainly. To keep the realized sequence for a given seed you must
sample `sample.int(length(remain_indices), 1)` from an *ascending* free list, and
removing from a sorted R vector is O(n) however you do it — the asymptotics are
inherent to the ordering, not to the implementation. The swap-with-last free list
gives O(batch) removal but leaves the list unordered, so a given seed picks
different centres from the second iteration onward. Sampling remains uniform over
the unclaimed voxels, so this is a re-realisation, not a behaviour change: over 60
seeds the searchlight count was 466.7 ± 8.6 before and 463.6 ± 9.0 after, with
*exactly* the same voxels claimed in every run, and no defect (duplicate centre,
double-claimed voxel, out-of-mask coordinate) in either. The first searchlight of
a run is identical, since no removals have happened yet — which is why the
existing seeded test keeps passing. Logged as breaking in NEWS.

**What else the review found.** `drop_voxels` had unasserted preconditions: the
logical vector it replaced was idempotent and duplicate-safe, a swap-with-last
list is neither, and violating either silently desynchronises `n_free` from the
live set so that voxels go unclaimed. Live callers do satisfy the contract — the
reviewer verified the full invariant set over 2,616 randomised runs across ten
mask shapes — but it now rests on `sphere_at_cpp` never emitting a duplicate
coordinate, so the contract is checked rather than assumed. Also fixed: the new
`ClusteredNeuroVec[` bounds check rejected `6.5` while accepting `2.7` on a
6-voxel axis, since it tested before truncating (it now truncates first, as R's
own array indexing does); the matrix subset passed a negative `m` through to R's
drop-these-rows semantics where the scalar loop errored; and `mapf`'s all-integer
`centre_lin` arithmetic would overflow to `NA` past 2^31 voxels, where the old
coordinate gather had no ceiling.

The review's mutation testing caught 23 of 23 genuine mutations across all five
sites, but flagged four weak assertions, all since replaced: a reproducibility
test that could not fail (same seed twice passes for any RNG-driven
implementation — it now also asserts that a *different* seed gives a different
partition), a `skip_if` that would silently retire the cluster-0 case, a
dead-code check that only tested the R wrappers and not the native registration
table, and the slab case above.

**Verification.** The golden harness was extended to 9,163 recorded outputs
covering all five sites and captured against the pre-change build: 9,149 identical,
and all 14 differences are `random_searchlight` realisations (`n`, `n_covered`,
`total_voxels`) — the invariant components matched. The free-list data structure
was separately fuzzed over 14,992 removal steps against a naive set reference,
checking set equality, the live count, that `pos` is a bijection onto
`1..n_free`, and that `free[pos[v]] == v`. `tests/testthat/test-vectorized-loops.R`
adds ~410 assertions, each comparing against a reference implementation of the
replaced loop.

Two measurement traps caught here, both from timing large allocations:
`dense_array_5d` first appeared to be 15–22× faster and `ClusteredNeuroVec[` only
1.9×; with `gcFirst = TRUE` and a median of seven runs the honest figures are 1.7×
and 39×. The reference implementation even measured *faster* at a larger size in
one run, which is the tell. Anything allocating hundreds of MB needs gc control
before its ratio means anything.

**Dead code.** `radius_search_3d_{nonisotropic,direct,precomputed}`,
`local_spheres` and `kernel_filt_3d_cpp` were compiled and registered but called
from nowhere in `R/`, `tests/` or `inst/`. Removed. `local_sphere` (singular) was
**kept**: it is the reference implementation that the offset-template equivalence
test in `test-roi-series-fastpaths.R` compares against.

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

## Status

All four phases are done, `use_cpp = FALSE` is resolved (deprecated), and the
dead compiled code is gone. What remains is one low-value item:

1. **`.sphere_offset_cache` sizing** — bound by total bytes rather than entry
   count. A pathological radius can park hundreds of MB in the namespace; see
   the note at the end of Phase 2.

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
signature-only stubs that error when called. That accounts for all 10 failing
tests in the baseline — 8 `RNifti::writeNifti`, 1 `RNifti::niftiHeader`, 1
`bigstatsr::as_FBM` — and none touch `roi.R`, `searchlight.R`, `series()` or
`spat_filter.R`. The set is byte-identical before and after every change on this
branch, which is how "no regressions" was checked.

Before merging, re-run the golden sweep and the full suite against real
dependencies. Nothing here depends on the shimmed packages, but those 10 errors
should drop to 0 and that is worth confirming.

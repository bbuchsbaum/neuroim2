# neuroim2 0.18.0

## Bug Fixes

* **Breaking:** Removed the exported S4 generic `scale()`. It had no methods whatsoever, so attaching neuroim2 masked `base::scale()` and broke it for *every* input, not just neuroimaging objects (`scale(c(1, 2, 3))` errored with "unable to find an inherited method"). The package itself already called `base::scale()` internally to work around the mask. `scale()` now resolves to base R again; use `scale_series()` for the neuroimaging-specific scaling generic (GitHub #16).

## NIfTI I/O rewritten on a compiled core

The read and write paths no longer go through `readBin`/`writeBin`. A single
translation unit (`src/nifti_data_io.cpp`) converts between the file's stored
type and `double` in one pass, for plain and gzipped files alike, and both
`read_vol()`/`read_vec()` and `write_vol()`/`write_vec()` now call it. On this
machine, against `nibabel` 5.4.2 reading the same files:

```
                                   before     after    nibabel
read 3-D int16 .nii  (14 MB)        0.46 s    0.049 s    0.027 s
read 3-D float32 .nii.gz            0.45 s    0.220 s    0.247 s
read 4-D int16 .nii  (56 MB)        3.09 s    0.371 s    0.107 s
read 4-D int16 .nii.gz              2.81 s    0.750 s    0.585 s
```

Two things were behind the 4-D case in particular. `read_mapped_vols()` built a
236 MB vector of `double` indices via `outer()` simply to enumerate a
*contiguous* range, gathered it element by element through the memory map, and
handed back a time-by-voxel matrix that `DenseNeuroVec()` then transposed
straight back. It now reads whole volumes sequentially and returns them
voxels-by-volumes, the layout the object stores. And building the object no
longer copies the payload: an S4 class extending `array` *is* that array with
the slots as attributes, so the read path sets them in place, the same
technique the ROI constructors already used. The equivalence to `new()` is
pinned by `identical()` in the test suite.

`read_mapped_vols()` is internal, but its return layout changed from
`[volumes x voxels]` to `[voxels x volumes]`. Two callers were also missing
`scl_slope`/`scl_inter` scaling entirely; both are fixed.

## Every NIfTI datatype the standard defines

`int8`, `uint16`, `uint32`, `int64` and `uint64` are read and written. They
were rejected before with *"Unsupported NIfTI data-type code"* — `uint16` in
particular is not exotic.

`int32` images no longer lose voxels holding `-2147483648`. That bit pattern is
R's `NA_integer_`, so routing the payload through an R integer vector turned it
into `NA` and dropped it from every downstream statistic; the compiled reader
lands in `double` and never sees an R integer.

Complex (`DT_COMPLEX64/128/256`), colour (`DT_RGB24`, `DT_RGBA32`),
`DT_FLOAT128` and 1-bit `DT_BINARY` are still not read, but the error now names
the type and says why: a `NeuroVol` holds one real number per voxel, and
inventing a projection — a modulus, a luminance — would be answering a question
the caller did not ask.

## Writing preserves the header

`as_nifti_header()` built every header from defaults plus the object's
geometry, so a read-write round trip silently discarded everything else. It now
starts from the header the image was read from. Concretely, a round trip used
to lose:

* **the repetition time** (`pixdim[4]` reset to 0) and `xyzt_units`,
* **`sform_code`**, rewritten from whatever it was to 1 — so an image in MNI
  space came back labelled scanner-anatomical,
* `descrip`, `aux_file`, `intent_code`, `intent_name`, `cal_min`/`cal_max`, the
  slice-timing fields and `toffset`.

All of them survive now. The source header travels on a new `header` slot on
`NeuroObj`, which every image class inherits; `header(x)` reads it, and that
method now accepts images as its documentation always claimed.

`write_vol()`/`write_vec()` also default to the source file's datatype rather
than `FLOAT`, when the values are still exactly representable in it. An
unmodified `int16` image round-trips byte for byte instead of doubling in size;
computed data that no longer fits is promoted to `FLOAT` rather than quietly
quantised.

### Integer output is scaled, not truncated

`write_vol(x, f, data_type = "SHORT")` went through
`writeBin(as.integer(x), ...)`, which truncates toward zero without a warning:
data spanning ±3.7 came back with a maximum error of **0.998**. `scl_slope` and
`scl_inter` are now derived to fit the target type, giving **5.4e-05** on the
same data — the same precision `nibabel` produces. Data that is already
integral and in range is still written unscaled, so masks, label volumes and
atlases stay exact.

## NIfTI-2

Read and written, plain or gzipped. `write_vol()`/`write_vec()` take
`version = 2` or `format = "NIFTI2"`, and
select it automatically when a dimension exceeds 32767, which NIfTI-1's 16-bit
`dim` field cannot represent. Files written this way are read correctly by
`nibabel`, including the affine, TR and units.

## Conformance and diagnostics

* **`scl_slope = NaN` no longer aborts the read.** It is one of the three legal
  spellings of "this file is not scaled", alongside `0` and an absent field, and
  it is what `nibabel`'s in-memory header carries by default — so neuroim2 was
  refusing to open files the reference implementation writes. An intercept
  recorded next to a dead slope is dropped with it.
* **Singleton axes are preserved.** The number of axes is taken from `dim[0]`,
  as the standard says, instead of being inferred by stopping at the first
  extent equal to 1. A 64 × 64 × 1 single-slice volume used to load as 64 × 64,
  and a one-volume run as a 3-D volume.
* **5-D files keep their shape.** `read_nifti_header()` used to fold a
  degenerate fourth axis and warn, so `read_header()` and `read_hyper_vec()`
  could not see what the file actually declared. The fold now happens only
  inside `read_vec()`, where it is a documented convenience; a genuinely 5-D
  file is refused there with a pointer to `read_hyper_vec()`.
* **`qform_code == 0` and `sform_code == 0` now follow METHOD 1.** Both
  transforms being marked unknown means the quaternion is not a statement about
  where the image is; the affine is derived from `pixdim` with the image centred
  on the origin, which is what the standard specifies and what `nibabel`
  returns. This also gives ANALYZE 7.5 files their correct affine instead of
  `diag(2, 2, 2)` with a zero translation.
* **A truncated file says what is missing**: how many elements the header
  describes, at what offset, and how many bytes the file actually supplies. It
  used to fail with `argument is of length zero`. A file that is not an image at
  all now reports its leading `sizeof_hdr` field and what that field must be.
* `intent_name` is read as the fixed-width byte field it is. Reading it as a C
  string ran past the field whenever those 16 bytes held no null terminator.

## Filtering

`gaussian_blur()` accepts a `NeuroVec` and smooths it volume by volume, with
the argument checks, kernel and mask resolved once instead of per volume.

The separable kernel also got two exact optimisations. Its `y` and `z` passes
accumulate over whole contiguous rows and slices rather than walking a column
with stride `d0` (or `d0 * d1`), which touched a fresh cache line per
iteration; the additions for each output element still happen in the same
order, so the result is bit-identical. And when the mask covers every voxel —
what `gaussian_blur(vol)` with no mask means — the mask-insulation denominator
is the zero-padded convolution of a constant, which factorises into three
per-axis edge corrections, replacing three full-volume passes with `d0 + d1 +
d2` additions.

```
                                          before     after   nibabel/scipy
3-D, FWHM 6 mm, no mask                   0.60 s    0.21 s      0.13 s
4-D, FWHM 6 mm, 200 volumes               2.04 s    0.97 s      0.57 s
```

Equivalence to the dense kernel is checked over 432 configurations of
dimension, window, sigma, spacing, mask and `normalize` in the test suite.

## Documentation

* Clarified that `gaussian_blur()`'s `sigma` is expressed in the spatial units of the image (millimetres for typical NIfTI data), **not** voxels, while `window` remains a voxel half-width. Mixing the two silently under-smooths: an FWHM converted to sigma in voxels applies roughly half the intended smoothing with no error or warning. Added a Details section on units and an FWHM → sigma conversion example (GitHub #23).
* `vignette("large-data")` and `?read_vec` now recommend the masked and
  memory-mapped paths for anything large, rather than presenting them as a
  fallback for when memory runs out. For the scattered access that searchlight,
  ROI and connectivity work perform they are also *faster*: 5,000 voxel time
  courses from a 64 × 64 × 36 × 200 run measured 44 ms with `mask=` and 144 ms
  via `mode = "mmap"`, against 4.8 s for loading the image and then indexing it.

## Other changes

* `as.array()` on a dense image returns a plain array. It used to hand back the
  payload still carrying the `space`, `header` and `class` attributes, because
  the slots of an S4 class extending `array` *are* attributes on that array.
* `linear_access()` on a `FileBackedNeuroVec` gathers in one vectorised step
  instead of looping over the request in R, and applies data scaling, which it
  previously ignored.
* `dev/bench/nibabel/` holds the conformance and timing harnesses these numbers
  come from; `dev/nibabel-gap-analysis.md` is the write-up. 30 probes covering
  the standard now agree with `nibabel` on 28, with the remaining two refused by
  design.

# neuroim2 0.17.0

## Testing

* The test suite now runs clean: 2969 passing, no failures and no errors.
  Previously `R CMD check` reported 7 failures and 2 errors on every platform.
* `test-plot-registration-qc.R` called `ggplot2::get_labs()`, which only exists
  from ggplot2 3.5.2 while DESCRIPTION sets no version floor, so it errored on
  older ggplot2. Tests now read plot titles through a version-agnostic helper.
* The vdiffr golden-image tests have moved to `dev/visual-snapshots/`. They
  compare SVG text byte-for-byte, which encodes the font metrics of the machine
  that produced the snapshot, so a single stored image cannot match the
  Windows, macOS and Linux check matrix — and testthat deletes any `_snaps/`
  directory whose test file did not run, so they could not simply be skipped.
  `tests/testthat/test-plot-structure.R` covers the same plotting calls with
  platform-independent structural assertions (panel counts, layer counts, data
  shape, argument validation) and runs on every check. See the README in
  `dev/visual-snapshots/` for reviewing intentional visual changes.

## Documentation

* The vignettes have been rewritten as nine articles in three tiers: *Learn*
  (`neuroim2`, `spaces-and-coordinates`, `volumes-and-vectors`,
  `reading-and-writing`), *Do* (`regions-and-searchlights`,
  `resampling-and-orientation`, `smoothing-and-filtering`, `visualization`) and
  *Scale* (`large-data`), replacing the previous fourteen. Old article URLs
  redirect to their successors on the package website; `vignette()` calls using
  the old names will need updating.
* Vignette examples no longer take time series from `global_mask_v4.nii`, which
  is a binary mask repeated across four timepoints — every "time series" drawn
  from it was a constant vector. Anything with a time axis now uses
  `simulate_fmri()`, so time courses have real temporal structure and the
  smoothing, searchlight and ROI examples demonstrate measurable effects.
* Vignette YAML placed `css` and `includes` at the document root, where
  `html_vignette` ignores them, so the package stylesheet never reached the
  rendered articles. They are now nested under the output format.

## Performance

* Searchlight and ROI extraction are substantially faster. Each of the three
  changes below was landed and verified as output-preserving on its own —
  compared against a build of the previous version across 8,625 configurations
  of `spherical_roi()`, `spherical_roi_set()`, `series()` and the four
  searchlight iterators — before the correctness fixes listed under *Bug Fixes*
  deliberately changed some of those outputs. The speed-ups and the behaviour
  changes are independent.

  - ROI objects are no longer built with `new()`. `ROIVolWindow` reaches a basic
    type through two levels of `contains=`, so the default `initialize()` walked
    that chain with `callNextMethod()` and ran `validObject()` at each level,
    recursing into the `NeuroSpace` slot's nested axis tree — about 450 us per
    object regardless of ROI size, which was roughly 70% of an exhaustive
    searchlight's runtime. An internal constructor clones a cached prototype and
    asserts the class invariants once, producing an `identical()` object about
    11x cheaper.
  - The voxel offsets of a spherical neighbourhood depend only on the radius and
    the voxel spacing, so they are computed once per `(radius, spacing)` pair and
    reused; per centre the work is a translate-and-clip. Candidate windows are
    also sized per axis rather than from `min(spacing)`, which scanned about 3x
    too many candidates on anisotropic voxels.
  - `series()` gathers in a single compiled pass.
    `series(DenseNeuroVec, matrix)` previously looped over timepoints,
    rebuilding a double index vector each pass (~17x faster);
    `series(AbstractSparseNeuroVec, integer)` allocated a zero matrix and
    scattered into it (~1.8x faster, close to memory-bandwidth-bound). The
    sparse fast path is gated on the concrete `SparseNeuroVec` class so
    subclasses that override `matricized_access()` keep their hook.

  - The searchlight iterators do the whole per-centre job in one compiled pass --
    translate the cached offset template, clip to the volume, apply the mask,
    gather values, locate the centre row -- with everything invariant hoisted out
    of the iterator closure. `searchlight_coords()` in particular used to build a
    full `ROIVolWindow` per centre and then discard all but its coordinates.
    `spherical_roi_set()` expands every centre in one batched call, so
    `searchlight(eager = TRUE)` benefits too.
  - The iterators remain lazy. Materialising a whole-brain radius-8 mm
    searchlight at once would be about 0.9 GB of coordinates, or 1.5 GB as ROI
    objects, so the work went into making the per-element call cheap rather than
    into emitting everything up front.

  On a 91x109x91 volume at 2 mm with a 290,137-voxel mask and radius-8 mm
  searchlights, per element and projected over the whole brain:

  ```
  searchlight_coords()        97.0 -> 7.7 us     28.1 s -> 2.2 s
  searchlight(eager = FALSE)  92.3 -> 28.0 us    26.8 s -> 8.1 s
  spherical_roi()            660.0 -> 85.0 us   191.5 s -> 24.7 s
  ```

* `gaussian_blur()` now evaluates the kernel separably where that is cheaper. A
  Gaussian is a tensor product, so the full `(2 * window + 1)^3` kernel is
  equivalent to three 1-D passes --- 15x fewer multiplies at `window = 3`. The
  masked default (`normalize = TRUE`) is separable too, as the ratio of two
  separable convolutions: the mask-insulated blur is `blur(x * m) / blur(m)`,
  where `m` is the mask indicator.

  Which form is cheaper depends on both the window and the mask fraction --- the
  separable passes cover the whole volume, while the dense kernel visits mask
  voxels only --- so the two implementations are both kept and the cost is
  estimated per call. Measured on a 91x109x91 volume at 2 mm across mask
  fractions from 2% to 100%, the dispatch rule picks the slower kernel by more
  than 5% in 7 of 224 cases and never by more than 1.4x, and its total time is
  within 0.6% of always choosing the per-case winner.

  For a brain-mask-shaped input (about a quarter of the bounding box) and for
  unmasked whole-volume smoothing:

  ```
                       mask 25%   whole volume
  window = 1              1.7x         6.8x
  window = 2              2.9x         9.6x
  window = 3              5.0x        17.3x
  ```

  Speed-ups are larger with `normalize = FALSE` (3.9x / 5.5x / 9.4x masked,
  15x / 19x / 32x unmasked), which needs half as many passes. The whole-volume
  column is the case `gaussian_blur(vol)` hits when no mask is supplied.

  Results are unchanged to within floating-point summation order: over 3,200
  configurations of dimension, spacing, sigma, window, mask shape and
  `normalize`, the largest absolute difference from the previous kernel was
  2.6e-16, and where the two differ most in relative terms the separable result
  is the closer of the two to a high-precision reference. Code comparing blurred
  volumes with `identical()` rather than `all.equal()` may see the difference.
  The separable path needs working memory the dense one did not: two vectors of
  `prod(dim)` doubles, or three plus a byte per voxel when normalizing. Measured
  peak-RSS delta on a 200x200x200 volume at `window = 1`, `normalize = TRUE`:
  153 MB against the dense path's 77 MB at a 20% mask, 167 against 139 at 50%,
  and 190 against 230 at a full mask --- where the separable form is the cheaper
  of the two, because the dense path's voxel-coordinate matrix costs 24 bytes per
  mask voxel. The scratch is `malloc`ed rather than allocated by R, so it is not
  counted against `mem.maxVSize` and exhaustion surfaces as a C++ allocation
  failure rather than R's "cannot allocate vector of size".

  One input does change materially, and only there: at `window >= 23` the dense
  kernel's three-axis weight product underflows to exactly zero in the far
  corners, so an `Inf` in the data multiplied to `NaN` and poisoned the output.
  The separable form applies the three factors in separate passes, none of which
  underflow, and propagates the `Inf`. Neither answer is meaningful --- the input
  is not finite --- but they differ.

* Five scalar R loops over voxels are now vectorised. Each was verified against
  a reference that reproduces the loop it replaced, and against a build of the
  previous version across 9,163 recorded outputs.

  ```
                                                       before    after
  ClusteredNeuroVec [i, j, k, m]   10,000 voxels     0.194 s   0.005 s   39x
  mapf(NeuroVol, Kernel)           64^3, 3x3x3        3.81 s   0.096 s   40x
                                   64^3, 5x5x5        4.10 s   0.288 s   14x
  random_searchlight()             289,148 voxels    16.8 s    3.6 s      4.7x
  NeuroHyperVec [i, j, k, l, m]    64,000 voxels      0.342 s  0.195 s    1.8x
  as.dense(NeuroHyperVec)          40^3 x 40 x 10     1.10 s   0.653 s    1.7x
  ```

  - `mapf()` evaluated `sum(x[cbind(ii, jj, kk)] * weights)` once per centre. A
    kernel tap is a fixed coordinate offset, so in linear-index terms it is a
    fixed scalar offset: the convolution is now one vectorised pass per tap ---
    a few dozen --- rather than one R iteration per centre, of which there can be
    a million. Whether a tap can leave the volume at all is decided by the
    extreme centres, so that is a scalar test per tap and not a per-centre one.
  - `random_searchlight()` rebuilt its list of unclaimed voxels with
    `remain_indices[remaining[remain_indices]]` on every iteration, which is
    quadratic in the mask size. It now keeps a free list whose removals cost
    O(batch): survivors from the tail are swapped into the holes left behind.
  - `as.dense()` on a `NeuroHyperVec` built a fresh volume-sized vector and a
    fresh 3-D array for every (feature, trial) pair; it now writes the in-mask
    positions of all trials of a feature in one scatter. What is left is
    allocating and filling the dense result itself, which is why the speed-up is
    modest --- that result is 195 MB for the largest case above.


## Bug Fixes

* **Breaking:** `random_searchlight()` returns a different set of searchlights
  for a given random seed. Sampling is still uniform over the unclaimed voxels,
  so the distribution of results is unchanged --- over 60 seeds the number of
  searchlights was 466.7 +/- 8.6 before and 463.6 +/- 9.0 after, with exactly the
  same voxels claimed in every run --- but the realized sequence differs because
  the free list is no longer held in ascending order. The first searchlight of a
  run is unaffected. Analyses that recorded a seed and expect the identical
  partition will need to re-run.

* **Breaking:** `mapf()` now clips the kernel at the volume boundary instead of
  failing there. Taps falling outside the volume contribute nothing. Previously
  a coordinate past the far face raised `subscript out of bounds`, and one before
  the near face made the gather return a *shorter* vector, which then recycled
  against the weights and produced a silently wrong value with only a recycling
  warning. Two cases are affected: any `mask` reaching within the kernel's radius
  of a face, which made the masked form unusable; and an unmasked volume thinner
  than twice the kernel half-width, because the centre range
  `hwidth:(dim - hwidth)` counts *backwards* when `dim < 2 * hwidth` and so
  yields centres inside the margin rather than none --- on a 10 x 3 x 10 slab with
  a 3x3x3 kernel, 86 voxels change. Interior results are unchanged.

* `mapf()`'s mask dimension check was written `!all.equal(dim(mask), dim(ovol))`.
  `all.equal()` returns a character description rather than `FALSE` on a
  mismatch, so the guard raised `invalid argument type` instead of its own
  message. It now reports "mask must have same dimensions as input volume".

* `ClusteredNeuroVec[i, j, k, m]` now rejects voxel coordinates outside the
  volume. They previously reached an `if (NA > 0)` and failed with "missing value
  where TRUE/FALSE needed".

* Removed `radius_search_3d_nonisotropic()`, `radius_search_3d_direct()`,
  `radius_search_3d_precomputed()`, `local_spheres()` and `kernel_filt_3d_cpp()`.
  All five were compiled and registered but called from nowhere in the package;
  `kernel_filt_3d_cpp()` was in any case a 2-D filter despite its name.

* **Breaking:** `nonzero = TRUE` now drops voxels whose value is `NA` or `NaN`,
  in every ROI builder and every searchlight iterator. Previously the filter was
  written `vals != 0`, which evaluates to `NA` for a missing voxel; an `NA`
  logical subscript emits an all-`NA` row, so a `ROIVolWindow` could come back
  with `NA NA NA` coordinates. Missing values also made the eager and lazy
  searchlight iterators disagree with each other: on a 3x3x3 volume with one `NA`
  voxel they differed on 21 of 25 searchlights. "Nonzero" cannot be established
  for a missing value, so it is dropped. The rule now lives in one place, in the
  compiled searchlight core.

* Fixed a `SIGFPE` crash in `index_to_grid()` for spaces whose running product of
  dimensions exceeds `int`: `index_to_grid(NeuroSpace(c(65536L, 65536L, 2L)), 1:3)`
  killed the R session. The stride accumulator wrapped to zero and the next
  division faulted. It is now 64-bit.

* `searchlight()` and `searchlight_coords()` again reject a `radius` smaller than
  the finest voxel dimension, as the single-ROI builders always have. The lazy
  iterators had stopped validating it and returned degenerate one-voxel windows
  instead.


* Fixed a segfault in `gaussian_blur()` when `mask` and `vol` had different
  dimensions. Nothing checked that they matched, and the kernel built its
  membership lookup by writing at `mask_idx` into a vector sized for the volume,
  so an index past the end was an unbounded write. `gaussian_blur()` now requires
  matching dimensions, and the kernel bounds-checks the index regardless.
  `enhance_stat_map()` gained the same dimension check.

* Fixed `gaussian_blur()` returning an all-zero volume for `sigma` above about
  1e203. The kernel carried a per-axis `sqrt(2 * pi * sigma)` factor that
  cancelled exactly in the normalization but cubed to `Inf` first, and `Inf/Inf`
  made every weight `NaN`. The factor is no longer formed, which leaves the
  normalized kernel unchanged for every finite `sigma`.

* Fixed a session-killing crash in `gaussian_weights()` for `window >= 645`,
  where `(2 * window + 1)^3` wrapped negative in `int` and allocated a vector of
  negative length.

* `gaussian_blur()` now rejects a `sigma` or `window` that is missing,
  non-finite, or not length one; previously `window = Inf` and `sigma = Inf` were
  accepted and produced an all-zero volume. A `window` larger than the volume is
  clamped to `max(dim(vol))`, which changes no result --- every tap that far out
  is out of bounds from every voxel --- and keeps absurd values out of the kernel
  builders. `enhance_stat_map()` now validates `spatial_sigma` and
  `intensity_sigma`, which were the two numeric arguments it did not check.

* **Breaking:** `spherical_roi()`, `cuboid_roi()` and `square_roi()` now share
  one definition of centroid handling, `nonzero` filtering and the
  centre/parent bookkeeping. Each previously implemented these separately and
  the three disagreed with one another; where they disagreed, at least one was
  wrong.

  - `nonzero = TRUE` now means what it documents. `cuboid_roi()` and
    `square_roi()` used to force the centre voxel back into the ROI after
    filtering, so a "nonzero" region could contain a zero value; they had done
    this only so that `parent_index` could be read back off the coordinate
    matrix. `parent_index` is now derived from the requested centroid directly,
    so the workaround is gone and a zero-valued centre is dropped like any
    other zero voxel.
  - `center_index` is the row of the centre voxel in `coords`, or `NA_integer_`
    when the centre is not part of the ROI. `spherical_roi()` used to fall back
    to `1L`, silently pointing at an unrelated voxel, and then computed
    `parent_index` *from that voxel* — so a filtered ROI reported the wrong
    parent index entirely.
  - `parent_index` is always the linear index of the requested centroid in the
    parent space, whatever the filtering, matching what the slot documents.
  - All three now validate the centroid identically: length 3 (or a 1x3
    matrix), `>= 1`, and within the volume. `cuboid_roi()` and `square_roi()`
    previously accepted out-of-bounds centroids and proceeded with a warning
    from `matrix()`.
  - A fractional centroid such as `c(5.7, 5, 5)` is truncated once and the same
    value is used both to build the neighbourhood and to locate the centre
    within it. Previously the grid was built from the truncated value while the
    centre was matched against the original, so the centre was never found.
  - `coords` is now an integer matrix with no dimnames, and rows are ordered
    lexicographically by `(x, y, z)` in every builder. `cuboid_roi()` and
    `square_roi()` returned `expand.grid()` output, which varies `x` fastest
    and so was ordered the opposite way, and carried `x`/`y`/`z` dimnames.

* **Breaking:** The `use_cpp` argument of `spherical_roi()` is deprecated and
  ignored; the compiled implementation is always used, and passing `FALSE`
  warns. The `use_cpp = FALSE` branch was wrong in two independent ways: it
  returned coordinates offset by **+1 voxel on every axis** (it produced
  1-based coordinates where the internal contract is 0-based, and the caller
  then added 1 regardless), and it sized its candidate cube with `round()`
  rather than `ceiling()` while comparing distances exclusively at the
  boundary, so it also dropped legitimate boundary voxels. Checked against
  brute-force enumeration over 384 geometries — every in-bounds voxel within
  the radius, computed directly — the compiled path was correct in all 384 and
  the fallback wrong in 35. The offset also made `spherical_roi()` throw
  `subscript out of bounds` for centroids near the volume edge; those calls now
  succeed.


* **Breaking:** Fixed `reorient()`, which never returned the orientation it was
  asked for. It applied a transform derived only from the target axis codes,
  ignoring the space's current orientation, and it interpreted those codes with
  the sense inverted, so `reorient(sp, c("R", "A", "S"))` produced an `LPI`
  space and asking for an image's existing orientation was not a no-op.
  Orientation codes are now read the way `affine_to_axcodes()` reports them —
  naming the direction each axis increases *towards* — and the transform is
  computed relative to the current orientation. Code that compensated for the
  old inversion by requesting the opposite codes will need updating.
* **Breaking:** Fixed a half-voxel offset in `index_to_coord()` and
  `coord_to_index()`. Both used a 0.5-voxel shift where `grid_to_coord()` and
  `coord_to_grid()` use the 1-based-to-0-based offset of 1, so the shortcut
  conversions disagreed with the equivalent two-step path by half a voxel and
  `coord_to_index()` could return a neighbouring voxel. Voxel `(1, 1, 1)` now
  maps exactly to `origin()` by either route. Raster extents produced by the
  plotting helpers shift by half a voxel as a consequence, placing voxel
  centres on their affine positions.
* Fixed an off-by-one in world-to-index conversion on anisotropic and
  large-origin grids. `grid_to_index()` truncates, and `NeuroSpace` stores
  affines to 7 significant figures, so a coordinate that is arithmetically an
  exact voxel centre could arrive as `2.9999985` and truncate to the voxel
  below — on the shipped EPI mask the error reached 1.6e-6 grid units.
  `coord_to_index()` now rounds to the nearest voxel centre, and
  `grid_to_index()` snaps values within 1e-4 of an integer. Genuinely
  fractional grid coordinates still truncate, matching R's array indexing.
  (The previous 0.5-voxel offset was partly masking this, since `trunc(x + 0.5)`
  rounds; correcting the offset alone exposed it.)
* Fixed `linear_access()` on sparse `NeuroVec` objects, which could return wrong values or error when there were more masked voxels than time points. The sparse data matrix is stored as `[time × voxel]`; the linear-index path now indexes rows and columns in the correct order.
* Fixed `ROIVol` arithmetic so values are aligned by voxel index, not by the original coordinate order, and fixed sparse `Summary` group methods so reductions include implicit structural zeros. `ROIVol` arithmetic now documents and enforces a same-support contract (identical space and voxel set; order may differ); missing voxels are not treated as zero.
* Hardened `simulate_fmri()` edge cases: zero FWHM values now disable the corresponding smoothing step, `n_time = 1` no longer trips AR loops, tiny or constant masks no longer produce `NA` heteroscedasticity fields, and scalar arguments now receive explicit validation.

## Improvements

* Faster data access for dense and sparse neuroimaging objects: `DenseNeuroVol` and `DenseNeuroVec` subsetting now indexes the backing array directly instead of materialising full spatial grids via `expand.grid()`; `ArrayLike3D` builds column-major linear indices with vector arithmetic; and `linear_access()`, `lookup()`, and `validate_indices()` use lighter bounds checks (`anyNA` + `range()` instead of allocating large logical vectors).
* `NeuroVec` arithmetic now uses matrix-level operations instead of per-volume S4 dispatch, `NeuroVec` comparisons work for sparse-backed inputs, and temporal `mean()` methods now honor `na.rm`.

# neuroim2 0.16.0

## New Features

* Added `enhance_stat_map()`, a display-oriented preprocessor for unsmoothed ("salt-and-pepper") statistical maps. It runs a selective median despike, an edge-preserving base smooth (guided/bilateral/gaussian), and a signal-gated unsharp pass so true clusters keep (or sharpen) their amplitude while the noise floor is denoised. `plot_overlay()` and `plot_ortho()` gained an `enhance` argument (`TRUE`, `FALSE`, or a list forwarded to `enhance_stat_map()`) to apply it inline.
* **Breaking:** `plot_overlay()` now returns a single assembled patchwork object by default (`assemble = TRUE`), honoring `ncol`/labels and suitable for `ggsave()`, with an overlay statistic colorbar (`colorbar = TRUE`, threshold marked). Pass `assemble = FALSE` for the previous behavior (draw a panel grid and return the per-slice ggplot list invisibly). `patchwork` is now a hard dependency.
* Added a publication `style = "report"` shared across `plot_overlay()`, `plot_ortho()`, and `plot_montage()` for a consistent look (dark brain tiles on a light card, bold/italic typography, a titled colorbar, brain-bounding-box cropping, and a smoothed background); `plot_overlay()` additionally adds a positive/negative legend strip. In report mode these functions return a single assembled patchwork object. The underlying features are individually toggleable via `crop` and `interpolate` (and `legend` on `plot_overlay()`).
* `plot_overlay()` now renders signed statistical maps correctly by default: overlays with both positive and negative values use a diverging palette and symmetric limits, so negatives are as visible as positives. Choosing a sequential `ov_cmap` for signed data now warns. A new `ov_alpha_mode = "ramp"` ramps opacity from the threshold to the cap, and `ov_cap` pins the magnitude scale.
* `plot_overlay()` gained `ov_alpha_mode = "soft"`: a nonlinear, self-tuning opacity curve (`alpha = clamp((|v|-lo)/(hi-lo),0,1)^gamma`). The knee `lo` defaults to the threshold, or — when unset — to the median in-mask magnitude (a robust noise-floor proxy), and `gamma` is auto-tuned from the value distribution so opacity rises rapidly away from zero while the noisy bulk stays transparent. Override the exponent with `alpha_gamma`.
* `plot_overlay()`, `plot_ortho()`, and `plot_montage()` now accept explicit numeric display ranges (e.g. `ov_range = c(-6, 6)`, `range = c(0, 1000)`) in addition to `"robust"`/`"data"`, for consistent scaling across panels and subjects. The overlay color limits and proportional-alpha denominator are now computed once over the displayed volume rather than per slice.
* `bilateral_filter()` and `bilateral_filter_4d()` now accept `range_scale` so callers can fix the intensity-kernel scale across observed and null maps instead of re-estimating it separately for each input.
* `bilateral_filter()`, internal vector bilateral filtering, and `bilateral_filter_4d()` now filter each center voxel using only in-mask, in-bounds neighbors with weight renormalization.
* `gaussian_blur()` now insulates the blur to the mask by default (new `normalize = TRUE` argument). Each in-mask output voxel is computed from in-mask neighbors only and the kernel is renormalized by the in-mask weight (a "smooth-in-mask" convolution, cf. AFNI `3dBlurInMask`). Previously the masked path read out-of-mask neighbor values into the convolution and normalized by the full kernel, so (1) out-of-mask `NaN`/`Inf` (e.g. brain-exterior values in first-level statistic maps) silently erased a `~window`-voxel shell of the masked region, and (2) finite exterior values (e.g. zero padding) biased in-mask edge voxels. Out-of-mask values --- finite or not --- can no longer affect in-mask outputs. Pass `normalize = FALSE` to restore the legacy full-kernel behavior (GitHub #22).

## Improvements

* `resolve_cmap()` now resolves any palette in `grDevices::hcl.pals()` (e.g. `"RdBu"`, `"Spectral"`, `"Reds"`) and the `"coolwarm"` diverging alias, instead of silently returning a viridis-like ramp. Unknown palette names now emit a warning rather than mis-coloring silently. Added the internal `is_diverging_cmap()` classifier.
* The 3D bilateral backend now guards zero or non-finite auto-estimated range scales, avoiding `NaN` outputs for singleton or constant masks.

## Testing

* Added tests for `enhance_stat_map()` covering impulse (salt-and-pepper) suppression, signal-peak preservation, noise-floor denoising, mask restriction, all three base methods, and the `enhance` plot arguments.
* Added tests for `resolve_cmap()` palette resolution and warnings, `is_diverging_cmap()`, numeric display ranges across `plot_overlay()`/`plot_ortho()`/`plot_montage()`, signed-map diverging defaults, `ov_alpha_mode = "ramp"`, and `assemble`/`colorbar` (including `ggsave()` round-trip).
* Added `tools/visual-qc-plots.R`, which renders labelled PNGs of the overlay scenarios (signed diverging vs. the sequential-palette bug, alpha modes, colorbar, and `enhance_stat_map()` de-speckling) over the bundled MNI template for human visual QC.
* Added regression tests for mask-normalized bilateral filtering, volume-boundary behavior, fixed `range_scale` parity with the default auto scale, singleton-mask stability, and 4D mask-normalized filtering.

# neuroim2 0.15.0

## New Features

* Added registration QC plotting helpers: `plot_checkerboard()` for alternating tiles from two registered volumes, and `plot_edge_overlay()` for comparing fixed and moving edge maps over a structural background.
* Added a dark plotting style via `theme_neuro(style = "dark")` and matching `style` arguments for montage, overlay, orthogonal, checkerboard, and edge-overlay plots.
* Added diverging colormap aliases `coldhot` and `blue-red` for signed statistical overlays.

## Improvements

* Registration QC plots now validate that all inputs share the same 3D `NeuroSpace` grid, not just matching array dimensions.
* QC panel layouts now validate `zlevels`, `along`, and `ncol` before drawing and provide clear errors for empty or invalid layouts.
* `title`, `subtitle`, and `caption` are now layout-level draw labels, preserving per-slice panel titles such as `z = 12`.
* `plot_overlay()` now validates same-grid inputs, supports `draw = FALSE`, uses consistent intensity limits across selected slices, hides repeated background legends, and supports symmetric overlay limits for signed maps.
* `plot_ortho()` now supports `draw = FALSE`, layout-level labels, dark styling, coordinate validation, and a cleaner no-legend three-plane layout.
* `scale_fill_neuro()` now squishes out-of-range values to the nearest color endpoint instead of censoring them to `NA`, preventing robust intensity limits from creating black/transparent holes in high-intensity anatomy.
* `plot_checkerboard()` now accepts `cmap` so registration checkerboards can use the same anatomical display palette as surrounding QC plots.
* `plot_edge_overlay()` now keeps all-zero edge slices transparent instead of tinting the full panel when edge limits collapse.

## Testing

* Added focused tests for registration QC plotting, including same-grid validation, invalid layout arguments, invisible `ggplot` return values, and `draw = TRUE` rendering.
* Added tests for `plot_overlay()` grid validation and no-draw return values, `plot_ortho(draw = FALSE)`, invalid orthogonal coordinates, scale squishing, and all-zero edge overlay transparency.

# neuroim2 0.14.0

## New Features

* `NeuroVecSeq` now supports matrix conversion and dense coercion helpers for easier interoperability with standard matrix workflows.
* `NeuroVec` now supports optional per-volume `volume_labels()` metadata across dense, sparse, mapped, file-backed, bigvec, and `NeuroVecSeq` backends.
* Named volume access is now supported via `vec[["label"]]`, with strict unique-match semantics, and `sub_vector()` now accepts character label vectors.
* `write_vec()` now round-trips per-volume labels through a custom NIfTI extension; `read_vec()` and low-footprint readers restore labels on load.
* Added AFNI-inspired masking helpers: `apply_mask()` to apply an existing 3D mask, `clip_level()` to estimate a foreground clip threshold or gradual clip map, and `automask()` to derive a brain-like mask from image intensities for `NeuroVol`, `NeuroVec`, sparse, mapped, and file-backed objects.

## Testing

* The test suite now runs clean: 2969 passing, no failures and no errors.
  Previously `R CMD check` reported 7 failures and 2 errors on every platform.
* `test-plot-registration-qc.R` called `ggplot2::get_labs()`, which only exists
  from ggplot2 3.5.2 while DESCRIPTION sets no version floor, so it errored on
  older ggplot2. Tests now read plot titles through a version-agnostic helper.
* The vdiffr golden-image tests have moved to `dev/visual-snapshots/`. They
  compare SVG text byte-for-byte, which encodes the font metrics of the machine
  that produced the snapshot, so a single stored image cannot match the
  Windows, macOS and Linux check matrix — and testthat deletes any `_snaps/`
  directory whose test file did not run, so they could not simply be skipped.
  `tests/testthat/test-plot-structure.R` covers the same plotting calls with
  platform-independent structural assertions (panel counts, layer counts, data
  shape, argument validation) and runs on every check. See the README in
  `dev/visual-snapshots/` for reviewing intentional visual changes.

## Documentation

* Added new introductory workflow and container vignettes and refocused the advanced volume and ROI vignettes around current package workflows.
* Standardized vignette theme/setup chunks, fixed vignette metadata and dependency declarations, and clarified return types in the `read_vec()`, `read_vol()`, and `read_image()` documentation.

## Improvements

* Fixed `NeuroVec` label inference for lists of `NeuroVol` inputs and refreshed the generated documentation for `volume_labels()` and related extractors.
* Improved 4D masking performance by moving representative-volume reduction for `clip_level()`/`automask()` into compiled code while preserving the existing `median` and `mean_abs` semantics.

## Testing

* Added differential tests for `clip_level()` against explicit voxelwise `median` and `mean_abs` reference volumes, plus a regression test for integer-valued histogram behavior.
* Added metamorphic tests for `mean_abs` sign-flip invariance and backend-parity tests covering `clip_level()`, `automask()`, and `apply_mask()` on dense, mapped, and file-backed vectors.

# neuroim2 0.12.0

## Bug Fixes

* `plot_overlay()` no longer renders the overlay flipped relative to the background. Previously the overlay was drawn via `grid::rasterGrob`, which ignored the slice's affine transform, so voxel `(1,1)` was placed at the top-left of the panel while the background placed it at its true world (mm) position. Overlays are now reoriented to match the background and span the full pixel-edge extent.
* `plot_overlay()` no longer errors with "NAs are not allowed in subscripted assignments" when the overlay contains `NA` voxels and `ov_thresh > 0`.
* `plot_overlay()` panel titles now use `x = `, `y = `, or `z = ` based on the slicing axis (`along`).

## New Features

* `NeuroVec` now supports optional per-volume `volume_labels()` metadata across dense, sparse, mapped, file-backed, bigvec, and `NeuroVecSeq` backends.
* Named volume access is now supported via `vec[["label"]]`, with strict unique-match semantics, and `sub_vector()` now accepts character label vectors.
* `write_vec()` now round-trips per-volume labels through a custom NIfTI extension; `read_vec()` and low-footprint readers restore labels on load.

## Testing

* Added focused tests for per-volume label access, concatenation semantics, `NeuroVecSeq`, and NIfTI label round-tripping.

# neuroim2 0.11.0

## Bug Fixes

* Fixed `NeuroSpace()` to derive `spacing()` and `origin()` from the affine matrix when constructed with `trans=`. Previously, `spacing()` returned `(1,1,1)` for spaces created from an explicit affine.
* Fixed `drop_dim()` for 3D-to-2D `NeuroSpace` objects to properly subset the affine matrix, preserving spatial transforms instead of lossy reconstruction from spacing/origin.

## Dependency Changes

* Replaced `crayon` with `cli` for all user-facing output. The `cli` package provides structured error messages, progress bars, and consistent ANSI formatting.
* Replaced `assertthat` with `cli::cli_abort()` across all 253 assertion sites, providing richer error messages with argument and class markup.
* Added `vdiffr` to Suggests for visual regression testing.

## Improvements

* All 28 `show()` methods now use a unified formatting style via internal `show_header()`/`show_rule()`/`show_field()` helpers.
* New `show()` methods for `DenseNeuroVol` and `NeuroSpace` (previously had no informative display).
* `random_searchlight()` and `searchlight(eager=TRUE)` now display a `cli` progress bar in interactive sessions.
* New `normalize_mask()` internal helper consolidates duplicated mask-coercion logic.

## Testing

* New `NeuroSpace` test suite (26 tests) covering construction, coordinate transforms, affine operations, and dimension manipulation.
* New NIfTI I/O round-trip tests (9 tests) verifying data/affine preservation across read-write cycles.
* New oblique affine regression tests (6 tests) for downsample, resample, and deoblique.
* New `vdiffr` plot snapshot tests (7 tests) for `plot()`, `plot_ortho()`, `plot_montage()`, and `plot_overlay()`.
* New shared test helper module with factory functions (`make_vol()`, `make_vec()`, `make_mask()`, etc.).

## Testing

* The test suite now runs clean: 2969 passing, no failures and no errors.
  Previously `R CMD check` reported 7 failures and 2 errors on every platform.
* `test-plot-registration-qc.R` called `ggplot2::get_labs()`, which only exists
  from ggplot2 3.5.2 while DESCRIPTION sets no version floor, so it errored on
  older ggplot2. Tests now read plot titles through a version-agnostic helper.
* The vdiffr golden-image tests have moved to `dev/visual-snapshots/`. They
  compare SVG text byte-for-byte, which encodes the font metrics of the machine
  that produced the snapshot, so a single stored image cannot match the
  Windows, macOS and Linux check matrix — and testthat deletes any `_snaps/`
  directory whose test file did not run, so they could not simply be skipped.
  `tests/testthat/test-plot-structure.R` covers the same plotting calls with
  platform-independent structural assertions (panel counts, layer counts, data
  shape, argument validation) and runs on every check. See the README in
  `dev/visual-snapshots/` for reviewing intentional visual changes.

## Documentation

* New "Coordinate Systems and Spatial Transforms" vignette explaining affine transforms, voxel/world coordinate conversion, orientation codes, and common gotchas.
* Consolidated `@rdname` method families, reducing man pages from 276 to 265.

# neuroim2 0.10.0

* Fixed `downsample()` for `DenseNeuroVol` and `DenseNeuroVec` so output `NeuroSpace` objects now carry a correctly rescaled affine transform. Previously, voxel dimensions could change while `trans()` still reflected the pre-downsample grid.

# neuroim2 0.9.1

* `plot_overlay()` gains an `ov_alpha_mode` argument: `"binary"` (default, existing behaviour) applies a uniform alpha to thresholded pixels, while `"proportional"` scales per-pixel alpha by the absolute overlay value for smoother blending. Internal helpers `matrix_to_colors()`, `matrix_to_rgba()`, and `matrix_to_raster_grob()` now accept an `alpha_map` argument to support this.
* Refactored orientation internals in `R/axis.R`: new helpers `.default_axcode_labels()`, `.validate_ornt()`, `ornt_transform()`, `inv_ornt_aff()`, `apply_orientation()`, `flip_axis()`, `io_orientation()`, `axcodes2ornt()`, and `ornt2axcodes()` provide a comprehensive NiBabel-compatible orientation API.
* Plot colorbar guides improved with better default labelling.
* Fixed dimension comments in the `ClusteredNeuroVec` vignette: searchlight `values()` dimensions were documented as N x T but are actually T x N (time points x neighbors).
* Removed stale `albersdown` dependency from vignette setup chunks; all vignettes now use `neuroim2::theme_neuro()` exclusively.
* Added `.ecosystem.yml` to `.Rbuildignore` to silence hidden-file NOTE.
* Added AFNI-style `deoblique()` for `NeuroSpace`/`NeuroVol`, with `gridset`/`newgrid` controls and default isotropic grid spacing equal to the minimum input voxel size. `NeuroVol` inputs are resampled to an axis-aligned deobliqued target space.
* Fixed `plot(NeuroVol)` and `plot(NeuroSlice)` memory blowups for oblique/sheared affines by rasterizing on pixel-grid coordinates instead of world-coordinate grids; added regression tests for oblique affine plotting.

# neuroim2 0.9.0

* Added 5D NIfTI support for hyper-vectors: new `read_hyper_vec()` reader returns `NeuroHyperVec`; `read_image(type = "auto")` now dispatches to `NeuroHyperVec` for 5D inputs (with optional spatial masking); and `write_vec()` now supports `NeuroHyperVec` so 5D NIfTI read/write round-trips are supported.
* Fixed `write_vec()` affine round-trip regression: NIfTI sform (direct affine) is now preferred over qform (quaternion-derived) on read, matching the convention used by FSL, FreeSurfer, and ANTs. Previously, the sform was silently replaced by the qform, causing world-coordinate drift for vector fields and warp images.
* Fixed `as_nifti_header()` to derive `qoffset` from the transform matrix translation column, ensuring internal header consistency between quaternion parameters and the sform.
* Relaxed `NeuroSpace` affine precision from 6 to 7 significant digits, matching NIfTI float32 precision and reducing cumulative round-trip truncation.

# neuroim2 0.8.7

* New public orientation utility API: `affine_to_orientation()`, `orientation_transform()`, `apply_orientation()`, `orientation_inverse_affine()`, `orientation_to_axcodes()`, `axcodes_to_orientation()`, and `affine_to_axcodes()`.
* New public space utility API: `output_aligned_space()` (with NiBabel-compatible alias `vox2out_vox()`) and `slice_to_volume_affine()` (alias `slice2volume()`), including support for `NeuroSpace`/`NeuroVol` inputs, >3D spatial handling via first 3 axes, and optional zero-based indexing compatibility.
* New public affine utility API: `apply_affine()`, `to_matvec()`, `from_matvec()`, `append_diag()`, `dot_reduce()`, `voxel_sizes()`, `obliquity()`, and `rescale_affine()`, with stricter validation and center-preserving affine rescaling semantics.
* Fixed `reorient(NeuroSpace, orient)` so returned spaces now carry the requested orientation axes (instead of preserving original axes).
* New `sub_clusters()` generic: subset a `ClusteredNeuroVol` or `ClusteredNeuroVec` by integer cluster ID, numeric, or character name (looked up in the label map). Returns a new object of the same class containing only the selected clusters.
* `scale_series()` for `DenseNeuroVec` is now ~10x faster by operating row-wise on the voxels-by-time matrix instead of double-transposing through `base::scale`.
* New dedicated `scale_series()` method for `SparseNeuroVec`: scales only masked voxels in-place on the T×K `@data` matrix, returns 0 (not NaN) for zero-variance voxels, and stays sparse.
* `as.dense()` is now an identity (no-copy) for `DenseNeuroVol`.
* Arithmetic ops (`+`, `-`, `*`, `/`, `^`) now work for `ClusteredNeuroVol` (with a warning that cluster structure is not preserved) and scalar ops for `DenseNeuroVol` and `SparseNeuroVol`.
* Logic ops (`&`, `|`) and negation (`!`) now work across `DenseNeuroVol`, `SparseNeuroVol`, and `LogicalNeuroVol`, returning `LogicalNeuroVol`.
* Compare ops (`>`, `<`, `==`, etc.) for `SparseNeuroVol` and `DenseNeuroVol` now correctly return `LogicalNeuroVol` instead of raw sparse vectors or mistyped volumes.
* Added NIfTI extension classes (`NiftiExtension`, `NiftiExtensionList`) and `read_image()` convenience wrapper.

# neuroim2 0.8.5

* Windows build fix: added `Makevars.win` to correctly link TBB libraries from RcppParallel without the non-existent `-lRcppParallel` flag.
* PDF manual fix: replaced Unicode characters (Greek letters, special symbols) in roxygen documentation with ASCII equivalents to resolve LaTeX errors during PDF generation.

# neuroim2 0.8.4

* I/O reliability: enforce full-length binary reads (detects truncated images) and ensure gz connections are cleaned up on error.
* Data correctness: apply per-volume slope+intercept scaling consistently across NeuroVol/NeuroVec/SparseNeuroVec loaders (treat NIfTI slope==0 as identity).
* Performance: SparseNeuroVecSource no longer materializes full 4D arrays; reads masked voxels via mmap for uncompressed files or streams volumes sequentially for gz files.
* Docs: fix missing Rd entries/aliases for several exported functions and S4 methods; remove vignette dependency on albersdown.

# neuroim2 0.8.3

* Arithmetic/comparison fixes: SparseNeuroVec now unions masks and keeps outputs sparse; sparse–sparse NeuroVol arithmetic returns SparseNeuroVol; numeric vs SparseNeuroVol comparison no longer errors.
* 4D bilateral filter now measures intensity variance across all timepoints in the mask and skips non-finite neighbours, eliminating spurious NaNs on constant or noisy inputs.
* Added regression tests covering the zero-window identity and constant-volume stability for the parallel 4D bilateral filter backend.
* New `meta_info()` helper returns a normalized list of basic header metadata from a filename or `FileMetaInfo` (dim, spacing, origin, trans, path, etc.), making 3D/4D image introspection simpler for new users.

# neuroim2 0.8.2

* README refreshed: CRAN/R-universe install, CI/coverage badges, website & cheatsheet links.
* Docs: `spherical_roi()` now cross-links to `spherical_roi_set()`; ROI vignette shows multi-ROI creation.
* SparseNeuroVec:
  - New validity checks to catch data/mask/space shape mismatches (#5).
  - Robust `as.matrix.SparseNeuroVec()` implementation (#2).
* New `resample_to()` wrapper for readable interpolation names; delegates to existing `resample()` methods.


# neuroim2 0.8.1

* Initial CRAN submission.

# neuroim2 0.17.0

## Bug Fixes

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

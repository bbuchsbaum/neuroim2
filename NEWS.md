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

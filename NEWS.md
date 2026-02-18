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

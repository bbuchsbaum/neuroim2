# neuroim2 0.8.2 (development)

* README refreshed: CRAN/R-universe install, CI/coverage badges, website & cheatsheet links.
* Docs: `spherical_roi()` now cross-links to `spherical_roi_set()`; ROI vignette shows multi-ROI creation.
* SparseNeuroVec:
  - New validity checks to catch data/mask/space shape mismatches (#5).
  - Robust `as.matrix.SparseNeuroVec()` implementation (#2).
* New `resample_to()` wrapper for readable interpolation names; delegates to existing `resample()` methods.


# neuroim2 0.8.1

* Initial CRAN submission.

# neuroim2 vs. nilearn harness

nibabel is the reference for file I/O and does no image processing; nilearn is
the reference for the processing that neuroim2 and nilearn both do. This
harness covers only that overlap. It deliberately excludes nilearn's machine
learning (decoders, connectivity estimators, parcellation learners), which
neuroim2 does not attempt, and neuroim2's own out-of-core representations,
which nilearn has no equivalent for.

Findings are written up in `dev/nilearn-gap-analysis.md`.

## What is compared

| Area | neuroim2 | nilearn / scipy |
|---|---|---|
| Gaussian smoothing | `gaussian_blur()` | `image.smooth_img()` |
| Resampling | `resample()`, `resample_to()` | `image.resample_img()`, `resample_to_img()` |
| Canonical reorientation | `as_canonical()` | `image.reorder_img()` |
| Connected components | `conn_comp()` | `scipy.ndimage.label()` |
| Brain masking from data | `automask()` | `masking.compute_epi_mask()` |
| Spherical ROI signals | `spherical_roi()` + `series()` | `maskers.NiftiSpheresMasker` |
| Parcel signals | `ClusteredNeuroVec()` | `maskers.NiftiLabelsMasker` |
| Masked extraction | `series(vec, idx)` | `maskers.NiftiMasker` |
| Searchlight neighbourhoods | `searchlight_coords()` | `sklearn.neighbors` graph, as `decoding.SearchLight` builds it |
| Utilities | `mean()`, `[[`, `>`, `as.mask()`, `concat()` | `mean_img`, `index_img`, `threshold_img`, `binarize_img`, `concat_imgs` |

Smoothing is measured by its **point-spread function** — smooth an impulse,
measure the FWHM that comes out — rather than by comparing the arguments that
went in. That is the only way to see whether a call delivers the smoothing it
was asked for.

## Requirements

Python with `nilearn`, `nibabel`, `numpy`, `scipy` and `scikit-learn`, and an R
with neuroim2 installed **and a real `RNiftyReg`**, which the resampling and
reorientation probes need.

## Running

```sh
cd dev/bench/nilearn
export NEUROIM2_NL_DIR=/tmp/nl
python3 make_data.py      # shared inputs; both sides read the same files
python3 ref_py.py         # nilearn/scipy results + timings -> ref.json
Rscript  probe_r.R        # neuroim2 results + timings      -> res.json
python3 compare.py        # the report
```

Add `NEUROIM2_NL_SCALING=1` to both `ref_py.py` and `probe_r.R` for the
connected-component scaling probe. It runs `conn_comp()` at three volume sizes
up to 182x218x182, which currently takes about **15 minutes** at the largest —
that slowness is the finding, so it is opt-in rather than part of the default
run.

## Reading the numbers

Timings are medians of 1–3 repetitions after a warm-up, on one machine; the
ratios are the finding, not the milliseconds. Two comparisons are not
apples-to-apples unless you read the note next to them:

- **Smoothing** is timed twice on the neuroim2 side. `window = 2` is the call
  people write, but it evaluates a narrower kernel than nilearn does, so it is
  doing less work; `window = 4` is the matched-kernel comparison.
- **`conn_comp()`** returns cluster sizes and per-cluster voxel lists as well as
  the labelling, where `scipy.ndimage.label()` returns only the labelling. That
  accounts for some of the difference, but not for three orders of magnitude —
  the labelling itself is where the time goes.

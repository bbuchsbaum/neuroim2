# neuroim2 vs. nibabel harness

Two independent comparisons against `nibabel`, which is the de-facto reference
implementation for neuroimaging file I/O:

- **conformance** — can neuroim2 read what nibabel can read, and does it get the
  same numbers and the same affine?
- **timing** — how long do the equivalent read/write calls take?

Findings from the first run are written up in `dev/nibabel-gap-analysis.md`.

## Requirements

Python with `nibabel` and `numpy` (plus `scipy` for the smoothing comparison),
and an R with neuroim2 installed. The two sides never talk to each other
directly; they exchange JSON files.

## Conformance

```sh
cd dev/bench/nibabel
python3 make_probes.py                          # probes/ + probes/expected.json
Rscript  probe_r.R > probes/r_results.json
python3  compare.py                             # exits non-zero on any mismatch
```

`make_probes.py` writes 30 NIfTI files that each isolate one part of the
standard: every NIfTI-1 datatype code, the `scl_slope`/`scl_inter` conventions
(including the `0` and `NaN` "no scaling" spellings), qform-only / sform-only /
neither, `qfac = -1`, big-endian, `.hdr`/`.img` pairs, ANALYZE 7.5, NIfTI-2,
gzip, non-finite payloads, 5-D, a trailing singleton dimension, and a header
extension. Several are produced by patching raw header bytes after the fact,
because nibabel's writer normalises the values we want to probe.

`compare.py` compares only order-invariant summaries (min/max/sum, NaN and Inf
counts) plus element `[0]` and the affine. nibabel flattens C-order and R
flattens Fortran-order, so "the second element" is a different voxel in the two
languages and is not comparable.

## Timing

```sh
export NEUROIM2_BENCH_DIR=/tmp/bench
python3 make_bench_data.py
Rscript  bench_r.R
python3  bench_py.py
```

Both scripts report the median of 3–5 repetitions after a warm-up call, on a
182x218x182 3-D volume and a 64x64x36x200 4-D run. Run them back to back so the
page cache is in the same state for both.

Two things to keep in mind when reading the numbers. nibabel's `get_fdata()` is
the honest "read the whole image" comparison — `np.asanyarray(img.dataobj)`
returns a memmap without touching the data, so its near-zero timing is not a
read. And neuroim2 always materialises `double`, while nibabel keeps the native
dtype, so for int16 input neuroim2 moves 4x the bytes by construction.

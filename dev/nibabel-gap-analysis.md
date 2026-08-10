# Where neuroim2 stands against nibabel

Everything below was measured in one container against `nibabel` 5.4.2 / numpy
2.4.6 / scipy 1.17.1 and neuroim2 at `2c6c782`, built from source under R 4.3.3.
The harness is in `dev/bench/nibabel/` and re-runs from scratch. Absolute times
are specific to this machine; the ratios are what transfer.

---

## First, what the comparison actually means

nibabel and neuroim2 do not overlap as much as the goal implies, and pretending
otherwise would produce a bad roadmap. nibabel is deliberately an I/O and
header-semantics library: it reads and writes about a dozen formats, exposes the
affine and the header faithfully, and refuses to do image processing. neuroim2
reads one and a half formats and spends its 26k lines of R and 4.4k lines of C++
on things nibabel has no opinion about — ROIs, searchlights, sparse and
file-backed representations, filtering, clustering, plotting.

So "better than nibabel" splits cleanly into two claims with very different
difficulty:

- **On I/O, nibabel is the reference implementation and neuroim2 currently
  loses — by a lot, and for reasons that are specific and fixable.** This is the
  part of the goal that is real work, and it is where this document spends its
  time. Nothing here requires new science; it requires a compiled reader, a
  wider datatype table, and a header that survives a round trip.
- **On analysis, there is nothing to beat**, because nibabel does not compete.
  The relevant yardsticks there are nilearn and scipy.ndimage, and the existing
  `dev/optimization-roadmap.md` is already the right document for that axis.

The useful target is therefore: *match nibabel on read/write speed and NIfTI
conformance, keep the analysis layer that nibabel does not have, and be the
better choice for out-of-core work* — which, as measured below, neuroim2
already is.

---

## Performance

Median of 3–5 repetitions after a warm-up, warm page cache, both languages run
back to back. Volume: 182x218x182. Run: 64x64x36x200 (56 MB int16).

| Operation | neuroim2 | nibabel | |
|---|---:|---:|---|
| read header | 2 ms | <1 ms | ~5x slower |
| read 3-D int16 `.nii` (14 MB) | 0.46 s | 0.024 s | **19x slower** |
| read 3-D float32 `.nii` (28 MB) | 0.28 s | 0.025 s | **11x slower** |
| read 3-D float32 `.nii.gz` | 0.45 s | 0.21 s | 2.1x slower |
| read 4-D int16 `.nii` (56 MB) | 3.09 s | 0.094 s | **33x slower** |
| read 4-D int16 `.nii.gz` | 2.81 s | 0.53 s | 5.3x slower |
| read one sub-volume out of a 4-D file | 8 ms | 1 ms | 8x slower |
| write 3-D float32 `.nii` | 0.22 s | 0.014 s | **16x slower** |
| write 3-D float32 `.nii.gz` | 1.23 s | 0.93 s | 1.3x slower |
| Gaussian smooth, FWHM 6 mm, 3-D | 0.60 s | 0.12 s | 5x slower |
| Gaussian smooth, FWHM 6 mm, 4-D (200 vols) | 2.04 s | 0.48 s | 4.3x slower |

The nibabel column uses `get_fdata()`, which is the honest full-read comparison.
`np.asanyarray(img.dataobj)` returns a memmap without touching the data and
times at 0.000 s; that is not a read and is not in the table. The smoothing row
is generous to neuroim2 — scipy's default `truncate=4.0` evaluates an 11-tap
kernel where `window = 2` evaluates 5, so nibabel is doing more work and is
still 5x faster.

### Where the 4-D read time goes

`read_mapped_vols()` (`R/binary_io.R:176`) is the whole story. To read a
contiguous range of a memory-mapped file it does this:

```r
idx_set <- as.vector(outer(seq_len(nels), (idx - 1L) * nels, "+"))
ret     <- .read_mmap(meta, idx_set)
mat     <- matrix(ret, nels, length(idx))
t(mat)                                   # -> [time x voxels]
```

and `DenseNeuroVec()` (`R/neurovec.R`) then sees a T x V matrix, finds
`ncol(data) == prod(dim(space)[1:3])`, and calls `t(data)` to put it back the
way it arrived. Measured on the reference run, out of 3.09 s total:

```
build the index vector      1.91 s     29.5 M doubles = 236 MB, to enumerate 1..N
scattered mmap gather       1.05 s     element-by-element through R's `[`
transpose  V x T -> T x V   0.85 s
transpose back in the ctor  0.38 s
```

Three of those four steps are pure loss, and the fourth should be sequential.
Replacing the whole thing with one `readBin` of `nels * nt` elements — still
pure R — gives a **bit-identical** object (`identical()` on both the data and
the `NeuroSpace`) in 1.02 s, a 3.6x speed-up for a change of a few lines.

### Where the rest of the time goes: `readBin`

The 3-D path (`R/neurovol.R:367`) has none of the above pathology — it is
already a single sequential read — and is still 19x slower than nibabel. That
one is `readBin` itself. Reading 7.2 M int16 into an R integer vector costs
0.175 s; a compiled `fread` loop that converts straight to `double` in one pass
costs 0.019 s and returns an `identical()` vector. R's connection layer is
roughly **9x** slower than the obvious C, and every read path in the package
goes through it.

Prototyped headroom, all verified `identical()` against current output:

| | now | prototype | |
|---|---:|---:|---|
| 3-D int16 payload | 0.46 s | 0.019 s | 24x — would beat nibabel's 0.024 s |
| 4-D payload, compiled | 3.09 s | 0.17 s | 18x — within 1.8x of nibabel |
| 4-D, sequential `readBin` only (pure R) | 3.09 s | 1.02 s | 3.6x, ~10 lines |

The conclusion worth acting on: **the I/O gap is not architectural.** A compiled
reader and writer put neuroim2 at parity with nibabel on uncompressed NIfTI, and
one of the two steps to get there is a few lines of R.

### A structural disadvantage worth naming

neuroim2 always materialises `double`. nibabel keeps the native dtype. For int16
input neuroim2 therefore moves and holds 4x the bytes, by construction, before
any inefficiency. A 56 MB run becomes 236 MB in memory. R's only cheaper numeric
type is 32-bit `integer`, whose `NA` collides with `INT_MIN` (see the
conformance findings), so this is a genuine ceiling rather than an oversight —
but it argues for making the sparse, masked and file-backed paths the
recommended way to touch large data, not the fallback.

### Where neuroim2 already wins

Scattered time-series extraction from a 4-D file that is not in memory —
the actual inner loop of MVPA, connectivity and any ROI analysis:

| 5,000 scattered voxels x 200 timepoints | |
|---|---:|
| neuroim2, `read_vec(mask=)` | **0.044 s** |
| neuroim2, mmap mode, open + `series()` | **0.144 s** |
| nibabel, full load then fancy index | 0.011 s |
| nibabel, `dataobj[i,j,k,:]` per voxel (out-of-core) | 2.99 s |

nibabel wins the middle row only by loading the entire image first, which stops
being an option at the sizes where this matters. For genuinely out-of-core
scattered access — nibabel's `ArrayProxy` re-entering `fileslice` per voxel —
neuroim2 is **21x faster**. This is a real, defensible advantage and it is worth
building on rather than treating as incidental.

---

## Conformance

`dev/bench/nibabel/make_probes.py` writes 30 NIfTI files that each isolate one
part of the standard, `probe_r.R` reads them through neuroim2, and `compare.py`
diffs the result against what nibabel gets. **16 of 30 agree today.**

What passes is worth stating first, because it is the load-bearing part:
`scl_slope`/`scl_inter` scaling in all its spellings including `0` and negative
slopes, qform-only and sform-only files, `qfac = -1`, big-endian, gzip,
`.hdr`/`.img` pairs, header extensions, NaN/Inf payloads, and float32/float64/
int16/uint8 data. The affine handling in particular is correct, including the
sform-over-qform precedence that FSL, FreeSurfer and ANTs use.

### Format and datatype coverage

| Gap | Detail |
|---|---|
| **7 of 12 NIfTI-1 datatypes rejected** | `int8` (256), `uint16` (512), `uint32` (768), `int64` (1024), `uint64` (1280), `complex64` (32), `RGB24` (128) all fail with *"Unsupported NIfTI data-type code"*. nibabel reads every one. `uint16` is not exotic — CT and several scanner exports use it. |
| **NIfTI-2 unsupported** | `read_header()` fails at *"header size is not 348"*. NIfTI-2 is how images past the 32k-per-axis limit are stored, and it is what CIFTI-2 is built on. |
| **ANALYZE 7.5 affine is wrong** | The origin field is ignored, so the affine comes back as `diag(2,2,2)` with a zero translation where nibabel returns the correct one. Legacy, but silently wrong is worse than unsupported. |
| Formats absent entirely | MGH/MGZ, GIFTI, CIFTI-2, FreeSurfer geometry/annot/label, PAR/REC, MINC, ECAT, TRK/TCK. AFNI BRIK/HEAD is read-only in both libraries. NIML reading is neuroim2-only; nibabel has no equivalent. |
| Write is NIfTI-1 or HDF5 only | `write_vol`/`write_vec` reject every other `format=` string, AFNI included. Same practical position as nibabel, which also only writes what it fully models. |

### Correctness

| Gap | Detail |
|---|---|
| **`scl_slope = NaN` is a hard error** | *"Invalid scale slope for volume 1: NaN"*. NaN is the documented "no scaling" spelling alongside `0` and is what nibabel's own in-memory header carries by default; files that use it exist and neuroim2 cannot open them. |
| **int32 loses `INT_MIN`** | A voxel holding `-2147483648` comes back as `NA`, because that bit pattern *is* R's `NA_integer_`. It then drops out of `min`, `sum` and every downstream statistic silently. Reading int32 as `double` rather than `integer` fixes it. |
| **5-D with a singleton 4th axis is collapsed** | A `(6,7,5,1,3)` image — the shape of a per-voxel vector field, `intent_code` 1007 — warns and becomes `(6,7,5,3)`, i.e. three timepoints. The warning is `read_nifti_header:` and offers the caller no way to opt out. |
| Trailing singleton dropped | `(6,7,5,1)` loads as `(6,7,5)`. A 1-volume run and a 3-D volume become indistinguishable. |
| `qform_code = 0` **and** `sform_code = 0` | neuroim2 uses the qform anyway; the spec says fall back to METHOD 1 (pixdim, origin at the centre), which is what nibabel does. |
| Truncated file | `argument is of length zero`, from somewhere internal. nibabel says *"Expected 1008 bytes, got 808 bytes from <path> - could the file be damaged?"*. |

### Write fidelity

Read a nibabel-written file, write it back with `write_vec()`, diff the two
headers with nibabel. Preserved: `dim`, `pixdim[1:3]`, the affine, and every
data value. Discarded:

```
pixdim[4]  (TR)      2.0  -> 0.0     <- fMRI timing, silently gone
xyzt_units          mm+s  -> mm
sform_code    4 (MNI152)  -> 1 (scanner_anat)   <- relabels the coordinate space
datatype        int16 (4) -> float32 (16)       <- file doubles in size
descrip, aux_file, intent_code, intent_name,
cal_min, cal_max, slice_code, slice_end,
slice_duration, toffset                 -> all zeroed
```

The first three are the serious ones. Losing the TR breaks any downstream tool
that reads timing from the header. Rewriting `sform_code` from 4 to 1 tells the
next tool that an MNI-space image is in scanner coordinates, which is a wrong
answer rather than a missing one. The cause is structural, not a list of
oversights: `as_nifti_header()` (`R/nifti_io.R:155`) builds a header from
`createNIfTIHeader()` defaults plus the object's geometry, and never consults
the header the image was read from — `pixdim[5:8]` is hardcoded to `0`,
`xyzt_units` to `2`, `qform_code`/`sform_code` to `1`.

### Integer output truncates instead of scaling

`write_vol(v, f, data_type = "SHORT")` on data spanning ±3.7 goes through
`writeBin(as.integer(els), ...)` (`R/binary_io.R:425`), which truncates toward
zero. Round-trip error: **0.998**. No warning. nibabel writing the same data as
int16 computes `scl_slope = 1.13e-4`, `scl_inter = 5.6e-5` to fill the range,
and round-trips at **5.4e-5** — four orders of magnitude better.

---

## Quality signals that are already good

Worth recording, because the gap list above is not the whole picture.

- **The test suite is healthy.** 4,233 passing assertions, 0 failures. The 12
  errors in this container are all traceable to the dependency stubs it needed
  (`RNiftyReg`, `bigstatsr`) and an old `ggplot2`, not to the package.
- **The affine layer is correct** on every probe that exercises it — qform,
  sform, precedence, `qfac`, quaternion round-trip.
- **The optimisation work in `dev/optimization-roadmap.md` is real** and its
  process (golden sweeps, adversarial review, permanent equivalence tests) is
  stronger than what most R packages do. The three landed phases hold up.
- **NIML reading** has no nibabel equivalent, and the HDF5 (`format = "H5"`)
  write path is a capability nibabel does not offer at all.
- The out-of-core story is genuinely ahead of nibabel's, as measured above.

---

## What to do, in order

Sequenced by (gap size) x (confidence the fix is contained). Items 1–4 are the
ones that change the answer to "is this better than nibabel at I/O".

**1. Sequential read for contiguous ranges.** `read_mapped_vols()` stops
building an index vector to enumerate `1..N`, and `DenseNeuroVec()` stops
transposing a matrix it was just handed transposed. 3.6x on every 4-D read, in
pure R, bit-identical output. This is the cheapest real win in the package
right now.

**2. Compiled reader and writer.** One `.cpp` translation unit: file offset,
element count, source dtype, endianness, signedness → `REALSXP`, and the
inverse. Retires `readBin`/`writeBin` from the hot path and subsumes item 1.
Expected ~15–20x on uncompressed reads and parity with nibabel; ~15x on writes.
Also the natural place to fix the `int32`/`INT_MIN` collision (read to `double`,
never through R `integer`) and to add the missing datatypes, since the dtype
switch lives in exactly one place.

**3. Full datatype coverage.** `int8`, `uint16`, `uint32`, `int64`, `uint64`
straight away — they are table entries once item 2 exists. `complex64/128` and
`RGB24/RGBA32` need a representation decision first (complex `NeuroVol`? a
3-channel one?) and should not block the rest.

**4. Header preservation on write.** Carry the source `FileMetaInfo` on the
object and have `as_nifti_header()` start from it instead of from
`createNIfTIHeader()` defaults, overriding only what the geometry actually
changed. Then add auto-scaling for integer output, and default `write_vol` to
the source datatype rather than `FLOAT`. Make the round-trip test the golden
one: read → write → read, header equal field by field. This is the single
change that most affects whether neuroim2 is safe in a pipeline with other
tools.

**5. NIfTI-2.** Contained once item 2 exists — a second header struct and a
version branch in `read_header()`. Needed for large images and a prerequisite
for CIFTI later.

**6. Diagnostics.** A truncated file should say how many bytes were expected and
how many arrived. A 5-D image should not be silently reinterpreted. `qform_code
= sform_code = 0` should follow METHOD 1 and say so. These are small and they
are the difference between a library people trust and one they work around.

**7. Native 4-D filtering.** `gaussian_blur` on a `NeuroVec` currently means an
R-level loop over volumes. nibabel's `smooth_image` handles 4-D directly and is
4.3x faster.

**8. Lean into the out-of-core advantage.** The 21x win on scattered
out-of-core access is the one place neuroim2 is not chasing. Making the masked
and mmap paths the documented default for large data — rather than an option
behind `mode=` — is a positioning change more than an engineering one, and it
is what the memory-footprint ceiling argues for anyway.

Items 1–2 together are what move the headline numbers; item 4 is what moves the
correctness argument. If only one thing gets done, do item 2 and fold item 1
into it.

## Reproducing

```sh
cd dev/bench/nibabel
python3 make_probes.py && Rscript probe_r.R > probes/r_results.json && python3 compare.py
NEUROIM2_BENCH_DIR=/tmp/bench python3 make_bench_data.py
NEUROIM2_BENCH_DIR=/tmp/bench Rscript bench_r.R
NEUROIM2_BENCH_DIR=/tmp/bench python3 bench_py.py
```

Caveats that apply to every number above. One machine, one run each, medians of
3–5 repetitions — the ratios are the finding, not the milliseconds. CRAN is
unreachable from this container, so `RNiftyReg` and `bigstatsr` were replaced
with error-on-call stubs; `RNifti` itself is real, which matters because it is
the only dependency on an I/O path. No profiling was done on the AFNI or NIML
readers. The prototypes in the performance section were measured standalone and
checked with `identical()` against current output, but none of them has been
through the golden sweep in `dev/bench/golden.R`, so they are evidence that the
gap is closable, not finished work.

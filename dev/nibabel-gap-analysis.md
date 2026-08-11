# Where neuroim2 stands against nibabel

Measured in one container against `nibabel` 5.4.2 / numpy 2.4.6 / scipy 1.17.1,
with neuroim2 built from source under R 4.3.3. The harness is in
`dev/bench/nibabel/` and re-runs from scratch. Absolute times are specific to
this machine; the ratios are what transfer.

**Status.** The eight items this document originally raised have been
implemented; see `NEWS.md` for what changed. The tables below are the current
state, with the numbers from the first pass kept alongside so the change is
visible. What remains is listed under *Still open*.

---

## First, what the comparison actually means

nibabel and neuroim2 do not overlap as much as the goal implies, and pretending
otherwise would produce a bad roadmap. nibabel is deliberately an I/O and
header-semantics library: it reads and writes about a dozen formats, exposes the
affine and the header faithfully, and refuses to do image processing. neuroim2
reads one and a half formats and spends its 26k lines of R and 4.7k lines of C++
on things nibabel has no opinion about — ROIs, searchlights, sparse and
file-backed representations, filtering, clustering, plotting.

So "better than nibabel" splits cleanly into two claims with very different
difficulty:

- **On I/O, nibabel is the reference implementation.** neuroim2 used to lose
  badly here, for reasons that turned out to be specific and fixable.
- **On analysis, there is nothing to beat**, because nibabel does not compete.
  The relevant yardsticks there are nilearn and scipy.ndimage, and
  `dev/optimization-roadmap.md` is the right document for that axis.

The useful target is therefore: *match nibabel on read/write speed and NIfTI
conformance, keep the analysis layer that nibabel does not have, and be the
better choice for out-of-core work.*

---

## Performance

Median of 3–5 repetitions after a warm-up, warm page cache, both languages run
back to back. Volume: 182x218x182. Run: 64x64x36x200 (56 MB int16).

| Operation | was | now | nibabel | |
|---|---:|---:|---:|---|
| read header | 2 ms | 2 ms | <1 ms | |
| read 3-D int16 `.nii` (14 MB) | 0.46 s | **0.049 s** | 0.027 s | 1.8x slower |
| read 3-D float32 `.nii` (28 MB) | 0.28 s | **0.053 s** | 0.028 s | 1.9x slower |
| read 3-D float32 `.nii.gz` | 0.45 s | **0.220 s** | 0.247 s | **parity** |
| read 4-D int16 `.nii` (56 MB) | 3.09 s | **0.371 s** | 0.107 s | 3.5x slower |
| read 4-D int16 `.nii.gz` | 2.81 s | **0.750 s** | 0.585 s | 1.3x slower |
| read one sub-volume out of a 4-D file | 8 ms | **3 ms** | ~1 ms | |
| write 3-D float32 `.nii` | 0.22 s | 0.20 s | 0.07 s | ~3x slower |
| write 3-D float32 `.nii.gz` | 1.23 s | 1.20 s | 1.05 s | 1.1x slower |
| Gaussian smooth, FWHM 6 mm, 3-D | 0.60 s | **0.21 s** | 0.13 s | 1.6x slower |
| Gaussian smooth, FWHM 6 mm, 4-D (200 vols) | 2.04 s | **0.97 s** | 0.57 s | 1.7x slower |

The nibabel column uses `get_fdata()`, which is the honest full-read comparison.
`np.asanyarray(img.dataobj)` returns a memmap without touching the data and
times at 0.000 s; that is not a read and is not in the table.

Two caveats on specific rows. The **write** row is the least stable measurement
here: it is dominated by creating a fresh file, and both libraries' numbers move
by 3–5x with the state of the container's filesystem — nibabel measured
anywhere between 0.016 s and 0.090 s across runs, neuroim2 between 0.05 s and
0.30 s. Writing repeatedly to the *same* path, where file creation is not in the
measurement, the compiled writer moves 28 MB in 0.026 s against numpy's 0.017 s.
Treat "~3x" as the honest reading and the table entry as illustrative. The
**smoothing** row is generous to neuroim2 in the other direction: scipy's
default `truncate = 4.0` evaluates an 11-tap kernel where `window = 2` evaluates
5, so nibabel is doing more work per voxel.

### What the read path used to do

`read_mapped_vols()` was the whole 4-D story. To read a contiguous range of a
memory-mapped file it built a 236 MB vector of `double` indices with `outer()`
simply to enumerate `1..N`, gathered them one at a time through R's `[`,
materialised a time-by-voxel matrix, and handed it to `DenseNeuroVec()`, which
transposed it straight back. Measured, out of 3.09 s: 1.91 s building the index
vector, 1.05 s on the scattered gather, and 1.23 s on the two transposes.

The 3-D path had none of that pathology — it was already a single sequential
read — and was still 19x slower than nibabel. That one was `readBin` itself:
reading 7.2 M int16 into an R integer vector cost 0.175 s where a compiled
`fread`-and-convert loop cost 0.019 s. R's connection layer is roughly 9x slower
than the obvious C, and every read path went through it.

Both are gone. `src/nifti_data_io.cpp` converts between the stored type and
`double` in one pass, for plain and gzipped files alike, and the read path
adopts the resulting vector as the object's payload in place rather than letting
`new()` copy it.

### A structural disadvantage worth naming

neuroim2 always materialises `double`. nibabel keeps the native dtype. For int16
input neuroim2 therefore moves and holds 4x the bytes, by construction, before
any inefficiency: a 56 MB run becomes 236 MB in memory. R's only cheaper numeric
type is 32-bit `integer`, whose `NA` collides with `INT_MIN`, so this is a
genuine ceiling rather than an oversight. It is why the masked and file-backed
paths are now the documented recommendation for large data rather than a
fallback.

### Where neuroim2 wins

Scattered time-series extraction from a 4-D file that is not in memory — the
actual inner loop of MVPA, connectivity and any ROI analysis:

| 5,000 scattered voxels x 200 timepoints | |
|---|---:|
| neuroim2, `read_vec(mask=)` | **0.044 s** |
| neuroim2, mmap mode, open + `series()` | **0.144 s** |
| nibabel, full load then fancy index | 0.011 s |
| nibabel, `dataobj[i,j,k,:]` per voxel (out-of-core) | 2.99 s |

nibabel wins the middle row only by loading the entire image first, which stops
being an option at the sizes where this matters. For genuinely out-of-core
scattered access — nibabel's `ArrayProxy` re-entering `fileslice` per voxel —
neuroim2 is **21x faster**.

---

## Conformance

`dev/bench/nibabel/make_probes.py` writes 30 NIfTI files that each isolate one
part of the standard, `probe_r.R` reads them through neuroim2, and `compare.py`
diffs the result against what nibabel gets.

**28 of 30 agree; the other 2 are refused by design. It was 16 of 30.**

What passes: every fixed-width datatype code (`uint8`, `int8`, `int16`,
`uint16`, `int32`, `uint32`, `int64`, `uint64`, `float32`, `float64`), the
`scl_slope`/`scl_inter` conventions in all their spellings including `0` and
`NaN`, qform-only and sform-only files, `qfac = -1`, both `*_code` fields zero,
big-endian, gzip, `.hdr`/`.img` pairs, ANALYZE 7.5, NIfTI-2, header extensions,
non-finite payloads, 5-D, and a trailing singleton dimension.

The two refusals are complex (`DT_COMPLEX64`) and colour (`DT_RGB24`) data. A
`NeuroVol` holds one real number per voxel; inventing a projection — a modulus,
a luminance — would answer a question the caller did not ask. The error names
the type and says so, which is the deliberate outcome rather than a gap.

### Write fidelity

Read a nibabel-written file with a fully populated header, write it back with
`write_vec()`, diff the two headers with nibabel: **every field matches, and the
image payload is byte-identical on disk.** Previously the round trip preserved
only `dim`, `pixdim[1:3]`, the affine and the values, and discarded the TR,
`xyzt_units`, `sform_code` (rewriting an MNI-space image as scanner-anatomical),
`descrip`, `aux_file`, the intent fields, `cal_min`/`cal_max`, the slice-timing
fields and `toffset`, while promoting int16 to float32.

Integer output: writing data spanning ±3.7 as `SHORT` had a round-trip error of
**0.998**, because `writeBin(as.integer(x))` truncates toward zero without a
warning. It is now **5.4e-05** — the same precision nibabel achieves, by the same
means, deriving `scl_slope`/`scl_inter` to fit the target range. Data that is
already integral and in range is still written unscaled, so masks and label
volumes stay exact.

---

## Quality signals

- **The test suite is healthy.** 4,447 passing assertions, 0 failures — 4,233
  before, plus 214 in the new `test-io-conformance.R`. The 12 errors in this
  container are all traceable to the dependency stubs it needs (`RNiftyReg`,
  `bigstatsr`) and an old `ggplot2`, not to the package.
- **The affine layer is correct** on every probe that exercises it — qform,
  sform, precedence, `qfac`, quaternion round-trip, and now the METHOD 1
  fallback.
- **NIML reading** has no nibabel equivalent, and the HDF5 (`format = "H5"`)
  write path is a capability nibabel does not offer at all.
- The out-of-core story is genuinely ahead of nibabel's, as measured above.

---

## Still open

Ordered by what would matter most next.

**1. The 4-D read is 3.5x off nibabel, where the 3-D read is 1.8x.** The
compiled gather itself is 0.163 s of the 0.371 s; the rest is spread across
object construction, `NeuroSpace`, volume labels and dispatch. Worth a profile
before guessing.

**2. Formats nibabel has and neuroim2 does not**: MGH/MGZ, GIFTI, CIFTI-2,
FreeSurfer geometry/annot/label, PAR/REC, MINC, ECAT, TRK/TCK. NIfTI-2 landing
makes CIFTI-2 tractable, since it is a NIfTI-2 container with an XML extension —
and the extension machinery already exists. MGH/MGZ is the cheapest of the rest
and the most used.

**3. Complex and colour data.** Currently a clear refusal. Supporting them means
deciding what a complex or 3-channel `NeuroVol` *is*, which is a class-design
question, not an I/O one.

**4. Writing `.hdr`/`.img` pairs.** Reading them works; `as_nifti_header()` has
an `oneFile` argument but no exported route to it.

**5. The header does not survive derivation.** `read_vol(f)` and `write_vol()`
round-trip faithfully, but `v * 2` returns an object with an empty `header`
slot, because the ~20 `Arith`/`Compare` methods construct fresh objects. nibabel
carries `img.header` through. Whether a derived image *should* claim the
source's `descrip` and `cal_min` is a real question; the current answer is
conservative and undocumented, which is the part to fix.

**6. `enhance_stat_map()` and the bilateral filters** were not examined; only
`gaussian_blur` was profiled.

**7. Threading.** The 4-D filter loop is embarrassingly parallel and RcppParallel
is already a dependency. Untouched, because doing it safely needs a C++ driver
that allocates nothing from R inside the parallel region.

## Reproducing

```sh
cd dev/bench/nibabel
python3 make_probes.py && Rscript probe_r.R > probes/r_results.json && python3 compare.py
NEUROIM2_BENCH_DIR=/tmp/bench python3 make_bench_data.py
NEUROIM2_BENCH_DIR=/tmp/bench Rscript bench_r.R
NEUROIM2_BENCH_DIR=/tmp/bench python3 bench_py.py
```

Caveats that apply to every number above. One machine, one run each, medians of
3–5 repetitions — the ratios are the finding, not the milliseconds, and the
write rows are volatile for the reason given above. CRAN is unreachable from
this container, so `RNiftyReg` and `bigstatsr` were replaced with error-on-call
stubs; `RNifti` itself is real, though nothing on the I/O path uses it any more.
The AFNI and NIML readers were not profiled or probed.

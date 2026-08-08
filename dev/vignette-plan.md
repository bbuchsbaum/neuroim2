# Vignette rewrite plan

A proposal to take the neuroim2 article set from 14 overlapping documents to 9
deliberate ones, in an order that matches how people actually learn the package.

Status: proposal. Nothing in `vignettes/` has been changed.

---

## 1. Diagnosis

The current set is 14 files, ~3,200 lines. The prose is mostly clear and the
coverage is broad. What holds it back is not writing quality — it is five
structural problems, listed worst first.

### 1.1 The demo data contradicts what the text says it shows

`inst/extdata/global_mask_v4.nii` is the "4D fMRI time series" in nine of the
fourteen articles. Its header and voxel values:

| file | dim | datatype | distinct values |
|:--|:--|:--|:--|
| `global_mask_v4.nii` | 64 x 64 x 25 x **4** | float32 | **`{0, 1}`** |
| `global_mask2.nii.gz` | 64 x 64 x 25 | uint8 | `{0, 1}` |
| `global_mask2_v5.nii.gz` | 64 x 64 x 25 x 5 | float32 | `{0, 1}` |
| `mni_downsampled.nii.gz` | 48 x 57 x 48 | float32 | continuous, 0–9533 |

It is a binary mask replicated across four timepoints. Every "time series" the
articles extract from it is a constant vector of length 4. Concretely:

- `Overview.Rmd` builds an ROI, calls `series_roi()`, takes `rowMeans()`, and
  presents the result — four identical numbers — as "the typical `neuroim2`
  workflow".
- `VolumesAndVectors.Rmd` shows `series(vec, 12, 12, 12)` as a voxel time
  course. It is `0 0 0 0`.
- `AnalysisWorkflows.Rmd` computes parcel means over the mask. They are all 1.
- `Smoothing.Rmd` — the entire article — applies `gaussian_blur()`,
  `guided_filter()`, `bilateral_filter()`, `laplace_enhance()`,
  `bilateral_filter_4d()` and `cgb_filter()` to this file. Edge-preserving
  filters are demonstrated on an image with one edge and two intensities.
  `cgb_filter()` builds its graph from **time-series correlations** across four
  identical constant timepoints.

Meanwhile `simulate_fmri()` — exported, documented, purpose-built, with AR(1)
temporal structure, spatial smoothness, heteroscedasticity and latent
components — appears in **zero** vignettes. `mni_downsampled.nii.gz`, the only
image with real intensity structure, appears in two.

This is the single highest-leverage fix in the whole set. Nothing else on this
list changes the reader's experience as much.

### 1.2 The articles narrate their own edit history

Nine sentences tell the reader what the document used to be:

```
ImageVolumes.Rmd:29     This article is now the detailed follow-on to ...
ImageVolumes.Rmd:212    ... use vignette("Resampling"), which now owns that topic directly.
NeuroVector.Rmd:32      This article is now the deeper companion to ...
NeuroVector.Rmd:158     This article now complements the core path rather than replacing it.
regionOfInterest.Rmd:32 This article is now the detailed ROI reference to read after ...
regionOfInterest.Rmd:64 The introductory ROI ... workflow now lives in ...
regionOfInterest.Rmd:220 The basic 4D ROI extraction workflow ... now lives in ...
regionOfInterest.Rmd:252 Searchlight workflows ... are now covered in ...
```

Three of the fourteen titles begin with "Advanced", and all three articles are
defined by what was moved out of them rather than by a topic they own. Readers
should never be able to tell that a refactor happened.

### 1.3 The overview is a table of contents

`Overview.Rmd` contains three separate "where to go next" lists — a numbered
reading path near the top, inline pointers throughout, and a two-part link dump
at the end. Roughly a quarter of the article is navigation. Its actual teaching
content is reproduced almost verbatim in `VolumesAndVectors.Rmd`. A reader is
handed a five-item reading list before they have seen a single brain.

### 1.4 The ordering front-loads a question the reader cannot have yet

The current path is Overview → **ChoosingBackends** → coordinate-systems →
VolumesAndVectors → Resampling. Backend selection is a scaling decision. It is
asked second, before `NeuroVol` and `NeuroVec` have been properly introduced,
and the article's own advice is "start with `read_vec()`, move to mmap or
filebacked when memory forces the issue" — which is an argument for it not
being second.

### 1.5 Ten of fourteen articles contain no figures

| figures | article |
|--:|:--|
| 18 | `slice-visualization.Rmd` |
| 4 | `ImageVolumes.Rmd` |
| 1 | `coordinate-systems.Rmd`, `pipelines.Rmd` |
| **0** | the other ten, including **`Smoothing.Rmd`** and **`Resampling.Rmd`** |

`Smoothing.Rmd` compares six filters and shows the reader the printed `dim()`
of each result. Filtering and resampling are visual operations; describing them
in prose is the wrong medium.

### 1.6 Smaller, but real

- **Test scaffolding in reader-facing chunks.** 25 `stopifnot()` calls across
  six articles. `VolumesAndVectors.Rmd` has 11 — nearly one assertion per code
  chunk, in a document whose job is to be inviting. These are regression tests
  wearing a vignette costume.
- **Duplication.** `spherical_roi() + series_roi()` is taught in five articles.
  Sparse conversion in four. Resampling and downsampling in four. `concat()`,
  `split_clusters()` and `ClusteredNeuroVec` each in three.
- **Style drift.** `#` vs `##` for section headings; `collapse = T` vs `TRUE`;
  4-space-indented chunk bodies in two files; `date: Sys.Date()` on some
  articles (churns output, defeats reproducible builds) and not others; the
  package name inside three titles; 14 near-copies of the setup chunk. In
  `clustered-neurovec.Rmd` the drift produced a bug — `theme_set()` is called
  twice on consecutive lines with two different themes.
- **Rd markup in Markdown.** `Cookbook.Rmd:78` and `:103` contain
  `\code{split_blocks()}`, which renders literally as `\code{split_blocks()}`.
- **Advice that contradicts the package.** `regionOfInterest.Rmd` recommends
  `foreach` + `doParallel` for parallel searchlights. Neither is a dependency;
  the package imports `future` and `future.apply`.
- **Shipped uncertainty.** `clustered-neurovec.Rmd`'s "Integration with existing
  workflows" section is an all-comment chunk containing
  `# (if scale_series is implemented for ClusteredNeuroVec)`.
- **A hole where I/O should be.** `nifti_io.R`, `afni_io.R`, `niml_io.R`,
  `meta_info.R`, `nifti_extensions.R`, `file_format.R` and `header_api.R` have
  essentially no article. `write_vol()` gets four lines. Reading and writing
  files is the first thing every user does.
- **Substantial exports never shown:** `simulate_fmri()`, `NeuroSlice`/`slice()`,
  `automask()`, `clip_level()`, `partition()`, `bounds()`, `enhance_stat_map()`,
  `mapToColors()`, `bootstrap_searchlight()`, `resampled_searchlight()`,
  `BigNeuroVec`, the `*_shape()` generators, `NeuroVecSeq` (named once, never
  constructed), `NeuroHyperVec` (named once, never constructed).

---

## 2. Proposed structure

Nine articles in three tiers. Tier 1 is a sequence and should be read in order;
tiers 2 and 3 are entered on demand.

### Tier 1 — Learn (read in order, ~45 minutes)

**1. `neuroim2.Rmd` — "neuroim2: a tour"**
Named for the package so `vignette("neuroim2")` works and pkgdown lists it
first. One continuous arc, no branching, no reading list until the last
paragraph: read a real anatomical image, **look at it**, read a 4D series,
place an ROI, pull its time course, plot it, write a result to disk. The reader
sees a brain within the first thirty seconds and finishes with one complete
thing they did themselves. Target ~120 lines, three figures, exactly three
forward links at the end.

**2. `spaces-and-coordinates.Rmd`** *(from `coordinate-systems.Rmd`)*
The best article in the set; keep roughly 85% of it. This is the concept that
makes a `NeuroVol` different from an array, and everything downstream — ROIs at
MNI coordinates, overlays, resampling — depends on it. Trim the resampling and
`deoblique()` sections, which belong to article 6, to cross-references.

**3. `volumes-and-vectors.Rmd`** *(merges `VolumesAndVectors` + `ImageVolumes` +
`NeuroVector`)*
The container story, told once: construction from arrays and matrices, masks and
`LogicalNeuroVol`, indexing, `series()`, matrix views, `concat()` /
`sub_vector()` / `split_blocks()`, arithmetic and `Summary` methods. Absorbs the
`Cookbook` recipes for temporal reduction, block splitting and concatenation.

**4. `reading-and-writing.Rmd`** *(new)*
NIfTI read and write, gzip, `read_header()` and what the header actually
carries, sform vs qform priority and what that means for `trans()`, metadata and
extensions, multi-file reads and `read_vol_list()`, `write_vec()`, AFNI/NIML
formats. Currently the largest gap between what the package does and what the
articles cover.

### Tier 2 — Do

**5. `regions-and-searchlights.Rmd`** *(merges `AnalysisWorkflows` +
`regionOfInterest` + `pipelines` + Cookbook's split recipes)*
The full arc in one place: constructors (`spherical_roi`, `spherical_roi_set`,
`cuboid_roi`, `square_roi`, kernels, patches) → extraction (`series_roi`,
`values`, `coords`) → set operations → searchlights (standard, random,
clustered, bootstrap) → split/map/reduce (`split_clusters`, `split_reduce`,
`vols`, `vectors`, `slices`) → `conn_comp` on thresholded maps. Parallelism
example rewritten with `future.apply`, matching the package's own dependencies.

**6. `resampling-and-orientation.Rmd`** *(from `Resampling`, plus the
`reorient`/`deoblique` material from `coordinate-systems` and `ImageVolumes`)*
`resample_to()`, `downsample()`, `reorient()`, `deoblique()`, `as_canonical()`.
With before/after figures — this article currently has none, and grid changes
are exactly what a picture explains best.

**7. `smoothing-and-filtering.Rmd`** *(rewrite of `Smoothing`)*
Same six filters, but on `mni_downsampled` for the 3D spatial filters and on
`simulate_fmri()` output for the 4D ones, with a before/after panel per filter
and a single side-by-side comparison at the end. This is the article that
changes most: same structure, entirely new evidence.

**8. `visualization.Rmd`** *(from `slice-visualization`)*
Already strong. Two changes: move the 90-line hidden data-fabrication chunk into
`vignettes/_common.R` so the article opens on its actual subject, and add the
palette / `mapToColors()` / `enhance_stat_map()` material that has no home today.

### Tier 3 — Scale

**9. `large-data.Rmd`** *(merges `ChoosingBackends` + `clustered-neurovec` +
Cookbook's `as_mmap` recipe)*
"What to do when it does not fit." One escalation ladder — dense → sparse →
mapped → file-backed → `BigNeuroVec` → parcel-level `ClusteredNeuroVec` — with
the decision table, real memory numbers, and the cost of each move. Moving this
from position 2 to position 9 also lets it absorb the parcellation material and
become a substantial article instead of a thin table.

### Disposition of the current fourteen

| current | disposition |
|:--|:--|
| `Overview.Rmd` | rewritten as `neuroim2.Rmd` (1); loses the three link lists |
| `coordinate-systems.Rmd` | → `spaces-and-coordinates.Rmd` (2), trimmed |
| `VolumesAndVectors.Rmd` | → `volumes-and-vectors.Rmd` (3), backbone |
| `ImageVolumes.Rmd` | dissolved into (3) and (6) |
| `NeuroVector.Rmd` | dissolved into (3) and (9) |
| — | **new** `reading-and-writing.Rmd` (4) |
| `AnalysisWorkflows.Rmd` | → `regions-and-searchlights.Rmd` (5), backbone |
| `regionOfInterest.Rmd` | dissolved into (5) |
| `pipelines.Rmd` | dissolved into (5) — 80 lines, 3 examples, not an article |
| `Resampling.Rmd` | → `resampling-and-orientation.Rmd` (6) |
| `Smoothing.Rmd` | → `smoothing-and-filtering.Rmd` (7), rewritten |
| `slice-visualization.Rmd` | → `visualization.Rmd` (8), tightened |
| `ChoosingBackends.Rmd` | → `large-data.Rmd` (9), backbone |
| `clustered-neurovec.Rmd` | dissolved into (9) |
| `Cookbook.Rmd` | dissolved — every recipe has a topical home (see below) |

Cookbook recipes, redistributed: temporal reduction → (3); `split_blocks` → (3);
`concat` for vols and vecs → (3); `as_mmap` → (9); `mapf` + `Kernel` → (5);
`split_clusters` → (5); `split_reduce` → (5); `conn_comp` → (5).

*Judgment call worth your input:* dissolving `Cookbook` is what I would do — a
cookbook tends to become the place things go when nobody decides where they
belong, which is how it accumulated eight unrelated recipes. But keeping a
tenth "Recipes" page is defensible if you want a home for short answers that
outgrow a man page. Say the word and I will keep it.

### Proposed `_pkgdown.yml` articles block

```yaml
articles:
- title: Learn
  desc: Read these four in order. About 45 minutes end to end.
  contents:
  - neuroim2
  - spaces-and-coordinates
  - volumes-and-vectors
  - reading-and-writing
- title: Do
  desc: Task-focused articles. Read the one you need.
  contents:
  - regions-and-searchlights
  - resampling-and-orientation
  - smoothing-and-filtering
  - visualization
- title: Scale
  desc: When the data stops fitting in memory.
  contents:
  - large-data
```

Renaming costs you the old URLs and the old `vignette("Overview")` calls. On the
website that is fully covered by a pkgdown `redirects:` block, which I would add
for all fourteen old paths. In R the old names simply disappear, which is a
NEWS.md line in the next release. If you would rather not pay even that, the
whole plan works with the existing filenames — only the titles and contents
change. Renaming is my recommendation, not a requirement.

---

## 3. Cross-cutting changes

These apply to all nine articles and are, collectively, most of the work.

**3.1 Fix the demo data — three rules, no exceptions.**
- Intensity, filtering, visualization, resampling → `mni_downsampled.nii.gz`.
- Mask and support examples → `global_mask2.nii.gz` (which is what it is).
- Anything with a time axis → `simulate_fmri(mask, n_time = 40, seed = 1)`.

Generate at knit time rather than shipping a new data file: a 64 x 64 x 25 x 40
float series is ~16 MB uncompressed, well past what belongs in `inst/extdata`.
Step 0 of implementation is timing `simulate_fmri()` on the shipped mask; if the
knit cost is material, the fallback is to simulate on a downsampled grid. One
helper in `vignettes/_common.R` means every article gets it in a single line.

**3.2 Move the assertions to the test suite.**
All 25 `stopifnot()` calls come out of the reader-facing chunks and go into
`tests/testthat/test-vignette-examples.R`. This keeps every guarantee those
assertions currently provide — it moves them somewhere they run on every check
rather than somewhere they interrupt a first-time reader. Nothing is deleted
without a replacement.

**3.3 One shared setup.**
`vignettes/_common.R` holds the knitr options, the theme call and the demo-data
helpers. Each article sources it in one line. This removes 14 drifting copies of
the setup chunk and the double-`theme_set()` bug in `clustered-neurovec.Rmd`.

**3.4 A style contract.**
`##` for top-level sections (H1 is the title); sentence-case titles with no
"with neuroim2"; no `date:` field; `collapse = TRUE`; no indented chunk bodies;
no Rd markup in prose; `fig.alt` on every figure.

**3.5 Every article ends the same way.** A short "you can now…" summary, at most
three forward links, and the two or three `?topic` pages that matter. No article
*opens* with a link list.

**3.6 A figure budget.** Minimum one per article. Smoothing gets before/after
per filter; resampling gets a grid comparison; regions gets an ROI drawn on
anatomy; the tour gets three.

**3.7 Length targets.** ~800–1,400 words per Tier-1 article; the four-article
Tier 1 readable in 45 minutes and the whole set in under 90. The current set is
denser in words and thinner in evidence; the rewrite inverts that.

---

## 4. Sequencing

Each phase leaves the package in a shippable state.

| phase | work | why first |
|:--|:--|:--|
| 0 | Time `simulate_fmri()` on the shipped mask; write `vignettes/_common.R`; move the 25 assertions into `tests/` | unblocks everything; the data decision gates the rest |
| 1 | Rewrite `Smoothing` and `Resampling` on real data with figures | largest quality gain per line changed; both are self-contained |
| 2 | Build Tier 1: the tour, spaces, containers, I/O | the reader path; the merges land here |
| 3 | Build Tier 2 merges: regions-and-searchlights, visualization | biggest consolidation (four files into one) |
| 4 | Build Tier 3: large-data | absorbs backends + clustered |
| 5 | `_pkgdown.yml` reorder + redirects, NEWS entry, style sweep, final read-through | cleanup |

Phase 1 alone is worth doing even if nothing else happens: it fixes the two
articles where the current examples demonstrate the least.

---

## 5. Specimen

The tonal difference, on the opening of article 1.

**Now** — `Overview.Rmd`, first 250 words: a paragraph of scope-setting, a
`read_vol()` call printing `dim`/`spacing`/`origin`, then a numbered list of
five other vignettes to read, then a sentence recommending which of the five to
read first.

**Proposed** — `neuroim2.Rmd`:

> ```r
> library(neuroim2)
> anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
> plot(anat)
> ```
>
> *[figure: 3 x 3 axial montage]*
>
> That is a brain, read from a NIfTI file in one line. It is also an ordinary R
> array — `anat[24, 28, 24]` works, `mean(anat)` works, `anat > 0` gives you a
> mask. The difference is that it knows where it is:
>
> ```r
> spacing(anat)   # 4.02 4.02 4.02  — millimetres, not indices
> coord_to_grid(anat, c(-34, -28, 10))
> ```
>
> Ask for MNI coordinate (-34, -28, 10) and you get the voxel that sits there.
> Everything else in this package follows from that one fact, and the rest of
> this article is one pass through it: read an image, look at it, pull a time
> course out of a region, write the answer back to disk.

Same information, arrived at by doing something. The reading list moves to the
last paragraph, where a reader can act on it.

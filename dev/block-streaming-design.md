# Decoded access and storage-aware iteration — contract v2

This supersedes the BS-01 version-1 read_block/execution-engine proposal after
user review. Version 1 was a completed design artifact, not an implemented or
validated runtime architecture. Its completion must not authorize its superseded
work. Mote owns status; see [delivery index](block-streaming-delivery.md).

## Corrected diagnosis

- `linear_access(x, i, ...)` already exists (`R/all_generic.R:130`). Generic
  `series(NeuroVec, integer)` at `R/neurovec.R:727` computes linear indices and
  uses `x[idx]`; ArrayLike4D routes that to backend `linear_access`. A second
  public read primitive is unnecessary for the decoding fix.
- Mapped `linear_access` returns `x@filemap[idx]` without decoding
  (`R/mapped_neurovec.R:212`). Construction at line 168 also fails to retain
  source scale metadata. Fix retention and application together. File-backed
  access at `R/filebacked_neurovec.R:183` uses already-scaled
  `read_mapped_vols`. Source scaling policies live in `R/meta_info_api.R:115`.
- Sparse write #36 and constructor-orientation #31 are separate defects;
  neither is fixed by an additional read abstraction. They are not dependencies
  of the mapped decoding correction unless a focused test proves otherwise.
- `test-backend-parity.R` has five tests, but `test-io-conformance.R` has 27,
  including scaled SHORT fixtures at lines 115 and 452. The gap is coverage of
  those encodings across backends and public accessor paths, not a complete
  absence of scaling fixtures or general I/O assurance.
- `split_reduce(NeuroVec, factor, function)` at `R/common.R:474` reads voxel
  groups through `series`; grouping over time takes a full matrix path.
  `split_blocks` at `R/neurovec.R:1229` already returns lazy deflist blocks over
  time via `sub_vector`. Existing out-of-core computation is possible; its
  costs and uniformity need improvement.

Source was checked at HEAD 4de0c999857947342fcab453829bafb3c7bbbf1a in the
concurrently dirty plotting checkout on 2026-09-25. These are source findings,
not new numerical test results or verified current GitHub issue dispositions.

## Deliver now: existing accessor contract and parity

Document `linear_access` as returning decoded intensity values in the source's
declared units: file slope/intercept applied exactly once, with index order and
duplicates preserved. This does not imply the units are physically known.
Already decoded dense/sparse/FBM values are not scaled again. Sparse positions
outside stored support remain structural zero; stored missing values remain
missing. Reuse existing slope-zero/nonfinite normalization policy rather than
introducing an unrelated metadata-policy migration here.

Preserve public generic signatures, matrix orientations and specialized fast
paths. The generic represents paired linear elements, not a Cartesian block.
Do not translate N pairs into an N-by-N intermediate. Do not route through a
new primitive or make broad accessor refactoring a prerequisite for the fix.

Generate backend × datatype × accessor tests. Core encodings are FLOAT,
SHORT with nontrivial slope/intercept, and UINT8 (package name UBYTE).
Core accessors are `[`, `series`, `linear_access`, `as.matrix`, and
`split_reduce`. Reference both `read_vec(mode="normal")` and independently
specified raw values/scaling so a common bug cannot establish correctness.
Reuse existing fixtures and retain their coverage.

| Representation | Existing contract to preserve / test applicability |
| --- | --- |
| DenseNeuroVec | Series T-by-V; matrix V-by-T; direct array and compiled gather paths. |
| SparseNeuroVec | Stored T-by-K; series T-by-selected-V; matrix full V-by-T with structural zeros. Preserve lazy subclass accessor overrides. |
| BigNeuroVec | FBM T-by-K with sparse hooks. Existing construction can be eager: separate preparation from read/iteration claims. |
| MappedNeuroVec | Volume-major raw payload; fix retained metadata and decoded values on every inherited public read path. |
| FileBackedNeuroVec | Reads full selected volumes before gather; values are already decoded. Cost matters independently of numerical parity. |
| NeuroVecSeq | Mixed components; compare time order, labels, scaling and support with explicit references. Existing square-shape orientation inference is a test target, not a reason for a new engine. |
| ClusteredNeuroVec | Separate NeuroObj/ArrayLike4D class, not NeuroVec. Own series and T-by-C matrix methods; outside support is NA. No explicit linear_access registration found in inspected source. Mark unsupported accessor combinations instead of inventing APIs. |
| NeuroHyperVec | Own series and linear_access in R/neurohypervec.R:319,354; 5D spatial/trial/feature indexing, stored feature-by-trial-by-voxel data. Test with a 5D reference; never flatten it into a 4D backend for parity. |
| NeuroBucket | Declared list-backed NeuroVec with incomplete dedicated access. Audit applicability explicitly; do not advertise an inherited method as supported without a working case. |

Unsupported class/accessor combinations are explicit entries, not silent skips.
Clustered and sparse references must represent the same support semantics;
class-specific orientations are normalized deliberately, never guessed by
dimensions. Cover square and asymmetric fixtures, repeated/out-of-order indices,
selected volumes, sequence boundaries, zeros and missing values. Other encoding
and format coverage remains in the existing I/O suite.

## Deliver next: one storage-aware lazy iterator

Proposed API: `vec_blocks(x, along = "auto", budget = 256 * 1024^2)`.
Budget is bytes of iterator-owned working buffers. Explicit `along` accepts
`"time"` or `"voxel"` for supported 4D inputs. This is storage access, not a
callback scheduler: there is no f argument, map/reduce policy or output writer.

Build on the existing deflist pattern in `split_blocks`, with small internal
backend access hints. Construction inspects metadata but reads no payload.
It yields block records containing decoded values, original voxel/time indices,
actual chosen axis, spatial geometry and available volume labels. A 4D block's
values are a plain T-selected-by-V-selected matrix; existing `as.matrix`
orientations do not change. Consumers reconstruct coordinates from the record,
not from block ordinal or an assumed traversal.

Choose direction before size:

| Storage | Initial auto direction | Reason / condition |
| --- | --- | --- |
| Dense volume-major array | Time | Contiguous volumes; retain specialized indexing. |
| NIfTI mmap | Time | File layout is volume-major; scattered series access is possible but not automatically cheap for a full traversal. |
| File-backed NIfTI | Time | Existing reader loads selected whole volumes; voxel blocks over all times can reread the entire file per block. |
| Sparse T-by-K matrix / FBM | Voxel | Complete stored columns are contiguous; expand only selected full-grid positions and retain structural zeros. |
| Sequence | Component-local direction | Do not force different stores into one expensive global direction; each record reports its axis and original global times. |

These defaults must be supported by BS-03 measurements and read counters.
Expose chosen traversal and cost limitations. Explicit along requests are never
silently reinterpreted; disclose an expensive traversal or reject it when the
budget cannot be respected. Auto must not perform O(number-of-voxel-blocks ×
file-size) I/O for a full traversal of ordinary file-backed data.

V1 iterator targets the six primary 4D backends. Clustered and HyperVec remain
in the accessor applicability audit and tests; the iterator gives an explicit
unsupported-class error until a shape-preserving adapter is separately scoped.
This is a stated feature boundary, not omission of their existing APIs.

Use existing validated accessors/volume reads; no wholesale migration of
`series`, `linear_access` or `as.matrix`. Account for returned values, index
construction and reader scratch when choosing block length. If a minimum
legal block (e.g. one full volume) exceeds budget, fail before reading instead
of claiming that a small returned slice was cheap. New sub-volume I/O is not
mandatory for this first iterator; any supported limit must be documented.

Iteration is lazy and repeatable on an unchanged source; no intensity data is
modified. Verify deflist realization/caching behavior: the iterator itself must
not retain every previously yielded block. Caller-retained blocks, original
source storage, mapped page residency and R/OS overhead are outside its buffer
budget and must be disclosed when reporting RSS. Opening gzip caches and
constructing BigNeuroVec are separate memory/disk costs. No hidden full-file
copy or persistent temporary cache is justified by a budget argument.

For mixed sequences, validate spatial compatibility, preserve component support
semantics and labels, and emit global time coordinates without eager concatenation.
The iterator cannot infer whether a caller's algorithm needs complete time
courses or how partial reductions combine. Callers choose/validate traversal
and own accumulated state, result size, algorithms and concurrency.

## Explicit cuts and ownership

No public read_block, vec_map, vec_reduce, vec_apply_searchlight, output sink,
out="filebacked" path, writer migration, or internal parallel executor belongs
to this delivery. Keep searchlight_indices geometry-only as its NEWS contract
intends; rMVPA/fmri callers own their searchlight drivers and analysis policy.

RcppParallel cannot invoke arbitrary user R callbacks on worker threads. Future
multisession use would need qualified source reopening/serialization rather
than transferring live mmap/FBM handles. This plan implements neither mechanism
and makes no parallel-callback promise. Existing independent parallel code is
not removed or changed by this re-scope.

## Acceptance mapping

| Gate | Evidence | Tickets |
| --- | --- | --- |
| Decoding fix | Scaled SHORT regression fails old mmap path and passes fixed metadata + linear_access; inherited accessors agree | BS-08 |
| Class/API applicability | 4D/5D and clustered support/axis conventions documented and tested without invented APIs | BS-11 |
| Public parity | Generated encoding/backend/accessor matrix plus retained 27-test conformance suite | BS-02 |
| Storage cost | Actual bytes/volume reads for both directions; preparation separated from iteration | BS-03 |
| Laziness / budget | Metadata-only creation; bounded per-block buffers; no retained all-block cache; minimum-block refusal | BS-13/09 |
| Sequence composition | Component-local cost, exact global coordinates and no dense concatenation | BS-10 |
| Qualification | Value reconstruction, fixed-budget growth, linear auto I/O and existing fast-path baselines | BS-25 |
| Documentation / package | Runnable caller-owned loop, explicit unsupported cases, full tests and installed-artifact checks | BS-26/27 |

The near-term fix is independently deliverable; iterator work follows its own
measurement/qualification gates. No estimate of two or three sessions covers
all classes and API surfaces. This specification revision has not run runtime
tests or fixed the mapped bug; those remain open implementation tickets.

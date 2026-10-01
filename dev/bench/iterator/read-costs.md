# Storage traversal evidence (BS-03)

Run `Rscript dev/bench/iterator/read-costs.R` from the package root.
Measured 2026-09-25 on the working tree based on HEAD
`9a94840c1ad624473133d0acc4c8927df057f9b7` (with BS-08 mapped decoding).
Synthetic scaled SHORT: 16 × 16 × 8 × 24, eight blocks in each direction.
Raw log: `/tmp/bs03-costs.log`. These are small diagnostic workloads, not throughput
claims; first-use dispatch/compiler allocations affect elapsed/allocation totals.

| Storage | Time blocks allocated bytes | Voxel blocks allocated bytes | Time / voxel elapsed seconds |
|---|---:|---:|---:|
| Dense | 418384 | 1367648 | 0.000 / 0.015 |
| Mapped | 4365848 | 4574904 | 0.016 / 0.004 |
| File-backed | 7647320 | 18899128 | 0.222 / 0.218 |
| Sparse | 6021624 | 419952 | 0.019 / 0.000 |
| FBM | 5991760 | 856496 | 0.007 / 0.001 |

At the actual compiled volume-reader boundary, file-backed time traversal requested
24 volumes (98304 stored bytes); voxel traversal requested 192 (786432 bytes).
These count requested payload, not OS physical I/O: cache residency is independent.
Mapped pages are not counted at that boundary. Their volume-major file layout
supports time traversal, but this small warm-page timing does not establish a
mapped speed advantage. Dense chooses time by measured allocations and layout;
sparse/FBM choose voxel by measured allocations and T-by-K column layout.

Preparation seconds were dense .076, mapped .089, file-backed .001, sparse .002,
FBM .047. FBM construction is eager (prep_sparsenvec followed by as_FBM), and its
initial source allocations are outside the iterator. The fixture is uncompressed;
gzip source/cache preparation is likewise a separate pre-iteration activity and
is not qualified here. Mixed sequences must preserve each component's direction
rather than choose a single global direction (adapter remains BS-10).

## Budget implications

One time block returns V-by-selected-T values (then transposes to T-by-V); one
voxel block returns all-T-by-selected-V values. File-backed gather additionally
reads full requested volumes; a one-voxel/all-time block therefore requires the
whole decoded file in reader scratch. Its future adapter must budget that scratch
or reject such a request. At least a full spatial volume must fit for its auto
traversal. Sparse/FBM should use series over selected columns; time traversal uses
paired linear access and must budget its mapping/index buffers.

Dense primitive `.subset(x, indices)` was separately checked against independent
values: largest allocation 16432 bytes for 2048 doubles, versus 393216 bytes of
source payload. BS-13 can use this bounded primitive on the exact dense class
without extracting/coercing the entire data part. Arbitrary subclasses must not
inherit storage-cost claims that their overrides could invalidate.

Use conservative per-cell buffer accounting (including index/gather/transpose
scratch), not just returned matrix bytes. R object headers, allocator/GC overhead,
source storage, caller-retained blocks and mapped page residency must remain
explicitly distinct from this buffer budget. Timing alone cannot establish a
memory bound. BS-13 tests additionally profile individual realization with growing
sources and inspect deflist's non-memoising behavior.

# Accessor applicability (BS-11)

`V` is full-grid voxels, `K` stored voxels/clusters, `T` time/trials,
`F` features. This table describes supported calls, not shape-based inference.

| Class | `[` | `linear_access` | `series` | `as.matrix` | `split_reduce` |
|---|---|---|---|---|---|
| DenseNeuroVec | linear / 4D | paired linear values | T × selected V | V × T | voxel groups |
| SparseNeuroVec | linear / 4D | paired linear values | T × selected V | V × T | voxel groups |
| BigNeuroVec | linear / 4D | paired linear values | T × selected V | V × T | voxel groups |
| MappedNeuroVec | linear / 4D | decoded paired values | T × selected V | V × T | voxel groups |
| FileBackedNeuroVec | linear / 4D | decoded paired values | T × selected V | V × T | voxel groups |
| NeuroVecSeq | linear / 4D | component values in requested order | T × selected V | V × T | voxel groups |
| ClusteredNeuroVec | paired xyz plus times; full single volume | unsupported | xyz coordinates, T × selected V | T × K clusters | unsupported |
| NeuroHyperVec | Cartesian 5D, not single-index linear | paired 5D linear values | single xyz, F × T | unsupported | unsupported |
| NeuroBucket | unsupported | unsupported | unsupported | unsupported | unsupported |

Singleton series calls may drop dimensions. Clustered outside-mask positions are
NA; sparse/Big/HyperVec outside-mask positions are structural zero. Stored missing
values remain missing. Clustered matrix values cannot be compared to voxel
matrices without explicit broadcasting. HyperVec tests construct their reference
by independent spatial/trial/feature coordinates. Bucket has class inheritance
but no working read implementation; regression tests record this boundary.

All six ordinary 4D backends have T-by-V series contracts. The parity matrix
exercises square component results; NeuroVecSeq now keeps their declared
orientation instead of inferring it from dimensions and transposing them.

The parity matrix uses **voxel-group** split_reduce. Time-group split_reduce is
excluded explicitly: its implementation passes a V-by-T matrix to a method
expecting its factor along rows, so it is not a validated time-group API. This
pre-existing reduction defect does not change accessor parity or authorize an
execution engine. Constructor matrix orientation is avoided by constructing
sparse/Big objects from unambiguous 4D arrays.

Evidence: `tests/testthat/test-accessor-applicability.R` plus existing hypervec,
clustered, backend-parity and I/O tests. The iterator in BS-13 supports concrete
dense/sparse/Big storage; mapped/file-backed and sequence adapters remain BS-09
and BS-10. It rejects unsupported classes rather than coercing to dense.

# ClusteredNeuroVec Class

A class representing a 4D neuroimaging dataset where voxels are grouped
into clusters. Each cluster has a single time-series that is shared by
all voxels within that cluster.

## Slots

- `cvol`:

  A
  [`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md)
  object defining cluster assignments

- `ts`:

  A numeric matrix of dimensions T x K (time points x clusters)

- `cl_map`:

  An integer vector mapping each voxel to its cluster ID (0 for outside
  mask)

- `label`:

  A character string label for the object

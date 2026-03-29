# ROICoords

A class representing a region of interest (ROI) in a brain image,
defined by a set of coordinates. This class stores the geometric space
of the image and the coordinates of the voxels within the ROI.

## Slots

- `space`:

  An instance of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  representing the geometric space of the image data.

- `coords`:

  A `matrix` containing the coordinates of the voxels within the ROI.
  Each row represents a coordinate as, e.g. (i, j, k).

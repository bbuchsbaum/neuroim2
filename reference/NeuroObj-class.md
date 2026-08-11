# NeuroObj Class

Base class for all neuroimaging data objects with a Cartesian spatial
representation. This class provides a foundation for more specific
neuroimaging data structures.

## Slots

- `space`:

  An object of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  representing the geometry of the image object.

- `header`:

  A `list` holding the raw header of the file the object was read from,
  empty for objects built in memory. `NeuroSpace` records the geometry;
  this records everything else the file said – repetition time, units,
  intent, description, slice timing, the qform/sform codes – so that
  writing the object back does not have to invent it. Use
  [`header`](https://bbuchsbaum.github.io/neuroim2/reference/header-methods.md)
  to read it rather than touching the slot.

## See also

[`NeuroSpace-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md),
[`NeuroSlice-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSlice-class.md),
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)

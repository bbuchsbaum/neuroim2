# FileSource Class

Base class for representing a data source for images. The purpose of
this class is to provide a layer in between low level IO and image
loading functionality.

## Slots

- `meta_info`:

  An object of class
  [`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)
  containing meta information for the data source.

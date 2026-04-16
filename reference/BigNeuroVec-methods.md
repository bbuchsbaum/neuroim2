# Create a Memory-Mapped Neuroimaging Vector

Creates a BigNeuroVec object, which represents a large neuroimaging
vector using memory-mapped file storage. This allows working with
neuroimaging data that is too large to fit in memory.

## Usage

``` r
BigNeuroVec(
  data,
  space,
  mask,
  label = "",
  volume_labels = character(),
  type = c("double", "float", "integer"),
  backingfile = tempfile()
)
```

## Arguments

- data:

  The input data to be stored

- space:

  A NeuroSpace object defining the spatial properties

- mask:

  A logical mask indicating which voxels contain data

- label:

  Optional character string label for the vector

- volume_labels:

  Optional character vector of per-volume labels

- type:

  Storage type, one of "double", "float", or "integer"

- backingfile:

  Path to the file used for memory mapping (defaults to tempfile())

## Value

A new BigNeuroVec object that provides memory-efficient access to large
neuroimaging data through memory mapping. The object contains the
spatial properties, mask, and memory-mapped data storage.

## Examples

``` r
# \donttest{
# Load an example 4D brain image
example_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
example_4d_image <- read_vec(example_file)

# Create a mask (e.g., selecting voxels with values > 0)
mask <- array(as.vector(example_4d_image[,,,1]) > 0,
             dim = dim(example_4d_image)[1:3])

if(requireNamespace("bigstatsr", quietly = TRUE)) {
  # Create a BigNeuroVec with memory mapping
  big_vec <- BigNeuroVec(data = example_4d_image@.Data,
                         space = space(example_4d_image),
                         mask = mask,
                         label = "Example BigNeuroVec")
  print(big_vec)
}
#> <BigNeuroVec> [936.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Time Points   : 4
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
# }
```

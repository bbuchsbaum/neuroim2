# Compile spherical searchlights to full-volume linear indices

Returns the geometry of every spherical searchlight centred on a nonzero
voxel of `mask`, without constructing ROI objects, extracting analysis
data, or changing parallel execution state. Neighborhoods are compiled
sequentially by the same cached-offset and compiled clipping machinery
used by
[`searchlight_coords`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight_coords.md).

## Usage

``` r
searchlight_indices(mask, radius, nonzero = TRUE)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object defining the searchlight centres and, when `nonzero = TRUE`,
  the allowed neighborhood members.

- radius:

  A positive numeric scalar giving the spherical radius in millimetres.

- nonzero:

  A single logical value. If `TRUE` (the default), each neighborhood is
  restricted to finite, nonzero voxels of `mask`. It never changes the
  centres: every nonzero mask voxel is a centre.

## Value

A list-like object of class `searchlight_indices` with one integer
vector per centre. Every value is a stable, 1-based full-volume linear
index using R's column-major array order. Centre order is
`which(mask != 0)`. The result carries the following documented
attributes: `center_indices`, `space`, `radius`, and `nonzero`. The
object is eagerly compiled but contains indices only.

## Details

The full-volume index contract means an input whose first three
dimensions contain more than `.Machine$integer.max` voxels is rejected.
The function has no `cores` argument and never inspects or modifies
[`future::plan()`](https://future.futureverse.org/reference/plan.html).

## Examples

``` r
mask_data <- array(FALSE, c(7, 7, 7))
mask_data[2:6, 2:6, 2:6] <- TRUE
mask <- LogicalNeuroVol(mask_data, NeuroSpace(c(7, 7, 7)))

neighborhoods <- searchlight_indices(mask, radius = 2)
length(neighborhoods)
#> [1] 125
attr(neighborhoods, "center_indices")
#>   [1]  58  59  60  61  62  65  66  67  68  69  72  73  74  75  76  79  80  81
#>  [19]  82  83  86  87  88  89  90 107 108 109 110 111 114 115 116 117 118 121
#>  [37] 122 123 124 125 128 129 130 131 132 135 136 137 138 139 156 157 158 159
#>  [55] 160 163 164 165 166 167 170 171 172 173 174 177 178 179 180 181 184 185
#>  [73] 186 187 188 205 206 207 208 209 212 213 214 215 216 219 220 221 222 223
#>  [91] 226 227 228 229 230 233 234 235 236 237 254 255 256 257 258 261 262 263
#> [109] 264 265 268 269 270 271 272 275 276 277 278 279 282 283 284 285 286
index_to_grid(mask, neighborhoods[[1]])
#>       [,1] [,2] [,3]
#>  [1,]    2    2    2
#>  [2,]    2    2    3
#>  [3,]    2    2    4
#>  [4,]    2    3    2
#>  [5,]    2    3    3
#>  [6,]    2    4    2
#>  [7,]    3    2    2
#>  [8,]    3    2    3
#>  [9,]    3    3    2
#> [10,]    3    3    3
#> [11,]    4    2    2
```

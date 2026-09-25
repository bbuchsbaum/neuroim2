# Overlay fixed and moving edge maps on a background volume

Displays a structural/reference background with two edge channels
rendered in distinct colors. This is intended for registration QC where
fixed/template edges and moving/result edges should coincide.

## Usage

``` r
plot_edge_overlay(
  bgvol,
  fixed_edges,
  moving_edges,
  zlevels = NULL,
  along = 3L,
  bg_cmap = "grays",
  fixed_color = "#2bb8f0",
  moving_color = "#ff3b30",
  agree_color = "#ffd23f",
  bg_range = c("robust", "data"),
  edge_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  edge_thresh = 0,
  edge_alpha = 0.9,
  ncol = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = FALSE,
  style = c("light", "dark", "report"),
  labels = c("fixed", "moving"),
  legend = TRUE,
  bg_dim = 0.75,
  crop = TRUE,
  unit = c("index", "mm"),
  annotate = TRUE,
  assemble = TRUE,
  n_slices = 6L,
  thin = TRUE,
  compute_edges = FALSE,
  canvas = NULL,
  interpolate = TRUE,
  focus_brain = TRUE
)
```

## Arguments

- bgvol:

  Background 3D volume.

- fixed_edges:

  Edge map for the fixed/reference image on the same NeuroSpace grid as
  \`bgvol\`.

- moving_edges:

  Edge map for the moving/result image on the same NeuroSpace grid as
  \`bgvol\`.

- zlevels:

  Slices to plot: indices along \`along\` (`unit = "index"`) or world
  coordinates (`unit = "mm"`). `NULL` (default) picks `n_slices` slices
  spread over the brain.

- along:

  Native voxel-grid axis for slicing. Display orientation is inferred
  from the image affine.

- bg_cmap:

  Background palette.

- fixed_color, moving_color:

  Overlay colors for the two edge maps.

- agree_color:

  Colour for contour pixels present in both edge maps (where the images
  agree), keyed "both"; `NA` turns the agreement layer and its key entry
  off.

- bg_range, edge_range:

  "robust" or "data" intensity scaling.

- probs:

  Quantiles for robust scaling.

- edge_thresh:

  Values below this edge magnitude are transparent.

- edge_alpha:

  Global alpha for edge overlays.

- ncol:

  Number of columns (`NULL` = chosen to fill the canvas).

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if `TRUE`, also print the figure immediately (and return it
  invisibly). By default the figure is returned visibly.

- style:

  Visual style: `"light"`, `"dark"` or `"report"`.

- labels:

  Length-2 character vector naming the fixed and moving edge maps in the
  key (e.g. `c("MNI template", "registered EPI")`).

- legend:

  Logical; draw the colour key under the panels.

- bg_dim:

  Brightness multiplier (0–1) for the background so the edge colours
  read clearly.

- crop:

  Logical; crop panels to the head bounding box.

- unit:

  `"index"` or `"mm"`: how `zlevels` is interpreted. Panels are labelled
  in world coordinates.

- annotate:

  Logical; draw orientation letters on the first panel.

- assemble:

  Logical; return one patchwork figure (default) or the list of panel
  ggplots.

- n_slices:

  Number of automatically chosen slices.

- thin:

  Logical; thin edge bands to their ridge lines (1–2 pixels) so the two
  contours can be compared precisely. Default `TRUE`.

- compute_edges:

  Logical; if `TRUE`, `fixed_edges` and `moving_edges` are ordinary
  intensity images (e.g. the template and the registered image) and
  their edges are computed per slice as the in-plane gradient magnitude,
  keeping the strongest 12% of head voxels.

- canvas:

  Optional `c(width, height)` in inches to fit the layout to (default:
  the open device).

- interpolate:

  Logical; smooth the background and anti-alias the edge contours
  (default `TRUE`).

- focus_brain:

  Logical; fade edges that lie outside the (dilated) brain mask of
  `bgvol` (bright tissue connected to the centre of the head, which
  excludes scalp fat), so scalp and skull contours do not dominate the
  comparison. Default `TRUE`.

## Value

A figure (class `neuro_fig`, a patchwork whose layout is re-fitted to
the device it is drawn on) when `assemble = TRUE` or a named list of
panel ggplots (`assemble = FALSE`); invisibly when `draw = TRUE`.

## See also

Other plot_neuro:
[`plot_checkerboard()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_checkerboard.md),
[`plot_montage()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_montage.md),
[`plot_ortho()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_ortho.md),
[`plot_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)

## Examples

``` r
# \donttest{
fixed <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
moving <- fixed
moving[] <- fixed[c(2:dim(fixed)[1], 1), , ]   # a one-voxel shift
p <- plot_edge_overlay(fixed, fixed, moving, compute_edges = TRUE,
                       labels = c("template", "registered"))
ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
# }
```

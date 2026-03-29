# Simulate fMRI Data

Generates synthetic 4D fMRI data with realistic spatiotemporal
properties including temporal autocorrelation, spatial smoothness,
heteroscedasticity, and optional global signal fluctuations and latent
components.

## Usage

``` r
simulate_fmri(
  mask,
  n_time,
  TR = 2,
  spatial_fwhm = 6,
  ar_mean = 0.45,
  ar_sd = 0.08,
  noise_sd = 1,
  hetero_fwhm = 20,
  hetero_strength = 0.6,
  global_amp = 0.2,
  global_rho = 0.85,
  n_factors = 4,
  factor_fwhm = 12,
  factor_rho = 0.8,
  seed = NULL,
  return_centered = TRUE
)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object defining the brain mask region. Can be binary or continuous
  (non-zero values define the mask).

- n_time:

  Integer specifying the number of time points to simulate.

- TR:

  Numeric value for the repetition time in seconds (default = 2.0).
  Currently used only for metadata.

- spatial_fwhm:

  Numeric value specifying the spatial smoothness in mm (full width at
  half maximum) applied to each timepoint (default = 6).

- ar_mean:

  Numeric value for the mean of the AR(1) coefficient distribution
  across voxels (default = 0.45).

- ar_sd:

  Numeric value for the standard deviation of the AR(1) coefficient
  distribution (default = 0.08).

- noise_sd:

  Numeric value for the baseline noise standard deviation (default =
  1.0).

- hetero_fwhm:

  Numeric value for the spatial scale (FWHM in mm) of the
  heteroscedasticity field (default = 20).

- hetero_strength:

  Numeric value controlling the strength of spatial heteroscedasticity
  on log scale (default = 0.6).

- global_amp:

  Numeric value for the amplitude of global signal fluctuations as a
  fraction of median noise (default = 0.2). Set to 0 to disable.

- global_rho:

  Numeric value for the AR(1) coefficient of global signal (default =
  0.85).

- n_factors:

  Integer specifying the number of latent spatial components (default =
  4). Set to 0 to disable.

- factor_fwhm:

  Numeric value for the spatial smoothness (FWHM in mm) of latent
  component maps (default = 12).

- factor_rho:

  Numeric value for the AR(1) coefficient of latent component time
  courses (default = 0.8).

- seed:

  Integer seed for random number generation (default = NULL for no
  seed).

- return_centered:

  Logical indicating whether to center each voxel's time series to mean
  zero (default = TRUE).

## Value

A
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
object containing the simulated 4D fMRI data.

## Details

The simulation combines several realistic features:

- Voxel-wise AR(1) temporal autocorrelation with spatial variation

- Spatial smoothing applied to innovations for realistic spatial
  correlation

- Heteroscedastic noise with smooth spatial modulation

- Optional low-frequency global signal fluctuations

- Optional latent spatial components resembling resting-state networks

The spatial smoothing uses the package's optimized
[`gaussian_blur`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md)
function for efficiency.

## References

Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise
ratio and contrast-to-noise ratio for fMRI data. PloS one, 8(11),
e77089.

## Examples

``` r
# Create a simple spherical mask
dims <- c(32, 32, 20)
mask_array <- array(FALSE, dims)
center <- dims / 2
for (i in 1:dims[1]) {
  for (j in 1:dims[2]) {
    for (k in 1:dims[3]) {
      if (sum(((c(i,j,k) - center) / (dims/3))^2) <= 1) {
        mask_array[i,j,k] <- TRUE
      }
    }
  }
}

mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))

# Simulate 100 time points
sim_data <- simulate_fmri(mask, n_time = 100, seed = 42)

# Check dimensions
dim(sim_data)  # Should be c(32, 32, 20, 100)
#> [1]  32  32  20 100
```

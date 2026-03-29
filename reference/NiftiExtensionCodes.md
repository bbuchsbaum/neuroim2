# Known NIfTI Extension Codes

A named integer vector of registered NIfTI extension codes. These codes
identify the format/type of extension data.

## Usage

``` r
NiftiExtensionCodes
```

## Format

Named integer vector where names describe the extension type:

- ignore:

  0 - Unknown/private format (not recommended)

- DICOM:

  2 - DICOM format (attribute tags and values)

- AFNI:

  4 - AFNI group (ASCII XML attributes)

- comment:

  6 - Plain text comment

- XCEDE:

  8 - XCEDE format

- jimdiminfo:

  10 - JIM dimension info

- workflow_fwds:

  12 - Workflow forwards

- FreeSurfer:

  14 - FreeSurfer format

- pypickle:

  16 - Python pickle

- MiND_ident:

  18 - MiND identifier

- b_value:

  20 - B-value (diffusion)

- spherical_direction:

  22 - Spherical direction

- DT_component:

  24 - DT component

- SHC_degreeorder:

  26 - SHC degree order

- voxbo:

  28 - VoxBo format

- Caret:

  30 - Caret format

- CIFTI:

  32 - CIFTI format

- variable_frame_timing:

  34 - Variable frame timing

- eval:

  38 - Eval

- MATLAB:

  40 - MATLAB format

- Quantiphyse:

  42 - Quantiphyse

- MRS:

  44 - MRS NIfTI

## Examples

``` r
# Get the code for AFNI extensions
NiftiExtensionCodes["AFNI"]  # Returns 4
#> AFNI 
#>    4 

# Get the name for a code
names(NiftiExtensionCodes)[NiftiExtensionCodes == 4]  # Returns "AFNI"
#> [1] "AFNI"
```

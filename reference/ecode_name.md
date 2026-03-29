# Get Extension Code Name

Returns the name associated with a NIfTI extension code.

## Usage

``` r
ecode_name(ecode)
```

## Arguments

- ecode:

  Integer extension code.

## Value

Character string with the extension name, or "unknown" if not found.

## Examples

``` r
ecode_name(4L)   # Returns "AFNI"
#> [1] "AFNI"
ecode_name(6L)   # Returns "comment"
#> [1] "comment"
ecode_name(999L) # Returns "unknown"
#> [1] "unknown"
```

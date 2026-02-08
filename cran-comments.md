## Resubmission

This is a resubmission. In this version I have:

* Fixed Windows build failure: added `src/Makevars.win` to correctly link TBB libraries without the non-existent `-lRcppParallel` library flag.
* Fixed PDF manual generation: replaced Unicode characters (Greek letters Δ, σ and symbols ≈, ×, –) in roxygen documentation with ASCII equivalents to resolve LaTeX errors.

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

Local check:

- `R CMD check --as-cran neuroim2_0.8.5.tar.gz` (macOS Sonoma 14.3, R 4.5.1).

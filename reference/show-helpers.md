# Internal helpers for consistent show() formatting

All show() methods in the package should use these helpers to ensure a
uniform visual style. The output uses cli for section rules and plain
[`cat()`](https://rdrr.io/r/base/cat.html) for field lines.

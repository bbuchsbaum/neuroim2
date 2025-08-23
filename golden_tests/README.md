# Golden Tests for neuroim2

This directory contains golden tests for the neuroim2 package, designed to ensure semantic equivalence across different language implementations (R, Python, Rust).

## Overview

Golden tests focus on validating numeric outputs rather than implementation details. Each test specifies:
- Input data and parameters
- Expected numeric outputs with tolerances
- Semantic description of the algorithm
- Language-specific implementations

## Test Structure

```
golden_tests/
├── specs/                      # Test specifications in XML
│   └── core/
│       ├── spatial_reference/  # Spatial coordinate system tests
│       ├── volume_operations/  # 3D volume operations
│       ├── vector_operations/  # 4D time series operations
│       ├── spatial_algorithms/ # Algorithms like connected components
│       └── io_operations/      # File I/O tests
├── schema/                     # XML schema definition
└── validators/                 # Language-specific validators
    └── R/                      # R test validator
```

## The 10 Essential Tests

1. **NeuroSpace Construction and Coordinate Transformation** (`neurospace_construction.xml`)
   - Tests 3D space creation and voxel/world coordinate transformations

2. **NeuroVol Basic Construction** (`neurovol_construction.xml`)
   - Tests 3D volume creation from numeric arrays

3. **NeuroVol Arithmetic Operations** (`neurovol_arithmetic.xml`)
   - Tests element-wise operations between volumes and scalars

4. **NeuroVol Indexing and Subsetting** (`neurovol_indexing.xml`)
   - Tests various indexing patterns and dimension handling

5. **Sparse NeuroVol Operations** (`sparse_neurovol.xml`)
   - Tests sparse volume creation and operations

6. **NeuroVec Construction and Series Extraction** (`neurovec_construction.xml`)
   - Tests 4D data construction and time series extraction

7. **Connected Components** (`connected_components.xml`)
   - Tests clustering algorithm on binary volumes

8. **ROI Operations** (`roi_operations.xml`)
   - Tests region of interest creation and extraction

9. **Basic IO Read/Write Cycle** (`basic_io_cycle.xml`)
   - Tests data persistence and metadata preservation

10. **Concatenation Operations** (`concatenation.xml`)
    - Tests combining multiple 3D volumes into 4D

## Running Tests in R

```r
# Source the validator
source("golden_tests/validators/R/validate_golden_tests.R")

# Validate a single test
validate_golden_test("golden_tests/specs/core/spatial_reference/neurospace_construction.xml")

# Validate all tests
validate_all_golden_tests("golden_tests/specs")

# Use with testthat
test_golden()
```

## Adding New Language Implementations

1. Create a validator in your language under `validators/[language]/`
2. Parse XML test specifications
3. Implement numeric validation logic
4. Add your implementation to the `<implementations>` section of each test
5. Update `<propagation_status>` when complete

## Test Design Principles

- Use small, deterministic test data (e.g., 4x4x4 volumes)
- Focus on numeric equivalence with appropriate tolerances
- Include semantic descriptions for clarity
- Test both normal cases and edge cases
- Ensure tests are language-agnostic

## Tolerance Guidelines

- Exact values (integers): tolerance = 0
- Floating-point: tolerance = 1e-10 for deterministic operations
- Statistical operations: tolerance = 1e-6 or relative tolerance
- Iterative algorithms: looser tolerances as appropriate

## Status

All 10 essential tests have been implemented for R. Python and Rust implementations are pending.
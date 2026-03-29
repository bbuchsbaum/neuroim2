# File Format Operations for Neuroimaging Data

A collection of methods for handling neuroimaging file formats with
separate header and data files (e.g., ANALYZE, NIFTI). These methods
provide functionality for file name validation, extension handling, and
file path manipulation.

## File Format Structure

Neuroimaging formats often use paired files:

- A header file (e.g., '.hdr') containing metadata

- A data file (e.g., '.img') containing the actual image data

## Common Operations

- Validating file names against format specifications

- Converting between header and data file names

- Checking file existence and compatibility

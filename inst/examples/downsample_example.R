# Example: Fast downsampling of 4D functional neuroimaging data
library(neuroim2)

# Create a sample 4D functional image (64x64x32 voxels, 20 time points)
# In practice, this would be loaded from a file
set.seed(123)
data <- array(rnorm(64*64*32*20), dim = c(64, 64, 32, 20))
space <- NeuroSpace(dim = c(64, 64, 32, 20), 
                    origin = c(0, 0, 0),
                    spacing = c(3, 3, 3))  # 3mm isotropic voxels

# Create DenseNeuroVec object
nvec <- DenseNeuroVec(data, space)
cat("Original image:\n")
cat("  Dimensions:", dim(nvec), "\n")
cat("  Spacing:", spacing(nvec)[1:3], "mm\n")
cat("  Memory size:", format(object.size(nvec), units = "MB"), "\n\n")

# Method 1: Downsample by factor (0.5 = half resolution in each dimension)
nvec_factor <- downsample(nvec, factor = 0.5)
cat("Downsampled by factor 0.5:\n")
cat("  Dimensions:", dim(nvec_factor), "\n")
cat("  Spacing:", spacing(nvec_factor)[1:3], "mm\n")
cat("  Memory size:", format(object.size(nvec_factor), units = "MB"), "\n\n")

# Method 2: Downsample to target spacing (6mm voxels)
nvec_spacing <- downsample(nvec, spacing = c(6, 6, 6))
cat("Downsampled to 6mm spacing:\n")
cat("  Dimensions:", dim(nvec_spacing), "\n")
cat("  Spacing:", spacing(nvec_spacing)[1:3], "mm\n")
cat("  Memory size:", format(object.size(nvec_spacing), units = "MB"), "\n\n")

# Method 3: Downsample to specific output dimensions
nvec_dims <- downsample(nvec, outdim = c(16, 16, 8))
cat("Downsampled to 16x16x8 voxels:\n")
cat("  Dimensions:", dim(nvec_dims), "\n")
cat("  Spacing:", spacing(nvec_dims)[1:3], "mm\n")
cat("  Memory size:", format(object.size(nvec_dims), units = "MB"), "\n\n")

# Non-uniform downsampling (different factors per dimension)
nvec_nonuniform <- downsample(nvec, factor = c(0.5, 0.5, 0.25))
cat("Non-uniform downsampling (0.5, 0.5, 0.25):\n")
cat("  Dimensions:", dim(nvec_nonuniform), "\n")
cat("  Spacing:", spacing(nvec_nonuniform)[1:3], "mm\n")
cat("  Memory size:", format(object.size(nvec_nonuniform), units = "MB"), "\n")

# Note: The time dimension (4th dimension) is always preserved
cat("\nNote: Time dimension preserved in all cases (", dim(nvec)[4], "time points )\n")
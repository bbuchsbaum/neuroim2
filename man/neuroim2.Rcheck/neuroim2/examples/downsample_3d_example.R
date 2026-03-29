# Example: Fast downsampling of 3D brain volumes
library(neuroim2)

# Create a sample 3D brain volume (64x64x32 voxels)
# In practice, this would be loaded from a file
set.seed(123)
data <- array(rnorm(64*64*32), dim = c(64, 64, 32))
space <- NeuroSpace(dim = c(64, 64, 32), 
                    origin = c(0, 0, 0),
                    spacing = c(3, 3, 3))  # 3mm isotropic voxels

# Create DenseNeuroVol object
vol <- DenseNeuroVol(data, space)
cat("Original volume:\n")
cat("  Dimensions:", dim(vol), "\n")
cat("  Spacing:", spacing(vol), "mm\n")
cat("  Memory size:", format(object.size(vol), units = "MB"), "\n\n")

# Method 1: Downsample by factor (0.5 = half resolution in each dimension)
vol_factor <- downsample(vol, factor = 0.5)
cat("Downsampled by factor 0.5:\n")
cat("  Dimensions:", dim(vol_factor), "\n")
cat("  Spacing:", spacing(vol_factor), "mm\n")
cat("  Memory size:", format(object.size(vol_factor), units = "MB"), "\n\n")

# Method 2: Downsample to target spacing (6mm voxels)
vol_spacing <- downsample(vol, spacing = c(6, 6, 6))
cat("Downsampled to 6mm spacing:\n")
cat("  Dimensions:", dim(vol_spacing), "\n")
cat("  Spacing:", spacing(vol_spacing), "mm\n")
cat("  Memory size:", format(object.size(vol_spacing), units = "MB"), "\n\n")

# Method 3: Downsample to specific output dimensions
vol_dims <- downsample(vol, outdim = c(16, 16, 8))
cat("Downsampled to 16x16x8 voxels:\n")
cat("  Dimensions:", dim(vol_dims), "\n")
cat("  Spacing:", spacing(vol_dims), "mm\n")
cat("  Memory size:", format(object.size(vol_dims), units = "MB"), "\n\n")

# Non-uniform downsampling (different factors per dimension)
vol_nonuniform <- downsample(vol, factor = c(0.5, 0.5, 0.25))
cat("Non-uniform downsampling (0.5, 0.5, 0.25):\n")
cat("  Dimensions:", dim(vol_nonuniform), "\n")
cat("  Spacing:", spacing(vol_nonuniform), "mm\n")
cat("  Memory size:", format(object.size(vol_nonuniform), units = "MB"), "\n")

cat("\nNote: 3D downsampling reduces memory usage by the cubic factor of the downsampling.\n")
library(neuroim2)

# Create a simple 4D NeuroVec
# Dimensions: 3x3x3x3
data <- array(0, dim = c(3, 3, 3, 3))

# Fill with specific values to track the issue
# For time point 1: values 1-27
# For time point 2: values 101-127  
# For time point 3: values 201-227
for (t in 1:3) {
  offset <- (t-1) * 100
  for (k in 1:3) {
    for (j in 1:3) {
      for (i in 1:3) {
        linear_idx <- (k-1)*9 + (j-1)*3 + i
        data[i,j,k,t] <- linear_idx + offset
      }
    }
  }
}

# Print the values at coordinate [2,2,2] for all time points
cat("Direct array access data[2,2,2,]: ")
cat(data[2,2,2,], "\n")

# Create NeuroSpace and NeuroVec
space <- NeuroSpace(c(3, 3, 3, 3))
nv <- NeuroVec(data, space)

# Test series extraction at coordinate [2,2,2]
result <- series(nv, matrix(c(2,2,2), nrow=1))
cat("\nResult from series(nv, matrix(c(2,2,2), nrow=1)):\n")
cat("Length:", length(result), "\n")
cat("Values:", result, "\n")

# Expected: values at data[2,2,2,1:3] = 14, 114, 214
# Let's calculate the expected linear index
expected_idx <- (2-1)*9 + (2-1)*3 + 2  # = 9 + 3 + 2 = 14
cat("\nExpected linear index:", expected_idx, "\n")
cat("Expected values:", expected_idx, expected_idx+100, expected_idx+200, "\n")

# Let's trace through the series method implementation
cat("\n--- Debugging series method ---\n")
mat_input <- matrix(c(2,2,2), nrow=1)
cat("Input matrix:\n")
print(mat_input)

d4 <- dim(nv)[4]
cat("\nd4 (number of time points):", d4, "\n")

# The expanded matrix creation
expanded <- mat_input[rep(1:nrow(mat_input), each=d4),]
cat("\nExpanded coordinates (replicated d4 times):\n")
print(expanded)

# Add time indices
expanded <- cbind(expanded, 1:d4)
cat("\nExpanded with time indices:\n")
print(expanded)

# Extract values using the expanded indices
cat("\nExtracting values with nv[expanded]...\n")
vec <- nv[expanded]
cat("Extracted values:", vec, "\n")

# The final reshaping
result_matrix <- matrix(vec, d4, nrow(mat_input))
cat("\nFinal result matrix (d4 x nrow(i)):\n")
print(result_matrix)
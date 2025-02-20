#' Create a Memory-Mapped Neuroimaging Vector
#'
#' Creates a BigNeuroVec object, which represents a large neuroimaging vector using
#' memory-mapped file storage. This allows working with neuroimaging data that is
#' too large to fit in memory.
#'
#' @param data The input data to be stored
#' @param space A NeuroSpace object defining the spatial properties
#' @param mask A logical mask indicating which voxels contain data
#' @param label Optional character string label for the vector
#' @param type Storage type, one of "double", "float", or "integer"
#' @param backingfile Path to the file used for memory mapping (defaults to tempfile())
#'
#' @return A new BigNeuroVec object that provides memory-efficient access to large neuroimaging data through memory mapping.
#'         The object contains the spatial properties, mask, and memory-mapped data storage.
#'
#' @examples
#' \donttest{
#' # Load an example 4D brain image
#' example_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
#' example_4d_image <- read_vec(example_file)
#'
#' # Create a mask (e.g., selecting voxels with values > 0)
#' mask <- array(as.vector(example_4d_image[,,,1]) > 0,
#'              dim = dim(example_4d_image)[1:3])
#'
#' if(requireNamespace("bigstatsr", quietly = TRUE)) {
#'   # Create a BigNeuroVec with memory mapping
#'   big_vec <- BigNeuroVec(data = example_4d_image@.Data,
#'                          space = space(example_4d_image),
#'                          mask = mask,
#'                          label = "Example BigNeuroVec")
#'   print(big_vec)
#' }
#' }
#'
#' @rdname BigNeuroVec-methods
#' @export
BigNeuroVec <- function(data, space, mask, label = "", type = c("double", "float", "integer"), backingfile=tempfile()) {
  type <- match.arg(type)
  stopifnot(inherits(space, "NeuroSpace"))

  if (!requireNamespace("bigstatsr", quietly = TRUE)) {
    stop("Package 'bigstatsr' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  p <- prep_sparsenvec(data, space, mask)

  fbm <- bigstatsr::as_FBM(p$data, type=type, backingfile=backingfile)

  new("BigNeuroVec", space=p$space, mask=p$mask,
      map=IndexLookupVol(space(p$mask), as.integer(which(p$mask))), data=fbm, label=label)
}

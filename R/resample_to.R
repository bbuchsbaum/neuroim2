#' Resample an image with readable method names
#'
#' A convenience front-end to [resample()] that accepts human-friendly
#' method names and an engine switch. Internally delegates to the S4
#' `resample(source, target, interpolation = 0/1/3)` methods.
#'
#' @param source A `NeuroVol` (source image)
#' @param target A `NeuroVol` or `NeuroSpace` to match
#' @param method Interpolation method: `"nearest"`, `"linear"`, or `"cubic"`
#' @param engine Resampling engine. For now only `"internal"` is supported.
#' @param ... Reserved for future options
#' @return A `NeuroVol` in the target space
#' @examples
#' \donttest{
#' img <- read_vol(system.file("extdata","global_mask_v4.nii", package="neuroim2"))
#' sp  <- space(img)
#' sp2 <- NeuroSpace(sp@dim*2, sp@spacing/2, origin=sp@origin, trans=trans(img))
#' r1  <- resample_to(img, sp2, method = "linear")
#' }
#' @export
resample_to <- function(source, target,
                        method = c("nearest","linear","cubic"),
                        engine = c("internal"),
                        ...) {
  method <- match.arg(method)
  # Accept a broader set here so we can provide a clearer error below
  engine <- match.arg(engine, c("internal", "RNiftyReg"))

  if (engine != "internal") {
    stop("Only engine = 'internal' is supported at present.", call. = FALSE)
  }
  interp <- switch(method,
                   nearest = 0L,
                   linear  = 1L,
                   cubic   = 3L)
  resample(source, target, interpolation = interp)
}

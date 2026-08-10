# neuroim2 side of the nilearn comparison.
#
#     NEUROIM2_NL_DIR=/tmp/nl Rscript probe_r.R
#
# Writes res.json next to the inputs; compare.py diffs it against ref.json.
# Set NEUROIM2_NL_SCALING=1 to include the connected-component scaling probe,
# which takes about 15 minutes at the 1 mm size.

suppressMessages(library(neuroim2))

B <- Sys.getenv("NEUROIM2_NL_DIR", "/tmp/nl")
res <- list()

# `expr` must be re-evaluated on every repetition. Wrapping it in a closure
# would evaluate the promise once and then hand back the cached value, timing
# nothing; substitute/eval.parent forces the work each time.
tm <- function(expr, reps = 3) {
  invisible(eval.parent(substitute(expr)))
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(eval.parent(substitute(expr)))
    ts[i] <- proc.time()[["elapsed"]] - t0
  }
  median(ts)
}

vol <- read_vol(file.path(B, "vol.nii"))
run <- read_vec(file.path(B, "run.nii"))
mask <- read_vol(file.path(B, "mask.nii"))
labels <- read_vol(file.path(B, "labels.nii"))
truth <- as.array(mask) > 0

## ------------------------------------------------------------------ smoothing
## Impulse response, so the number reported is the smoothing actually applied
## rather than the argument that was asked for.
imp_sp <- NeuroSpace(c(41L, 41L, 41L), c(2, 2, 2))
imp <- array(0, c(41, 41, 41)); imp[21, 21, 21] <- 1
imp_vol <- NeuroVol(imp, imp_sp)
psf_fwhm <- function(o, spacing, centre) {
  a <- as.array(o); a[a < 0] <- 0
  prof <- apply(a, 1, sum)
  off <- (seq_along(prof) - centre) * spacing
  m <- sum(prof); mu <- sum(prof * off) / m
  2 * sqrt(2 * log(2)) * sqrt(sum(prof * (off - mu)^2) / m)
}
psf <- list()
for (sg in c(2, 3, 4)) {
  # The default derives the kernel half-width from sigma; window = 1 and 2 are
  # kept so the truncation this replaced stays visible and pinned.
  psf[[sprintf("sigma%g_default", sg)]] <- list(
    fwhm = psf_fwhm(gaussian_blur(imp_vol, sigma = sg, normalize = FALSE), 2, 21),
    requested = 2 * sqrt(2 * log(2)) * sg)
  for (w in 1:2) {
    psf[[sprintf("sigma%g_window%d", sg, w)]] <- list(
      fwhm = psf_fwhm(gaussian_blur(imp_vol, sigma = sg, window = w, normalize = FALSE), 2, 21),
      requested = 2 * sqrt(2 * log(2)) * sg)
  }
}
res$psf <- psf
# The default kernel now evaluates the same extent nilearn does, so this is a
# like-for-like comparison. `truncate = 4` matches scipy's default.
res$t_smooth_3d <- tm(gaussian_blur(vol, fwhm = 6))
res$t_smooth_4d <- tm(gaussian_blur(run, fwhm = 6), 1)

## ----------------------------------------------------------------- resampling
res$identity_resample_maxdiff <- max(abs(as.array(resample(vol, space(vol))) - as.array(vol)))

# A faithful downsample target: same world box and axis directions, 3 mm voxels.
tx <- trans(vol); tx[1:3, 1:3] <- tx[1:3, 1:3] %*% diag(rep(3 / 2, 3))
sp3 <- NeuroSpace(as.integer(ceiling(dim(vol) * 2 / 3)), c(3, 3, 3),
                  trans = tx, origin = tx[1:3, 4])
res$resample <- list()
for (kind in c("linear", "cubic")) {
  interp <- if (kind == "linear") 1L else 3L
  got <- as.vector(as.array(resample(vol, sp3, interpolation = interp)))
  ref <- as.vector(as.array(read_vol(file.path(B, sprintf("py_down_%s.nii", kind)))))
  res$resample[[kind]] <- list(cor = cor(got, ref),
                               rms = sqrt(mean((got - ref)^2)),
                               max_abs = max(abs(got - ref)))
}
res$t_resample <- tm(resample(vol, sp3, interpolation = 1L))

rl <- resample(labels, sp3, interpolation = 0L)
res$nearest_labels_in <- length(unique(as.vector(as.array(labels))))
res$nearest_labels_out <- length(unique(as.vector(as.array(rl))))

# The target space people reach for first, built from dims/spacing/origin alone.
# It mirrors an LAS source, so the result is near-empty; the point of the probe
# is that this is now said out loud rather than returned silently.
naive <- NeuroSpace(as.integer(ceiling(dim(vol) * 2 / 3)), c(3, 3, 3), origin = origin(vol))
warned <- FALSE
naive_out <- withCallingHandlers(
  resample(vol, naive),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
res$naive_target <- list(source_mean = mean(as.array(vol)),
                         resampled_mean = mean(as.array(naive_out)),
                         warned = warned)

## --------------------------------------------------------------- reorientation
odd <- read_vol(file.path(B, "odd.nii"))
cn <- as_canonical(odd)
a0 <- as.vector(as.array(odd)); a1 <- as.vector(as.array(cn))
res$reorder <- list(
  in_shape = dim(odd), out_shape = dim(cn),
  in_axcodes = paste(axcodes(space(odd)), collapse = ""),
  out_axcodes = paste(axcodes(space(cn)), collapse = ""),
  nonzero_in = sum(a0 != 0), nonzero_out = sum(a1 != 0),
  values_are_a_permutation = identical(sort(a0), sort(a1)),
  seconds = tm(as_canonical(odd), 1))

## --------------------------------------------------- connected components
res$cc <- list()
for (cn_ in c("6-connect", "18-connect", "26-connect")) {
  t0 <- proc.time()[["elapsed"]]
  cc <- suppressWarnings(conn_comp(vol, threshold = 110, connect = cn_))
  el <- proc.time()[["elapsed"]] - t0
  res$cc[[cn_]] <- list(n = num_clusters(cc$index),
                        largest = max(as.vector(as.array(cc$size))),
                        seconds = el)
}
res$cc_in_mask <- sum(as.array(vol) > 110)

if (nzchar(Sys.getenv("NEUROIM2_NL_SCALING"))) {
  res$cc_scaling <- list()
  for (nm in c("small", "typical", "hires")) {
    v <- read_vol(file.path(B, sprintf("scale_%s.nii", nm)))
    t0 <- proc.time()[["elapsed"]]
    cc <- suppressWarnings(conn_comp(v, threshold = 0.9, connect = "26-connect"))
    el <- proc.time()[["elapsed"]] - t0
    res$cc_scaling[[nm]] <- list(in_mask = sum(as.array(v) > 0.9),
                                 n = num_clusters(cc$index), seconds = el)
  }
}

## ------------------------------------------------------------------- masking
am <- automask(run)
aa <- as.array(am) > 0
res$epi_mask <- list(voxels = sum(aa),
                     dice = 2 * sum(aa & truth) / (sum(aa) + sum(truth)),
                     seconds = tm(automask(run), 1))

## --------------------------------------------- spheres, parcels, mask extract
cen <- as.matrix(utils::read.table(file.path(B, "sphere_centres_vox.txt")))
res$spheres <- list()
for (rad in c(6, 8)) {
  sphere_signal <- function() {
    vapply(seq_len(nrow(cen)), function(i) {
      roi <- spherical_roi(mask, cen[i, ], radius = rad, nonzero = TRUE)
      rowMeans(series(run, coords(roi)))
    }, numeric(dim(run)[4]))
  }
  sig <- sphere_signal()
  utils::write.table(sig, file.path(B, sprintf("r_spheres_%d.txt", rad)),
                     row.names = FALSE, col.names = FALSE)
  res$spheres[[as.character(rad)]] <- list(shape = dim(sig), seconds = tm(sphere_signal(), 1))
}

la <- as.array(labels)
parcel_signal <- function() {
  cl <- ClusteredNeuroVol(as.mask(labels), as.integer(la[la > 0]))
  ClusteredNeuroVec(run, cl, FUN = mean)@ts
}
ps <- parcel_signal()
utils::write.table(ps, file.path(B, "r_labels.txt"), row.names = FALSE, col.names = FALSE)
res$labels_signal <- list(shape = dim(ps), seconds = tm(parcel_signal(), 1))
res$masked_extract <- list(seconds = tm(series(run, which(truth)), 1))

## ------------------------------------------------ searchlight neighbourhoods
lv <- as.mask(mask)
res$searchlight <- list()
for (rad in c(6, 8)) {
  t0 <- proc.time()[["elapsed"]]
  sl <- searchlight_coords(lv, radius = rad, nonzero = TRUE)
  n <- vapply(seq_along(sl), function(i) nrow(sl[[i]]), integer(1))
  el <- proc.time()[["elapsed"]] - t0
  res$searchlight[[as.character(rad)]] <- list(centres = length(sl), mean = mean(n),
                                               median = as.integer(stats::median(n)),
                                               seconds = el)
}
# The two sibling iterators should agree about which voxels are centres.
res$centre_set <- list(searchlight = length(searchlight(lv, radius = 6)),
                       searchlight_coords = length(searchlight_coords(lv, radius = 6)),
                       in_mask = sum(as.array(lv)))

## ----------------------------------------------------------------- utilities
utils::write.table(as.vector(as.array(mean(run))), file.path(B, "r_mean.txt"),
                   row.names = FALSE, col.names = FALSE)
sub <- sub_vector(run, 1:20)
res$util <- list(mean_img = tm(mean(run)),
                 index_img = tm(run[[30]]),
                 threshold_img = tm(vol > 110),
                 binarize_img = tm(as.mask(vol > 110)),
                 math_img = tm(vol * 2 + 1),
                 concat_imgs = tm(concat(sub, sub, sub), 1))

## ---------------------------------------------------------------- write JSON
jnum <- function(v) {
  ifelse(is.finite(v), format(v, digits = 15),
         ifelse(is.nan(v) | is.na(v), "NaN", ifelse(v > 0, "Infinity", "-Infinity")))
}
to_json <- function(x) {
  if (is.list(x)) {
    if (!length(x)) return("{}")
    paste0("{", paste(sprintf('"%s":%s', names(x), vapply(x, to_json, character(1))),
                      collapse = ","), "}")
  } else if (is.logical(x) && length(x) == 1L) {
    if (isTRUE(x)) "true" else "false"
  } else if (is.character(x)) {
    if (length(x) == 1L) sprintf('"%s"', x)
    else paste0("[", paste(sprintf('"%s"', x), collapse = ","), "]")
  } else if (length(x) == 1L) {
    jnum(x)
  } else {
    paste0("[", paste(jnum(x), collapse = ","), "]")
  }
}
writeLines(to_json(res), file.path(B, "res.json"))
cat("wrote", file.path(B, "res.json"), "\n")

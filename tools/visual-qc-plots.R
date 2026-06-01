#!/usr/bin/env Rscript
# Visual QC for statistical-overlay plotting.
#
# Renders the scenarios exercised by the plotting fixes (issues #17/#18/#19 and
# enhance_stat_map) to labelled PNGs so a human (or agent) can eyeball that:
#   * signed maps show negatives as clearly as positives (diverging default),
#   * a sequential palette on signed data really does hide negatives (the bug),
#   * assemble=TRUE produces a single figure with an overlay colorbar,
#   * binary / proportional / ramp alpha modes look as intended,
#   * enhance_stat_map() de-speckles a salt-and-pepper map without flattening it.
#
# Usage:
#   Rscript tools/visual-qc-plots.R [output_dir]
# Default output_dir: tools/visual-qc-output

suppressMessages({
  if (file.exists("DESCRIPTION") && any(grepl("neuroim2", readLines("DESCRIPTION", n = 1)))) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(neuroim2)
  }
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "tools/visual-qc-output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(20240601)

# --- Realistic background + synthetic signed statistical map -----------------
anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
sp   <- space(anat)
d    <- dim(anat)

brain <- as.array(anat) > stats::quantile(as.array(anat)[as.array(anat) > 0], 0.20)
brain_idx <- which(brain)

# Signed t-like map: a positive cluster (left) and a negative cluster (right),
# on a unit-variance noise floor inside the brain.
make_signed <- function() {
  a <- array(0, d)
  a[brain_idx] <- stats::rnorm(length(brain_idx), 0, 1)
  a[12:19, 24:36, 16:32] <-  5   # positive
  a[30:37, 24:36, 16:32] <- -5   # negative
  a[!brain] <- 0
  NeuroVol(a, sp)
}

# Salt-and-pepper version for the enhance_stat_map demo.
make_noisy <- function() {
  a <- array(0, d)
  a[brain_idx] <- stats::rnorm(length(brain_idx), 0, 1)
  a[14:21, 26:34, 18:30] <- a[14:21, 26:34, 18:30] + 4   # real cluster
  spk <- sample(brain_idx, length(brain_idx) %/% 20)      # 5% impulses
  a[spk] <- 9 * sample(c(-1, 1), length(spk), replace = TRUE)
  a[!brain] <- 0
  NeuroVol(a, sp)
}

stat  <- make_signed()
noisy <- make_noisy()

zlevels <- round(seq(18, 30, length.out = 6))

save_overlay <- function(file, ..., w = 10, h = 3.6, dpi = 120) {
  dots <- list(...)
  page_bg <- if (identical(dots$style, "report")) "#f6f6f4" else "#141414"
  obj <- suppressWarnings(suppressMessages(plot_overlay(..., draw = FALSE)))
  path <- file.path(outdir, file)
  ggplot2::ggsave(path, obj, width = w, height = h, dpi = dpi, bg = page_bg)
  message("wrote ", path)
}

# 1. THE BUG: signed map forced through a sequential palette -> negatives vanish.
save_overlay("01_signed_sequential_BUG.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             bg_cmap = "grays", ov_cmap = "inferno", ov_thresh = 2.3,
             ov_symmetric = FALSE, style = "dark",
             title = "BUG: sequential 'inferno' on signed map (negatives hidden)")

# 2. THE FIX: default behavior auto-selects diverging + symmetric.
save_overlay("02_signed_diverging_default.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, style = "dark",
             title = "FIX: default diverging palette + symmetric limits + colorbar")

# 3. Explicit diverging palette (RdBu via hcl.pals, issue #17).
save_overlay("03_signed_RdBu.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             ov_cmap = "RdBu", ov_range = c(-6, 6), ov_thresh = 2.3, style = "dark",
             title = "RdBu (hcl.pals) with pinned ov_range = c(-6, 6)")

# 4-6. Alpha modes.
save_overlay("04_alpha_binary.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, ov_alpha_mode = "binary", style = "dark",
             title = "ov_alpha_mode = 'binary'")
save_overlay("05_alpha_proportional.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, ov_alpha_mode = "proportional", style = "dark",
             title = "ov_alpha_mode = 'proportional' (shared cap across panels)")
save_overlay("06_alpha_ramp.png",
             bgvol = anat, overlay = stat, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, ov_alpha_mode = "ramp", style = "dark",
             title = "ov_alpha_mode = 'ramp' (0 at threshold -> 1 at cap)")

# 7-8. enhance_stat_map de-speckling.
save_overlay("07_enhance_off_noisy.png",
             bgvol = anat, overlay = noisy, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, style = "dark",
             title = "Raw salt-and-pepper map (enhance = FALSE)")
save_overlay("08_enhance_on.png",
             bgvol = anat, overlay = noisy, zlevels = zlevels, ncol = 6,
             ov_thresh = 2.3, style = "dark", enhance = TRUE,
             title = "enhance = TRUE (despeckled, peaks preserved)")

# 9-10. 3x3 montages: default vs aggressive enhance, illustrating the
# "aggressiveness dial" (despike_k / radius / detail_gain + display threshold).
z9 <- round(seq(16, 32, length.out = 9))
save_overlay("09_montage_3x3_enhance_default.png",
             bgvol = anat, overlay = noisy, zlevels = z9, ncol = 3,
             ov_thresh = 2.3, style = "dark", enhance = TRUE,
             title = "3x3 montage, enhance = TRUE (default)",
             w = 9, h = 9.5)
save_overlay("10_montage_3x3_enhance_aggressive.png",
             bgvol = anat, overlay = noisy, zlevels = z9, ncol = 3,
             ov_thresh = 3, style = "dark",
             enhance = list(despike_k = 2.5, radius = 3, detail_gain = 0.8),
             title = "3x3 montage, enhance(despike_k=2.5, radius=3, detail_gain=0.8), thresh=3",
             w = 9, h = 9.5)

# 11. Publication "report" style: light card, cropped dark tiles, colorbar, legend.
save_overlay("11_report_style_3x3.png",
             bgvol = anat, overlay = noisy, zlevels = z9, ncol = 3,
             ov_thresh = 3, style = "report",
             enhance = list(despike_k = 2.5, radius = 3, detail_gain = 0.8),
             title = "3x3 Brain Activation Map",
             subtitle = "enhance(despike_k = 2.5, radius = 3, detail_gain = 0.8), thresh = 3",
             w = 9.5, h = 10)

message("\nVisual QC PNGs written to: ", normalizePath(outdir))

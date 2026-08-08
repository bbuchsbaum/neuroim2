## Capture / compare golden references for the functions being optimised.
## Usage:  Rscript golden.R capture <file>     Rscript golden.R compare <file>
suppressMessages(library(neuroim2))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; path <- args[2]
set.seed(20260808)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------------------------------------------------------------- fixtures
spacings <- list(c(1,1,1), c(2,2,2), c(1,1,4), c(0.8,0.8,3), c(3,1,2))
dims     <- list(c(12L,12L,12L), c(9L,14L,7L), c(20L,20L,20L))

mk_space <- function(d, sp) NeuroSpace(d, sp)
mk_vol <- function(d, sp, seed) {
  set.seed(seed)
  arr <- array(rnorm(prod(d)), d)
  arr[arr < -0.3] <- 0                       # make some voxels zero for nonzero=TRUE
  NeuroVol(arr, NeuroSpace(d, sp))
}

## centroids: corners, faces, interior, plus randoms
centroid_set <- function(d) {
  set.seed(7)
  fixed <- rbind(c(1,1,1), d, c(1, d[2], 1), c(d[1], 1, d[3]),
                 pmax(1L, d %/% 2L), c(2,2,2), c(d[1]-1, d[2]-1, d[3]-1))
  rnd <- cbind(sample.int(d[1], 6, TRUE), sample.int(d[2], 6, TRUE), sample.int(d[3], 6, TRUE))
  m <- rbind(fixed, rnd)
  storage.mode(m) <- "integer"
  unique(m)
}

## --------------------------------------------------- spherical_roi sweep
capture_spherical <- function() {
  out <- list()
  for (di in seq_along(dims)) for (si in seq_along(spacings)) {
    d <- dims[[di]]; sp <- spacings[[si]]
    vol <- mk_vol(d, sp, 100 + di*10 + si)
    spc <- mk_space(d, sp)
    cents <- centroid_set(d)
    for (radius in c(max(sp) + 1e-9, min(sp), 2.5, 4, 7, 13)) {
      if (radius < min(sp)) next
      for (nonzero in c(FALSE, TRUE)) for (use_cpp in c(TRUE, FALSE)) {
        for (target in c("space", "vol")) {
          bvol <- if (target == "space") spc else vol
          fill <- if (target == "space") NULL else NULL
          for (ci in seq_len(nrow(cents))) {
            key <- sprintf("sph|d%d|s%d|r%g|nz%d|cpp%d|%s|c%d",
                           di, si, radius, nonzero, use_cpp, target, ci)
            val <- tryCatch({
              r <- spherical_roi(bvol, cents[ci,], radius, fill = fill,
                                 nonzero = nonzero, use_cpp = use_cpp)
              list(ok = TRUE,
                   coords = unname(r@coords),
                   data = unname(as.numeric(r@.Data)),
                   ci = r@center_index, pi = r@parent_index,
                   cls = class(r)[1],
                   cmode = storage.mode(r@coords))
            }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
            out[[key]] <- val
          }
        }
      }
    }
  }
  ## explicit fill values
  d <- c(12L,12L,12L); sp <- c(2,2,2); spc <- mk_space(d, sp)
  for (f in list(1, 6, 0)) {
    r <- spherical_roi(spc, c(6L,6L,6L), 3.5, fill = f)
    out[[sprintf("sphfill|%g", f)]] <- list(ok=TRUE, coords=unname(r@coords),
      data=unname(as.numeric(r@.Data)), ci=r@center_index, pi=r@parent_index,
      cls=class(r)[1], cmode=storage.mode(r@coords))
  }
  ## matrix centroid form
  r <- spherical_roi(spc, matrix(c(6L,6L,6L), nrow=1), 3.5)
  out[["sphmat"]] <- list(ok=TRUE, coords=unname(r@coords), data=unname(as.numeric(r@.Data)),
                          ci=r@center_index, pi=r@parent_index, cls=class(r)[1],
                          cmode=storage.mode(r@coords))
  ## error paths
  for (nm in list(list(k="err_neg", c=c(0L,5L,5L), r=3.5),
                  list(k="err_oob", c=c(99L,5L,5L), r=3.5),
                  list(k="err_len", c=c(5L,5L), r=3.5),
                  list(k="err_smallr", c=c(5L,5L,5L), r=0.1))) {
    out[[nm$k]] <- tryCatch({ spherical_roi(spc, nm$c, nm$r); list(ok=TRUE, note="no error") },
                            error = function(e) list(ok=FALSE, msg=conditionMessage(e)))
  }
  out
}

## ------------------------------------------------ spherical_roi_set sweep
capture_roi_set <- function() {
  out <- list()
  d <- c(14L,11L,9L)
  for (si in seq_along(spacings)) {
    sp <- spacings[[si]]; spc <- mk_space(d, sp); vol <- mk_vol(d, sp, 500+si)
    cents <- centroid_set(d)
    for (radius in c(3, 6)) {
      if (radius < min(sp)) next
      for (nz in c(FALSE, TRUE)) for (target in c("space","vol")) {
        bvol <- if (target=="space") spc else vol
        rs <- spherical_roi_set(bvol, cents, radius, nonzero = nz)
        out[[sprintf("set|s%d|r%g|nz%d|%s", si, radius, nz, target)]] <-
          lapply(rs, function(r) list(coords=unname(r@coords), data=unname(as.numeric(r@.Data)),
                                      ci=r@center_index, pi=r@parent_index))
      }
    }
  }
  ## fill vector form
  spc <- mk_space(d, c(2,2,2)); cents <- centroid_set(d)
  rs <- spherical_roi_set(spc, cents, 4, fill = seq_len(nrow(cents)))
  out[["setfillvec"]] <- lapply(rs, function(r) list(coords=unname(r@coords),
    data=unname(as.numeric(r@.Data)), ci=r@center_index, pi=r@parent_index))
  rs <- spherical_roi_set(spc, cents, 4, fill = 9)
  out[["setfillscalar"]] <- lapply(rs, function(r) list(coords=unname(r@coords),
    data=unname(as.numeric(r@.Data)), ci=r@center_index, pi=r@parent_index))
  out
}

## ------------------------------------------------------- series() sweep
capture_series <- function() {
  out <- list()
  for (si in seq_along(spacings)) {
    d <- c(10L,9L,8L); sp <- spacings[[si]]; nt <- 7L
    set.seed(900+si)
    arr <- array(rnorm(prod(d)*nt), c(d, nt))
    dv <- NeuroVec(arr, NeuroSpace(c(d, nt), sp))
    set.seed(901+si)
    m <- array(runif(prod(d)) > 0.35, d)
    mask <- LogicalNeuroVol(m, NeuroSpace(d, sp))
    smat <- matrix(rnorm(sum(m)*nt), nt, sum(m))
    sv <- SparseNeuroVec(smat, NeuroSpace(c(d, nt), sp), mask)

    cds <- centroid_set(d)
    out[[sprintf("ser_dense_mat|%d", si)]] <- unname(series(dv, cds))
    out[[sprintf("ser_dense_int|%d", si)]] <- unname(series(dv, as.integer(c(1, 5, 77, prod(d)))))
    out[[sprintf("ser_dense_num|%d", si)]] <- unname(series(dv, c(3, 9, 100)))
    out[[sprintf("ser_dense_ijk|%d", si)]] <- unname(series(dv, 4L, 4L, 4L))
    out[[sprintf("ser_dense_1vox|%d", si)]] <- unname(series(dv, 42L))
    out[[sprintf("ser_dense_mask|%d", si)]] <- unname(series(dv, mask))
    out[[sprintf("ser_dense_nodrop|%d", si)]] <- unname(series(dv, 42L, drop = FALSE))

    out[[sprintf("ser_sp_mat|%d", si)]] <- unname(series(sv, cds))
    out[[sprintf("ser_sp_int|%d", si)]] <- unname(series(sv, as.integer(c(1, 5, 77, prod(d)))))
    out[[sprintf("ser_sp_num|%d", si)]] <- unname(series(sv, c(3, 9, 100)))
    out[[sprintf("ser_sp_ijk1|%d", si)]] <- unname(series(sv, 4L, 4L, 4L))
    out[[sprintf("ser_sp_ijkN|%d", si)]] <- unname(series(sv, c(1L,4L,7L), c(2L,4L,6L), c(3L,4L,5L)))
    out[[sprintf("ser_sp_1vox|%d", si)]] <- unname(series(sv, 42L))
    out[[sprintf("ser_sp_1vox_nodrop|%d", si)]] <- unname(series(sv, 42L, drop = FALSE))
    out[[sprintf("ser_sp_ijk1_nodrop|%d", si)]] <- unname(series(sv, 4L, 4L, 4L, drop = FALSE))
    out[[sprintf("ser_sp_mask|%d", si)]] <- unname(series(sv, mask))
    ## an all-out-of-mask request
    zero_idx <- which(!m)[1:5]
    out[[sprintf("ser_sp_allzero|%d", si)]] <- unname(series(sv, as.integer(zero_idx)))
    ## roi coord form
    roi <- spherical_roi(NeuroSpace(d, sp), c(5L,5L,4L), max(4, min(sp)))
    out[[sprintf("ser_dense_roi|%d", si)]] <- unname(series(dv, coords(roi)))
    out[[sprintf("ser_sp_roi|%d", si)]] <- unname(series(sv, coords(roi)))
    ## series_roi
    sr <- series_roi(dv, cds)
    out[[sprintf("ser_roi_obj|%d", si)]] <- list(d=unname(sr@.Data), c=unname(sr@coords))
    ## linear_access / [ on sparse
    out[[sprintf("sp_extract|%d", si)]] <- unname(sv[2:5, 3:6, 2:4, 1:3])
    out[[sprintf("sp_linacc|%d", si)]] <- unname(linear_access(sv, c(1, 500, 3000)))
    out[[sprintf("sp_asmat|%d", si)]] <- unname(as.matrix(sv))
  }
  out
}

## ------------------------------------------------------ searchlight sweep
capture_searchlight <- function() {
  out <- list()
  d <- c(10L,10L,10L)
  for (si in c(1,2,3)) {
    sp <- spacings[[si]]
    set.seed(300+si)
    m <- array(runif(prod(d)) > 0.3, d)
    mask <- LogicalNeuroVol(m, NeuroSpace(d, sp))
    for (nz in c(FALSE, TRUE)) {
      sl <- searchlight(mask, radius = max(4, min(sp)), eager = FALSE, nonzero = nz)
      out[[sprintf("sl_lazy|%d|nz%d", si, nz)]] <-
        lapply(seq_len(min(40, length(sl))), function(i) {
          r <- sl[[i]]; list(coords=unname(r@coords), ci=r@center_index, pi=r@parent_index,
                             mi=attr(r,"mask_index"), data=unname(as.numeric(r@.Data)))
        })
      sle <- searchlight(mask, radius = max(4, min(sp)), eager = TRUE, nonzero = nz)
      out[[sprintf("sl_eager|%d|nz%d", si, nz)]] <-
        lapply(seq_len(min(40, length(sle))), function(i) {
          r <- sle[[i]]; list(coords=unname(r@coords), ci=r@center_index, pi=r@parent_index,
                              mi=attr(r,"mask_index"), data=unname(as.numeric(r@.Data)))
        })
      slc <- searchlight_coords(mask, radius = max(4, min(sp)), nonzero = nz)
      out[[sprintf("sl_coords|%d|nz%d", si, nz)]] <-
        lapply(seq_len(min(40, length(slc))), function(i) unname(slc[[i]]))
    }
    set.seed(11); rsl <- random_searchlight(mask, radius = max(4, min(sp)))
    out[[sprintf("sl_random|%d", si)]] <-
      lapply(rsl, function(r) list(coords=unname(r@coords), ci=r@center_index,
                                   pi=r@parent_index, mi=attr(r,"mask_index")))
    set.seed(12); bsl <- resampled_searchlight(mask, radius = max(4, min(sp)), iter = 25)
    out[[sprintf("sl_resampled|%d", si)]] <-
      lapply(seq_len(25), function(i) { r <- bsl[[i]]
        list(coords=unname(r@coords), ci=r@center_index, pi=r@parent_index, mi=attr(r,"mask_index")) })
  }
  out
}

golden <- list(spherical = capture_spherical(),
               roi_set   = capture_roi_set(),
               series    = capture_series(),
               searchlight = capture_searchlight())

if (mode == "capture") {
  saveRDS(golden, path)
  cat(sprintf("captured: %s\n", paste(sprintf("%s=%d", names(golden), lengths(golden)), collapse=" ")))
} else {
  ref <- readRDS(path)
  bad <- 0L; checked <- 0L
  for (grp in names(ref)) {
    for (k in names(ref[[grp]])) {
      checked <- checked + 1L
      a <- ref[[grp]][[k]]; b <- golden[[grp]][[k]]
      cmp <- all.equal(a, b)
      if (!isTRUE(cmp)) {
        bad <- bad + 1L
        if (bad <= 25) cat(sprintf("MISMATCH [%s / %s]: %s\n", grp, k,
                                   paste(utils::head(cmp, 3), collapse=" | ")))
      }
    }
    extra <- setdiff(names(golden[[grp]]), names(ref[[grp]]))
    if (length(extra)) cat(sprintf("NOTE new keys in %s: %d\n", grp, length(extra)))
  }
  cat(sprintf("\n%d comparisons, %d mismatches\n", checked, bad))
  quit(status = if (bad > 0) 1 else 0)
}

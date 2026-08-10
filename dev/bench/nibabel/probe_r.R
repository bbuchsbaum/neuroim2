suppressMessages(library(neuroim2))
OUT <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "probes")

# Emits one JSON record per probe on stdout; compare.py does the comparison
# against nibabel's expected.json. Written by hand because jsonlite is not a
# dependency of this package.
files <- list.files(OUT, pattern = "\\.(nii|nii\\.gz|hdr)$", full.names = TRUE)
files <- files[!grepl("\\.img$", files)]

esc <- function(s) gsub('"', '\\\\"', s)
cat("[\n")
first <- TRUE
for (f in sort(files)) {
  name <- sub("\\.(nii\\.gz|nii|hdr)$", "", basename(f))
  res <- tryCatch({
    hdr <- read_header(f)
    nd <- length(dim(hdr))
    x <- if (nd >= 4) read_vec(f) else read_vol(f)
    dm <- dim(x)
    a <- as.numeric(as.vector(x@.Data))
    fin <- a[is.finite(a)]
    list(ok = TRUE,
         shape = dm,
         affine = as.vector(t(trans(x))),
         first = as.numeric(a[1]), second = as.numeric(a[2]),
         min = if (length(fin)) min(fin) else NA_real_,
         max = if (length(fin)) max(fin) else NA_real_,
         sum = if (length(fin)) sum(fin) else NA_real_,
         n_nan = sum(is.na(a) & !is.nan(a)) + sum(is.nan(a)),
         n_inf = sum(is.infinite(a)))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)),
     warning = function(w) list(ok = FALSE, msg = paste("WARNING:", conditionMessage(w))))

  if (!first) cat(",\n"); first <- FALSE
  if (isTRUE(res$ok)) {
    # Python's json module accepts the bare NaN / Infinity literals, so
    # non-finite values survive the round trip instead of becoming null.
    jnum <- function(v) {
      ifelse(is.finite(v), format(v, digits = 15),
             ifelse(is.na(v) & !is.nan(v), "NaN",
                    ifelse(is.nan(v), "NaN", ifelse(v > 0, "Infinity", "-Infinity"))))
    }
    num <- function(v) paste(jnum(v), collapse = ",")
    cat(sprintf('{"name":"%s","ok":true,"shape":[%s],"affine":[%s],"first":%s,"second":%s,"min":%s,"max":%s,"sum":%s,"n_nan":%d,"n_inf":%d}',
                name, paste(res$shape, collapse = ","), num(res$affine),
                num(res$first), num(res$second), num(res$min), num(res$max), num(res$sum),
                as.integer(res$n_nan), as.integer(res$n_inf)))
  } else {
    cat(sprintf('{"name":"%s","ok":false,"msg":"%s"}', name, esc(substr(res$msg, 1, 200))))
  }
}
cat("\n]\n")

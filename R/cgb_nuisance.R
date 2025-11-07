#' Build smooth time weights from motion/outlier metrics
#'
#' @description
#' Creates per-time-point weights \eqn{w_t \in [0, 1]} by smoothly combining
#' framewise displacement (FD), DVARS, and spike/outlier scores. Each series is
#' transformed through a soft logistic ramp so that values beyond the specified
#' thresholds receive progressively lower weights instead of hard 0/1 decisions.
#'
#' @param fd Optional numeric vector of framewise displacement values.
#' @param dvars Optional numeric vector of DVARS values.
#' @param spike Optional numeric vector with spike/outlier magnitudes.
#' @param fd_thr Threshold (in mm) where FD weights start to drop (default 0.5).
#' @param dvars_z Z-threshold applied to the standardized DVARS series (default 2.5).
#' @param spike_z Z-threshold applied to the standardized spike series (default 5).
#' @param fd_soft Logistic softness (in mm) controlling the slope around \code{fd_thr}.
#' @param dvars_soft Logistic softness for the DVARS z-scores.
#' @param combine Either \code{"min"} (take the minimum weight per TR) or
#'   \code{"prod"} (multiply all weights).
#'
#' @return Numeric vector of weights in \eqn{[0,1]} with length equal to the
#'   provided series. At least one of \code{fd}, \code{dvars}, or \code{spike}
#'   must be supplied.
#' @export
make_time_weights <- function(fd = NULL,
                              dvars = NULL,
                              spike = NULL,
                              fd_thr = 0.5,
                              dvars_z = 2.5,
                              spike_z = 5,
                              fd_soft = 0.1,
                              dvars_soft = 0.25,
                              combine = c("min", "prod")) {
  combine <- match.arg(combine)
  series_lengths <- c(length(fd), length(dvars), length(spike))
  series_lengths <- series_lengths[series_lengths > 0]
  if (!length(series_lengths)) {
    stop("Provide at least one of 'fd', 'dvars', or 'spike' to build weights.")
  }
  T_len <- unique(series_lengths)
  if (length(T_len) > 1) {
    stop("fd/dvars/spike must share the same length.")
  }
  T_len <- T_len[1]

  sigmoid <- function(x, thr, soft) {
    soft <- max(soft, .Machine$double.eps)
    1 / (1 + exp((x - thr) / soft))
  }

  pad_vec <- function(x) {
    if (length(x) == 0) {
      rep(NA_real_, T_len)
    } else {
      if (length(x) != T_len) stop("All series must have length ", T_len)
      as.numeric(x)
    }
  }

  fd_vals <- pad_vec(fd)
  dvars_vals <- pad_vec(dvars)
  spike_vals <- pad_vec(spike)

  if (all(is.na(fd_vals))) {
    w_fd <- rep(1, T_len)
  } else {
    fd_vals[!is.finite(fd_vals)] <- fd_thr
    w_fd <- sigmoid(fd_vals, fd_thr, fd_soft)
  }

  if (all(is.na(dvars_vals))) {
    w_dv <- rep(1, T_len)
  } else {
    dz <- as.numeric(base::scale(dvars_vals))
    dz[!is.finite(dz)] <- 0
    w_dv <- sigmoid(dz, dvars_z, dvars_soft)
  }

  if (all(is.na(spike_vals))) {
    w_sp <- rep(1, T_len)
  } else {
    sz <- as.numeric(base::scale(spike_vals))
    sz[!is.finite(sz)] <- 0
    w_sp <- 1 / (1 + exp((abs(sz) - spike_z) / 0.5))
  }

  w <- if (combine == "min") {
    pmin(w_fd, w_dv, w_sp)
  } else {
    w_fd * w_dv * w_sp
  }
  pmax(pmin(w, 1), 0)
}

#' Prepare weighted nuisance projectors for each run
#'
#' @description
#' Converts per-run confound matrices and time weights into orthonormal projectors
#' that can be consumed directly by the C++ graph builder. Each run produces
#' \eqn{Q_k} (columns spanning the weighted confound space) and \eqn{\sqrt{w_k}}
#' (per-time-point square-root weights). Supplying a \code{NULL} confound matrix
#' yields a zero-column projector, enabling pure time-weighting without
#' regression.
#'
#' @param confounds List of matrices (\eqn{T_k \times p_k}), or a single matrix
#'   reused for all runs. Each row corresponds to a time point.
#' @param time_weights Optional list of numeric vectors (length \eqn{T_k}) or a
#'   single vector reused for every run. If \code{NULL}, unit weights are used.
#' @param run_lengths Integer vector with the number of time points per run.
#'   Required when any run has both \code{confounds=NULL} and \code{time_weights=NULL}.
#' @param include_intercept Logical; prepend a column of ones before QR (default TRUE).
#' @param center_cols Logical; center each confound column before weighting (default TRUE).
#' @param scale_cols Logical; scale columns to unit variance before weighting (default FALSE).
#'
#' @return A list with elements \code{Q_list} (list of matrices) and
#'   \code{sqrtw_list} (list of numeric vectors). Each entry has the same length
#'   as \code{run_lengths}.
#' @export
prepare_confounds <- function(confounds,
                              time_weights = NULL,
                              run_lengths,
                              include_intercept = TRUE,
                              center_cols = TRUE,
                              scale_cols = FALSE) {
  if (missing(run_lengths)) {
    stop("'run_lengths' must be provided.")
  }
  run_lengths <- as.integer(run_lengths)
  K <- length(run_lengths)
  conf_list <- .cgb_normalize_confounds(confounds, run_lengths)
  weight_list <- .cgb_normalize_time_weights(time_weights, run_lengths)

  Q_list <- vector("list", K)
  sqrtw_list <- vector("list", K)

  for (k in seq_len(K)) {
    Tk <- run_lengths[k]
    wk <- weight_list[[k]]
    if (is.null(wk)) {
      wk <- rep(1, Tk)
    }
    wk[!is.finite(wk)] <- 0
    wk <- pmax(wk, 0)
    sw <- sqrt(wk)
    sqrtw_list[[k]] <- sw

    Ck <- conf_list[[k]]
    if (is.null(Ck)) {
      Ck <- matrix(numeric(0), nrow = Tk, ncol = 0)
    }
    if (!isTRUE(center_cols) || ncol(Ck) == 0) {
      centered <- Ck
    } else {
      centered <- base::scale(Ck, center = TRUE, scale = FALSE)
    }
    if (isTRUE(scale_cols) && ncol(centered) > 0) {
      centered <- base::scale(centered, center = FALSE, scale = TRUE)
    }
    if (include_intercept) {
      centered <- cbind(Intercept = 1, centered)
    }

    if (ncol(centered) == 0) {
      Q_list[[k]] <- matrix(0, nrow = Tk, ncol = 0)
    } else {
      Cw <- centered * sw
      qrobj <- qr(Cw, LAPACK = TRUE)
      rk <- qrobj$rank
      if (rk <= 0L) {
        Q_list[[k]] <- matrix(0, nrow = Tk, ncol = 0L)
      } else {
        Q_list[[k]] <- qr.Q(qrobj, complete = FALSE)[, seq_len(rk), drop = FALSE]
      }
    }
  }

  list(Q_list = Q_list, sqrtw_list = sqrtw_list)
}

# ---- Internal helpers ----------------------------------------------------

.cgb_normalize_time_weights <- function(time_weights, run_lengths) {
  K <- length(run_lengths)
  if (is.null(time_weights)) {
    return(vector("list", K))
  }
  if (!is.list(time_weights)) {
    if (K > 1) {
      stop("'time_weights' must be a list with one element per run.")
    }
    time_weights <- list(time_weights)
  }
  if (length(time_weights) == 1L && K > 1L) {
    time_weights <- rep(time_weights, K)
  }
  if (length(time_weights) != K) {
    stop("'time_weights' must have length ", K)
  }
  lapply(seq_len(K), function(k) {
    wk <- time_weights[[k]]
    if (is.null(wk)) {
      return(NULL)
    }
    wk <- as.numeric(wk)
    if (length(wk) != run_lengths[k]) {
      stop(sprintf("time_weights[[%d]] must have length %d", k, run_lengths[k]))
    }
    wk[!is.finite(wk)] <- 0
    pmax(pmin(wk, 1), 0)
  })
}

.cgb_normalize_confounds <- function(confounds, run_lengths) {
  K <- length(run_lengths)
  if (is.null(confounds)) {
    return(vector("list", K))
  }
  if (!is.list(confounds)) {
    if (K > 1) {
      stop("'confounds' must be a list with one entry per run.")
    }
    confounds <- list(confounds)
  }
  if (length(confounds) == 1L && K > 1L) {
    confounds <- rep(confounds, K)
  }
  if (length(confounds) != K) {
    stop("'confounds' must have length ", K)
  }
  lapply(seq_len(K), function(k) {
    ck <- confounds[[k]]
    if (is.null(ck)) {
      return(NULL)
    }
    ck <- as.matrix(ck)
    if (nrow(ck) != run_lengths[k]) {
      stop(sprintf("confounds[[%d]] must have %d rows", k, run_lengths[k]))
    }
    storage.mode(ck) <- "double"
    ck
  })
}

.cgb_compute_robust_weights <- function(arr,
                                        mask_idx,
                                        run_ends,
                                        method,
                                        c_value) {
  dims <- dim(arr)
  slice_xyz <- prod(dims[1:3])
  arr_vec <- as.numeric(arr)
  mask_idx <- as.integer(mask_idx)
  run_start <- c(1L, head(run_ends, -1L) + 1L)
  lapply(seq_along(run_ends), function(k) {
    start_t <- run_start[k]
    end_t <- run_ends[k]
    nt <- end_t - start_t + 1L
    if (nt <= 1L) {
      return(rep(1, nt))
    }
    prev <- arr_vec[mask_idx + (start_t - 1L) * slice_xyz]
    dvals <- numeric(nt)
    dvals[1] <- 0
    for (tt in 2:nt) {
      cur <- arr_vec[mask_idx + (start_t + tt - 1L) * slice_xyz]
      diff <- cur - prev
      valid <- is.finite(diff)
      if (any(valid)) {
        diff <- diff[valid]
        dvals[tt] <- sqrt(mean(diff * diff))
      } else {
        dvals[tt] <- 0
      }
      prev <- cur
    }
    .cgb_robust_profile(dvals, method, c_value)
  })
}

.cgb_robust_profile <- function(vals, method, c_value) {
  vals[!is.finite(vals)] <- 0
  med <- stats::median(vals, na.rm = TRUE)
  madv <- stats::mad(vals, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv <= 0) {
    madv <- stats::sd(vals, na.rm = TRUE)
    if (!is.finite(madv) || madv <= 0) {
      madv <- 1
    }
  }
  z <- (vals - med) / madv
  z[!is.finite(z)] <- 0
  if (method == "huber") {
    absz <- abs(z)
    w <- ifelse(absz <= c_value, 1, c_value / absz)
  } else {
    u <- abs(z) / c_value
    w <- ifelse(u < 1, (1 - u^2)^2, 0)
  }
  w[is.na(w)] <- 0
  w <- pmax(pmin(w, 1), 0)
  if (length(w) > 1) {
    w[1] <- max(w[1], w[2])
  } else {
    w[1] <- 1
  }
  w
}

#' Build a correlation-guided bilateral (CGB) graph
#'
#' @description
#' Computes a sparse row-stochastic graph whose weights combine spatial proximity
#' and pooled local time-series correlations. Supports optional censoring weights,
#' nuisance regression via weighted QR projectors, leave-one-run-out graph
#' construction, and robust down-weighting of high-DVARS volumes.
#'
#' @param runs A \code{\linkS4class{NeuroVec}} or a list of \code{NeuroVec}
#'   objects (typically one per run).
#' @param mask Optional \code{\linkS4class{LogicalNeuroVol}}/\code{NeuroVol} or
#'   logical array defining in-mask voxels. Defaults to all in-mask voxels.
#' @param window Integer half-width of the cubic spatial neighborhood
#'   (e.g., \code{1} yields a \(3\times3\times3\) window).
#' @param spatial_sigma Spatial Gaussian sigma in mm.
#' @param corr_map Mapping from pooled correlation to affinity; one of
#'   \code{"power"}, \code{"exp"}, or \code{"soft"}. The \code{"power"} and
#'   \code{"soft"} mappings rectify negative correlations, whereas
#'   \code{"exp"} preserves them (useful for sharpening more than smoothing).
#' @param corr_param Parameter for the chosen \code{corr_map}
#'   (gamma, tau, or r0 respectively).
#' @param topk Keep the strongest \code{k} neighbors after masking (0 keeps all).
#' @param leave_one_out Logical; if \code{TRUE} and multiple runs are provided,
#'   returns a list of graphs where run \code{u} excludes its own correlations.
#' @param run_weights Optional numeric weights per run used in Fisher-z pooling.
#'   Defaults to \eqn{n_k - 3} (usable frames minus three) when omitted.
#' @param add_self Logical; always inject a tiny self-edge before normalization.
#' @param time_weights Optional list (or single vector) of per-run time weights
#'   \eqn{w_t \in [0,1]} applied before correlation estimation. An intercept is
#'   always included so correlations are computed on weighted, demeaned series.
#' @param confounds Optional list (or single matrix) of per-run confound
#'   regressors to project out prior to correlation estimation.
#' @param robust One of \code{"none"}, \code{"huber"}, or \code{"tukey"}; when
#'   not \code{"none"} an additional DVARS-style reweighting is applied.
#' @param robust_c Tuning constant for the robust weights (Huber/Tukey).
#'
#' @return A list containing \code{row_ptr}, \code{col_ind}, \code{val},
#'   \code{dims3d}, and \code{mask_idx}, or (if \code{leave_one_out=TRUE})
#'   a list of such graphs.
#' @export
cgb_make_graph <- function(runs,
                           mask = NULL,
                           window = 1L,
                           spatial_sigma = 2,
                           corr_map = c("power", "exp", "soft"),
                           corr_param = 2,
                           topk = 16L,
                           leave_one_out = FALSE,
                           run_weights = NULL,
                           add_self = TRUE,
                           time_weights = NULL,
                           confounds = NULL,
                           robust = c("none", "huber", "tukey"),
                           robust_c = 1.345) {
  corr_map <- match.arg(corr_map)
  robust <- match.arg(robust)

  vec_list <- if (inherits(runs, "NeuroVec")) list(runs) else runs
  if (!is.list(vec_list) || !length(vec_list)) {
    stop("'runs' must be a NeuroVec or a non-empty list of NeuroVec objects.")
  }
  if (!all(vapply(vec_list, inherits, logical(1), "NeuroVec"))) {
    stop("All elements of 'runs' must be NeuroVec objects.")
  }
  K <- length(vec_list)

  dims_first <- dim(vec_list[[1]])
  if (length(dims_first) != 4) {
    stop("Each run must be a 4D NeuroVec.")
  }
  spatial_space <- drop_dim(space(vec_list[[1]]))
  spatial_dims <- dims_first[1:3]

  for (idx in seq_along(vec_list)) {
    cur <- vec_list[[idx]]
    if (!identical(spatial_dims, dim(cur)[1:3])) {
      stop("All runs must share identical spatial dimensions.")
    }
    if (!identical(spatial_space, drop_dim(space(cur)))) {
      stop("All runs must share the same spatial NeuroSpace (spacing/origin).")
    }
  }

  mask <- .cgb_prepare_mask(mask, spatial_space)
  mask_arr <- as.array(mask)
  if (!identical(dim(mask_arr), spatial_dims)) {
    stop("Mask dimensions do not match the spatial dimensions of 'runs'.")
  }
  mask_idx <- which(mask_arr != 0)
  if (!length(mask_idx)) {
    stop("Mask is empty; cannot build a graph.")
  }

  nt_each <- vapply(vec_list, function(x) dim(x)[4], integer(1))
  total_t <- sum(nt_each)
  arr <- array(0, dim = c(spatial_dims, total_t))
  cursor <- 0L
  for (idx in seq_along(vec_list)) {
    block <- as.array(vec_list[[idx]])
    block_len <- nt_each[idx]
    arr[,,, (cursor + seq_len(block_len))] <- block
    cursor <- cursor + block_len
  }
  storage.mode(arr) <- "double"

  run_ends <- cumsum(nt_each)
  spacing3 <- spacing(spatial_space)
  if (length(spacing3) != 3L) {
    stop("Spatial spacing must have length 3.")
  }

  corr_mode <- switch(corr_map,
                      power = 0L,
                      exp = 1L,
                      soft = 2L)
  if (!is.null(run_weights)) {
    if (length(run_weights) == 1L && K > 1L) {
      run_weights <- rep(run_weights, K)
    }
    if (length(run_weights) != K) {
      stop("'run_weights' must have one entry per run.")
    }
    run_weights <- as.numeric(run_weights)
  }

  needs_projection <- !is.null(confounds) ||
    !is.null(time_weights) ||
    robust != "none"

  build_fun <- if (needs_projection) {
    conf_list <- .cgb_normalize_confounds(confounds, nt_each)
    weight_list <- .cgb_normalize_time_weights(time_weights, nt_each)
    # If the only requested projection is an implicit intercept (i.e., no actual
    # confound columns provided, no time weights, and no robust reweighting),
    # then baseline correlation already handles mean-centering. In that case,
    # reuse the baseline builder to guarantee identical results.
    only_intercept <- (all(vapply(conf_list, function(m) is.null(m) || ncol(m) == 0L, logical(1))) &&
                       all(vapply(weight_list, is.null, logical(1))) &&
                       robust == "none")
    if (only_intercept) {
      function(loo_idx) {
        build_cgb_graph_cpp(
          arr,
          as.integer(mask_idx),
          as.integer(run_ends),
          as.integer(window),
          spatial_sigma,
          spacing3,
          as.integer(corr_mode),
          corr_param,
          as.integer(topk),
          as.integer(loo_idx),
          if (is.null(run_weights)) numeric() else run_weights,
          isTRUE(add_self)
        )
      }
    } else {
    if (robust != "none") {
      robust_list <- .cgb_compute_robust_weights(arr,
                                                 mask_idx,
                                                 run_ends,
                                                 robust,
                                                 robust_c)
      weight_list <- Map(function(base, extra) {
        if (is.null(base)) {
          extra
        } else {
          if (length(base) != length(extra)) {
            stop("Robust weighting length mismatch.")
          }
          base * extra
        }
      }, weight_list, robust_list)
    }
    prep <- prepare_confounds(conf_list,
                              time_weights = weight_list,
                              run_lengths = nt_each,
                              include_intercept = TRUE,
                              center_cols = TRUE,
                              scale_cols = FALSE)

    function(loo_idx) {
      build_cgb_graph_nuis_cpp(
        arr,
        as.integer(mask_idx),
        as.integer(run_ends),
        as.integer(window),
        spatial_sigma,
        spacing3,
        as.integer(corr_mode),
        corr_param,
        as.integer(topk),
        as.integer(loo_idx),
        if (is.null(run_weights)) numeric() else run_weights,
        prep$Q_list,
        prep$sqrtw_list,
        isTRUE(add_self)
      )
    }
    }
  } else {
    function(loo_idx) {
      build_cgb_graph_cpp(
        arr,
        as.integer(mask_idx),
        as.integer(run_ends),
        as.integer(window),
        spatial_sigma,
        spacing3,
        as.integer(corr_mode),
        corr_param,
        as.integer(topk),
        as.integer(loo_idx),
        if (is.null(run_weights)) numeric() else run_weights,
        isTRUE(add_self)
      )
    }
  }

  if (isTRUE(leave_one_out) && K > 1L) {
    lapply(seq_len(K), function(k) build_fun(k - 1L))
  } else {
    build_fun(-1L)
  }
}

#' Apply a precomputed CGB graph to volumetric data
#'
#' @param x A \code{\linkS4class{NeuroVec}} (4D) or \code{\linkS4class{NeuroVol}} (3D).
#' @param graph Graph list returned by \code{\link{cgb_make_graph}}.
#' @param passes Number of smoothing passes (>=1). Each pass multiplies by
#'   \code{W}; if \code{lambda < 1} a simple diffusion blend
#'   \code{(1 - lambda)I + lambda W} is applied per pass.
#' @param lambda Blend factor in \eqn{[0,1]} controlling diffusion strength.
#'
#' @return Smoothed object of the same class as \code{x}.
#' @export
cgb_smooth <- function(x, graph, passes = 1L, lambda = 1) {
  .cgb_validate_graph(graph)
  passes <- as.integer(passes)
  if (passes < 1L) stop("'passes' must be >= 1.")
  if (!is.numeric(lambda) || lambda < 0 || lambda > 1) {
    stop("'lambda' must lie in [0, 1].")
  }

  if (inherits(x, "NeuroVol")) {
    arr <- as.array(x)
    arr4 <- array(arr, dim = c(dim(arr), 1L))
    out <- apply_cgb_graph_cpp(arr4,
                               as.integer(graph$row_ptr),
                               as.integer(graph$col_ind),
                               graph$val,
                               as.integer(graph$mask_idx),
                               passes,
                               lambda)
    out <- array(out, dim = dim(arr))
    NeuroVol(out, space(x))
  } else if (inherits(x, "NeuroVec")) {
    arr <- as.array(x)
    out <- apply_cgb_graph_cpp(arr,
                               as.integer(graph$row_ptr),
                               as.integer(graph$col_ind),
                               graph$val,
                               as.integer(graph$mask_idx),
                               passes,
                               lambda)
    DenseNeuroVec(out, space(x))
  } else {
    stop("'x' must be a NeuroVol or NeuroVec.")
  }
}

#' Leave-one-run-out smoothing helper
#'
#' @param runs List of \code{\linkS4class{NeuroVec}} objects (one per run).
#' @param graphs List of graphs returned by \code{cgb_make_graph(..., leave_one_out=TRUE)}.
#' @param passes,lambda See \code{\link{cgb_smooth}}.
#'
#' @return A list of smoothed \code{NeuroVec} objects, one per run.
#' @export
cgb_smooth_loro <- function(runs, graphs, passes = 1L, lambda = 1) {
  if (!is.list(runs) || !is.list(graphs) || length(runs) != length(graphs)) {
    stop("'runs' and 'graphs' must be lists of equal length.")
  }
  Map(function(v, g) cgb_smooth(v, g, passes = passes, lambda = lambda),
      runs, graphs)
}

#' Correlation-guided bilateral filtering (convenience wrapper)
#'
#' @description
#' High-level interface that builds a correlation-guided bilateral (CGB) graph
#' with sensible defaults (similar to the bilateral filter interface) and
#' immediately applies it to smooth the data.
#'
#' @param runs A \code{\linkS4class{NeuroVec}} or a list of \code{NeuroVec}.
#' @param mask Optional \code{\linkS4class{LogicalNeuroVol}}/\code{NeuroVol} or
#'   logical array for spatial masking. Defaults to in-mask voxels.
#' @param spatial_sigma Spatial Gaussian sigma in mm. Used both for weighting
#'   and, when \code{window} is \code{NULL}, to auto-choose the neighborhood size.
#' @param window Integer half-width of the cubic neighborhood. If \code{NULL},
#'   it is computed as \code{ceiling(2 * spatial_sigma / min(spacing))} and at
#'   least 1.
#' @param corr_map Mapping from pooled correlation to affinity; one of
#'   \code{"power"}, \code{"exp"}, or \code{"soft"}. Defaults to \code{"power"}.
#' @param corr_param Parameter for \code{corr_map} (gamma/tau/r0 respectively).
#' @param topk Keep strongest \code{k} neighbors (0 keeps all). Defaults to 16.
#' @param passes Number of smoothing passes (>=1). Defaults to 1.
#' @param lambda Blend factor in [0,1] per pass. Defaults to 1 (pure diffusion).
#' @param leave_one_out If \code{TRUE} and multiple runs are supplied, builds
#'   LORO graphs and returns a list of smoothed runs.
#' @param run_weights Optional numeric weights per run for Fisher-z pooling.
#' @param add_self Logical; add a tiny self-edge before normalization.
#' @param time_weights Optional list (or single vector) of per-run time weights.
#' @param confounds Optional list (or single matrix) of per-run confounds.
#' @param robust One of \code{"none"}, \code{"huber"}, or \code{"tukey"}.
#' @param robust_c Tuning constant for robust weights.
#' @param return_graph Logical; if \code{TRUE}, also return the graph(s)
#'   alongside the smoothed data.
#'
#' @return If \code{leave_one_out=FALSE}, a smoothed \code{NeuroVec}. If
#'   \code{leave_one_out=TRUE}, a list of smoothed \code{NeuroVec}. When
#'   \code{return_graph=TRUE}, returns a list with elements \code{result} and
#'   \code{graph} (single object or lists accordingly).
#'
#' @examples
#' \donttest{
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Auto window from spatial_sigma and spacing, single pass
#' out <- cgb_filter(vec, mask, spatial_sigma = 3, window = NULL, topk = 16)
#'
#' # Stronger diffusion with two passes and lambda < 1
#' out2 <- cgb_filter(vec, mask, spatial_sigma = 4, window = NULL,
#'                    passes = 2, lambda = 0.7)
#' }
#'
#' @export
cgb_filter <- function(runs,
                       mask = NULL,
                       spatial_sigma = 2,
                       window = NULL,
                       corr_map = c("power", "exp", "soft"),
                       corr_param = 2,
                       topk = 16L,
                       passes = 1L,
                       lambda = 1,
                       leave_one_out = FALSE,
                       run_weights = NULL,
                       add_self = TRUE,
                       time_weights = NULL,
                       confounds = NULL,
                       robust = c("none", "huber", "tukey"),
                       robust_c = 1.345,
                       return_graph = FALSE) {

  corr_map <- match.arg(corr_map)
  robust <- match.arg(robust)

  # Normalize runs to list and extract spatial spacing
  vec_list <- if (inherits(runs, "NeuroVec")) list(runs) else runs
  if (!is.list(vec_list) || !length(vec_list)) {
    stop("'runs' must be a NeuroVec or a non-empty list of NeuroVec objects.")
  }
  if (!all(vapply(vec_list, inherits, logical(1), "NeuroVec"))) {
    stop("All elements of 'runs' must be NeuroVec objects.")
  }

  spatial_space <- drop_dim(space(vec_list[[1]]))
  spacing3 <- spacing(spatial_space)
  if (length(spacing3) != 3L) {
    stop("Spatial spacing must have length 3.")
  }

  # Auto-compute window if not provided
  if (is.null(window) || is.na(window)) {
    min_sp <- min(spacing3)
    window <- max(1L, as.integer(ceiling(2 * spatial_sigma / min_sp)))
  }

  graph <- cgb_make_graph(runs,
                          mask = mask,
                          window = as.integer(window),
                          spatial_sigma = spatial_sigma,
                          corr_map = corr_map,
                          corr_param = corr_param,
                          topk = as.integer(topk),
                          leave_one_out = leave_one_out,
                          run_weights = run_weights,
                          add_self = add_self,
                          time_weights = time_weights,
                          confounds = confounds,
                          robust = robust,
                          robust_c = robust_c)

  if (isTRUE(leave_one_out) && is.list(graph)) {
    # Multiple runs with LORO graphs
    res <- cgb_smooth_loro(vec_list, graph, passes = passes, lambda = lambda)
    if (isTRUE(return_graph)) {
      return(list(result = res, graph = graph))
    } else {
      return(res)
    }
  } else {
    # Single graph workflow
    if (inherits(runs, "NeuroVec")) {
      sm <- cgb_smooth(runs, graph, passes = passes, lambda = lambda)
    } else {
      # List input without LORO: build a single graph on all runs and
      # apply the same graph to each run.
      sm <- lapply(runs, function(v) cgb_smooth(v, graph, passes = passes, lambda = lambda))
    }
    if (isTRUE(return_graph)) {
      return(list(result = sm, graph = graph))
    } else {
      return(sm)
    }
  }
}

# ---- Internal helpers ----------------------------------------------------

.cgb_prepare_mask <- function(mask, spatial_space) {
  if (is.null(mask)) {
    LogicalNeuroVol(array(TRUE, dim(spatial_space)), spatial_space)
  } else if (inherits(mask, "LogicalNeuroVol") || inherits(mask, "NeuroVol")) {
    if (!identical(space(mask), spatial_space)) {
      stop("Mask space does not match the spatial space of 'runs'.")
    }
    mask
  } else if (is.array(mask)) {
    LogicalNeuroVol(mask, spatial_space)
  } else {
    stop("'mask' must be a LogicalNeuroVol, NeuroVol, or logical array.")
  }
}

.cgb_validate_graph <- function(graph) {
  required <- c("row_ptr", "col_ind", "val", "mask_idx")
  if (!is.list(graph) || !all(required %in% names(graph))) {
    stop("Graph must be a list returned by 'cgb_make_graph'.")
  }
}

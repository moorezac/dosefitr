#' Huber IRWLS
#' @description
#' Apply Iteratively Re-weighted Least Squares with a Huber loss to reduce
#'    the influence of outliers on model fits.
#' @param data Data frame containing the fitting data.
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#' @param args A named list of arguments to pass to `minpack.lm::nlsLM()`.
#' @param model_func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param k TODO: description.
#' @param outlier_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `huber`: Logical; use Huber robust regression? (default: FALSE).
#'     \item `huber_iter`: Maximum iterations for Huber (default: 100).
#'     \item `huber_k`: Huber tuning constant (default: 1.345).
#'     \item `huber_tol`: Convergence tolerance for Huber (default: 1e-6).
#'     \item `rout`: Logical; use ROUT outlier detection? (default: FALSE).
#'     \item `rout_q`: False discovery rate for ROUT (default: 1e-3).
#'     \item `rout_scale`: Scale estimator for ROUT, either "mad" or "quantile" (default: "mad").
#'   }
#' @return
#' Updated nls fit object with robust parameter estimates.
#' @details
#' This implementation is faster than performing a new fit at each iteration.
#' @keywords internal
huber_irwls <- function(
  obj,
  args,
  data,
  model_func,
  k = 1.345,
  outlier_opts = NULL
) {
  # Faster version
  if (is.null(outlier_opts)) {
    iter_max <- 100L
    k <- 1.345
    tol <- 1e-06
  } else {
    iter_max <- outlier_opts$huber_iter
    k <- outlier_opts$huber_k
    tol <- outlier_opts$huber_tol
  }

  # Extract initial parameters
  coef_curr <- coef(obj)
  coef_names <- names(coef_curr)

  m <- obj$m
  env <- m$getEnv()
  grad_fun <- m$gradient
  x <- get("x", env)
  y <- get("y", env)

  if ("group" %in% names(env)) {
    group <- get("group", env)
    group <- unique(group)
  } else {
    group <- NULL
  }

  huber_psi <- function(r, k) {
    # More robust scale estimation
    s <- median(abs(r - median(r))) * 1.4826
    # Ensure scale is never too small
    s <- max(s, 1e-6)

    u <- r / s
    w <- ifelse(abs(u) <= k, 1, k / abs(u))

    # Bound weights away from 0 and infinity
    w[!is.finite(w)] <- 1
    w <- pmax(w, 0.01) # Lower bound
    w <- pmin(w, 100) # Upper bound (more aggressive than 99th percentile)

    return(w)
  }

  damp <- 0.5

  for (i in seq_len(iter_max)) {
    # Compute residuals
    y_hat <- numeric(length(y))
    if (!is.null(group)) {
      for (j in group) {
        idx <- which(env$group == j)
        y_hat[idx] <- do.call(
          model_func,
          append(as.list(coef_curr), list(x = x[idx], group = j))
        )
      }
    } else {
      y_hat <- do.call(model_func, append(as.list(coef_curr), list(x = x)))
    }

    # This occurs more often without the damping step
    if (any(!is.finite(y_hat))) {
      break
    }

    res <- y - y_hat

    w <- huber_psi(res, k)

    X <- grad_fun()
    Xw <- X * sqrt(w)
    rw <- res * sqrt(w)

    XtX <- crossprod(Xw)
    Xtr <- crossprod(Xw, rw)

    lambda <- 1e-8
    step <- tryCatch(
      solve(XtX + diag(ncol(XtX)), Xtr),
      error = function(e) qr.solve(XtX + diag(ncol(XtX)), Xtr)
    )

    coef_new <- coef_curr + (damp * step)
    coef_new <- as.vector(coef_new)
    names(coef_new) <- coef_names

    # Check convergence
    value <- sqrt(sum((coef_new - coef_curr)^2))
    if (value < tol) {
      break
    }

    coef_curr <- coef_new
  }

  # Re-fit with calculated
  weighted_args <- args
  weighted_args$weights <- w

  # Re-fit
  new_fit <- attempt_fit(args = weighted_args, data = data, func = model_func)

  if (!is.null(new_fit$fit$error)) {
    warning("Error in Huber iteration, returning original object")
    return(obj)
  } else {
    return(new_fit$obj)
  }
}


#' Find Outlier Indices Using ROUT Method
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#' @param outlier_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `huber`: Logical; use Huber robust regression? (default: FALSE).
#'     \item `huber_iter`: Maximum iterations for Huber (default: 100).
#'     \item `huber_k`: Huber tuning constant (default: 1.345).
#'     \item `huber_tol`: Convergence tolerance for Huber (default: 1e-6).
#'     \item `rout`: Logical; use ROUT outlier detection? (default: FALSE).
#'     \item `rout_q`: False discovery rate for ROUT (default: 1e-3).
#'     \item `rout_scale`: Scale estimator for ROUT, either "mad" or "quantile" (default: "mad").
#'   }
#' @description
#' Identify outliers in model residuals using the ROUT (Robust regression and
#'    Outlier removal) method, which controls the False Discovery Rate (FDR).
#' @return
#' Numeric vector of outlier indices, or NULL if no outliers detected.
#' @details
#' The ROUT method uses a sequential testing procedure that controls the FDR
#'    at the specified level `q`. It's more appropriate than traditional methods
#'    when testing multiple points for outlier status.
#' @examples
#' \dontrun{
#' # Identify outliers in fitted model
#' outlier_indices <- find_rout_indices(
#'   obj = fitted_model,
#'   rout_opts = list(q = 0.01, scale_method = "mad")
#' )
#' }
#' @export
find_rout_indices <- function(
  obj,
  outlier_opts = NULL
) {
  if (is.null(outlier_opts)) {
    q <- 1e-3
    scale_method <- "mad"
  } else {
    q <- outlier_opts$rout_q
    scale_method <- outlier_opts$rout_scale
  }

  residuals_vec <- resid(obj)
  n <- length(residuals_vec)

  # Calculate scale factor
  if (scale_method == "mad") {
    scale_factor <- mad(residuals_vec)
  } else if (scale_method == "quantile") {
    scale_factor <- quantile(
      residuals_vec,
      probs = pnorm(1, 0, 1) - pnorm(-1, 0, 1)
    ) *
      (n / (n - length(coef(obj))))
  } else {
    stop("scale_method must be 'mad' or 'quantile'")
  }

  # Sort residuals and get indices
  abs_res_sorted <- sort(abs(residuals_vec), index.return = TRUE)$x
  indices_sorted <- sort(abs(residuals_vec), index.return = TRUE)$ix

  # Calculate FDR-adjusted p-values
  alphas <- q * seq(from = n, to = 1, by = -1) / n
  p_values <- 2 *
    pt(
      q = abs_res_sorted / scale_factor,
      df = n - length(coef(obj)),
      lower.tail = FALSE
    )

  # Find outliers
  indices_fdr <- which(p_values < alphas)

  if (length(indices_fdr) == 0) {
    return(NULL)
  } else {
    return(indices_sorted[seq(from = min(indices_fdr), to = n, by = 1)])
  }
}

#' Calculate Growth Rate Inhibition Metrics
#' @description
#' Computes GR (growth rate inhibition) values and optional NDR (normalised
#'    drug response) values from time-course data. Supports input as either raw
#'    data or pre-fitted model objects.
#' @param data A data frame containing time-course measurements. Required if
#'   \code{model_list} is NULL.
#' @param cap Logical; if TRUE (default), GR and NDR values are capped at 1.
#' @param return_all_columns Logical; if TRUE, returns all intermediate columns
#'   used in calculations. If FALSE (default), returns only original columns
#'   plus \code{gr} (and \code{ndr} if applicable).
#' @param model_list A named list of model results, typically the
#'   \code{model_list} element returned by \code{\link{fit_concentrations}}.
#'   Each element should contain \code{curve} (the fitted model) and
#'   \code{data} (the subset of data).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is \code{NA}.
#' @return
#' A data frame containing the original data plus:
#'  \itemize{
#'    \item `gr`: Growth rate inhibition values, ranging from -1 (complete cell
#'      death) through 0 (cytostasis) to 1 (no effect).
#'    \item `ndr`: Normalised drug response values (only if
#'      \code{positive_control_name} is provided), ranging from -1 to 1.
#'  }
#' If \code{return_all_columns = TRUE}, additional columns include:
#'  \code{time_zero_values}, \code{negative_control_zero}, \code{normalised},
#'  \code{negative_control_norm}, and (if applicable) \code{positive_control_zero}
#'   and \code{positive_control_norm}.
#' @details
#' GR values are calculated using the formula from Hafner et al. (2016):
#' \deqn{GR = 2^{log_2(x(c)/x_0) / log_2(x_{ctrl}/x_0)} - 1}
#' where \eqn{x(c)} is the measured value at concentration \eqn{c}, \eqn{x_0}
#' is the value at time zero, and \eqn{x_{ctrl}} is the negative control value.
#' NDR values provide an additional normalisation against a positive control,
#' useful when a maximal effect reference is available.
#' Groups marked with \code{exclude = TRUE} in the input data will have their
#' GR and NDR values set to -1.
#' @references
#' Hafner, M., Niepel, M., Chung, M., & Sorger, P. K. (2016). Growth rate
#' inhibition metrics correct for confounders in measuring sensitivity to
#' cancer drugs. \emph{Nature Methods}, 13(6), 521-527.
#' @export
calc_inhibition_metrics <- function(
  data = NA,
  model_list = NA,
  x_var,
  y_var,
  treatment_column,
  concentration_column,
  negative_control_name,
  positive_control_name = NA,
  cap = TRUE,
  return_all_columns = FALSE
) {
  if (!all(is.na(model_list))) {
    bound <- bind_model_data(model_list)

    x_range <- range(bound$data$x, na.rm = TRUE)

    predictions <- mapply(
      bound$nls_model_list,
      bound$data_list,
      SIMPLIFY = FALSE,
      FUN = function(nls_model, d) {
        predict_for_group(
          nls_model = nls_model,
          data = d,
          concentration_column = concentration_column,
          treatment_column = treatment_column,
          range = x_range,
          on_missing = "zero"
        )
      }
    )
    data <- do.call(rbind, predictions)
  } else {
    data <- prepare_data(data, x_var, y_var)
  }
  rownames(data) <- NULL

  # Store original columns for later selection
  original_columns <- colnames(data)

  data <- data[order(data$x), ]

  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )
  # Calculate time zero values for each treatment-concentration group
  group_id <- paste(
    data[[treatment_column]],
    data[[concentration_column]],
    sep = "_"
  )
  data$time_zero_values <- stats::ave(
    data$y,
    group_id,
    FUN = function(y) rep(y[1], length(y))
  )

  has_positive <- !is.na(positive_control_name)
  negative_control_mask <- data[[treatment_column]] == negative_control_name
  positive_control_mask <- if (has_positive) {
    data[[treatment_column]] == positive_control_name
  } else {
    rep(FALSE, nrow(data))
  }

  data$negative_control_zero <- control_time_zero(
    data$time_zero_values,
    negative_control_mask,
    "Negative control group"
  )
  if (has_positive) {
    data$positive_control_zero <- control_time_zero(
      data$time_zero_values,
      positive_control_mask,
      "Positive control group"
    )
  }

  # Calculate normalized values
  data$normalised <- data$y / data$time_zero_values

  # Get control normalized values for each time point
  data$negative_control_norm <- control_norm_by_time(
    data$normalised,
    data$x,
    negative_control_mask
  )
  if (has_positive) {
    data$positive_control_norm <- control_norm_by_time(
      data$normalised,
      data$x,
      positive_control_mask
    )
  }

  # Remove last time point per group
  keep_rows <- unlist(lapply(
    split(seq_len(nrow(data)), group_id),
    function(idx) {
      if (length(idx) > 1) idx[-length(idx)] else integer(0)
    }
  ))
  data <- data[sort(keep_rows), , drop = FALSE]

  data$gr <- 2^(log2(data$normalised) / log2(data$negative_control_norm)) - 1
  if (cap) {
    data$gr[data$gr > 1] <- 1
  }

  if (has_positive) {
    data$ndr <- (1 -
      2^(log2(data$normalised) / log2(data$positive_control_norm))) /
      (1 -
        2^(log2(data$negative_control_norm) / log2(data$positive_control_norm)))
    data$ndr[data$ndr > 1] <- 1
    data$ndr[data$ndr < -1] <- -1
  }

  rownames(data) <- NULL

  if ("exclude" %in% colnames(data)) {
    mask_exclude <- data$exclude %in% TRUE
    data$gr[mask_exclude] <- -1
    if ("ndr" %in% colnames(data)) {
      data$ndr[mask_exclude] <- -1
    }
  }

  if (return_all_columns) {
    return(data)
  }

  cols_to_select <- c(original_columns, "gr")
  if (has_positive) {
    cols_to_select <- c(cols_to_select, "ndr")
  }

  data[, cols_to_select, drop = FALSE]
}

#' Calculate All Drug Response Metrics
#' @description
#' Computes summary drug response metrics (GR50, LGR score) for one or more
#'   inhibition metrics (GR and optionally NDR) from a list of fitted model
#'   objects. For each metric, fits a dose-response curve, extracts the
#'   half-maximal inhibitory concentration, and calculates a goodness-of-fit
#'   adjusted LGR (Longitudinal Growth Rate) score.
#' @param model_list A named list of model results, typically the
#'   \code{model_list} element returned by \code{\link{fit_concentrations}}.
#'   Each element should contain \code{curve} (the fitted model) and
#'   \code{data} (the subset of data).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR-based metrics are also calculated. Default is \code{NA}.
#' @return
#' A named list containing summary metrics for each inhibition metric computed.
#' For GR (and NDR if applicable), the list includes:
#' \itemize{
#'   \item \code{gr}: GR50 value as a formatted concentration string, or
#'     \code{NA} if it could not be found.
#'   \item \code{lgr_gr}: LGR score derived from GR values.
#'   \item \code{lgr_adj_gr}: Goodness-of-fit adjusted LGR score derived from
#'     GR values.
#'   \item \code{ndr}, \code{lgr_ndr}, \code{lgr_adj_ndr}: The NDR
#'     equivalents (only if \code{positive_control_name} is provided).
#' }
#' @details
#' For each inhibition metric, the function:
#' \enumerate{
#'   \item Calls \code{\link{calc_inhibition_metrics}} to compute per-timepoint
#'     GR or NDR values from the model predictions.
#'   \item Fits a five-parameter sigmoid dose-response curve to the metric
#'     values across concentrations.
#'   \item Identifies the concentration at which the metric equals 0.5
#'     (GR50 or NDR50) using root-finding on the fitted curve.
#'   \item Computes per-concentration area under the curve (AUC) using the
#'     trapezoidal rule, then derives a normalised area over the curve (AOC)
#'     as a percentage of the total bounding area.
#'   \item Fits a second dose-response curve to the AOC values across
#'     concentrations and integrates it to yield a scalar LGR score.
#'   \item Adjusts the LGR score by the goodness-of-fit (R²) of the AOC
#'     curve fit, penalising poorly constrained curves.
#' }
#' The LGR score captures the overall magnitude and shape of the drug response
#' across the full concentration range, rather than relying on a single summary
#' statistic. Goodness-of-fit adjustment follows the approach used in the
#' DPRL pipeline.
#' @seealso \code{\link{calc_inhibition_metrics}}, \code{\link{fit_curve}},
#'   \code{\link{auc_trapezoid}}
#' @references
#' Hafner, M., Niepel, M., Chung, M., & Sorger, P. K. (2016). Growth rate
#' inhibition metrics correct for confounders in measuring sensitivity to
#' cancer drugs. \emph{Nature Methods}, 13(6), 521-527.
#' @export
calc_all_metrics <- function(
  model_list,
  x_var,
  concentration_column,
  treatment_column,
  negative_control_name,
  positive_control_name = NA
) {
  metrics <- if (is.na(positive_control_name)) "gr" else c("gr", "ndr")

  metric_list <- list()

  for (metric in metrics) {
    data <- prepare_metric_data(
      model_list = model_list,
      metric = metric,
      x_var = x_var,
      concentration_column = concentration_column,
      treatment_column = treatment_column,
      negative_control_name = negative_control_name,
      positive_control_name = positive_control_name,
      treatment_only = TRUE
    )$data

    # Fit on a copy with x/y set explicitly: calculate_lgr() below needs `data`
    # with x still holding time
    fit_data <- data
    fit_data$x <- fit_data[[concentration_column]]
    fit_data$y <- fit_data[[metric]]

    fitted_model <- fit_curve(
      data = fit_data,
      x_var = concentration_column,
      y_var = metric,
      curve_options = list(
        model = "five_param_sigmoid_log",
        lower_bounds = list(bottom = -1),
        upper_bounds = list(top = 1)
      )
    )

    predicted <- predict_data(
      nls_model = fitted_model$curve$nls_model,
      lower_x = min(fitted_model$data$x),
      upper_x = max(fitted_model$data$x)
    )

    metric_50 <- try(
      find_x_for_y(
        nls_model = fitted_model$curve$nls_model,
        x_values = predicted$x,
        target = 0.5
      ),
      silent = TRUE
    )

    metric_list[[metric]] <- if (inherits(metric_50, "try-error")) {
      NA_character_
    } else {
      metric_50$root
    }

    lgr <- calculate_lgr(
      data = data,
      metric = metric,
      concentration_column = concentration_column
    )

    metric_list[[paste("lgr", metric, sep = "_")]] <- lgr$lgr
    metric_list[[paste("lgr_adj", metric, sep = "_")]] <- lgr$lgr_adj
  }

  metric_list
}


#' Calculate GR/NDR values from a model list and prepare them for plotting
#'
#' @param model_list,x_var,concentration_column,treatment_column Passed to
#'   [calc_inhibition_metrics()].
#' @param metric `"gr"` or `"ndr"`.
#' @param negative_control_name,positive_control_name Control labels.
#' @param treatment_only Keep only treatment rows (drop controls).
#' @return A list with `data` (finite `metric` values, with a `condition`
#'   column) and `info` (from [extract_concentration_info()], taken before
#'   filtering).
#' @keywords internal
prepare_metric_data <- function(
  model_list,
  metric,
  x_var,
  concentration_column,
  treatment_column,
  negative_control_name,
  positive_control_name = NA,
  treatment_only = FALSE
) {
  metric <- match.arg(metric, c("gr", "ndr"))

  if (identical(metric, "ndr") && is.na(positive_control_name)) {
    stop("'ndr' requires a positive control; none was supplied.")
  }

  data <- calc_inhibition_metrics(
    model_list = model_list,
    x_var = x_var,
    y_var = metric,
    concentration_column = concentration_column,
    treatment_column = treatment_column,
    negative_control_name = negative_control_name,
    positive_control_name = positive_control_name
  )

  # calc_inhibition_metrics() assigns conditions internally but drops the
  # column unless return_all_columns = TRUE, so re-attach it here.
  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )

  # Concentration metadata is taken before filtering, so the palette and
  # labels still cover controls and excluded doses.
  info <- extract_concentration_info(
    data,
    concentration_column,
    positive_control_name
  )

  # Downstream dose-response fits group by concentration only, so more than
  # one drug would be silently pooled into a single curve
  treatments <- unique(as.character(
    data[[treatment_column]][data$condition == "treatment"]
  ))
  if (length(treatments) > 1) {
    stop(sprintf(
      "Expected one treatment (besides controls) but found %d: %s. Subset `model_list` to one treatment at a time.",
      length(treatments),
      paste(treatments, collapse = ", ")
    ))
  }

  data <- data[is.finite(data[[metric]]), , drop = FALSE]

  if (treatment_only) {
    data <- data[data$condition == "treatment", , drop = FALSE]
  }

  if (!nrow(data)) {
    stop(sprintf(
      "No finite '%s' values remain%s.",
      metric,
      if (treatment_only) " for treatment conditions" else ""
    ))
  }

  rownames(data) <- NULL

  list(data = data, info = info)
}

#' Calculate the LGR score from per-timepoint GR/NDR values
#'
#' @param data Data with `x` (time), `metric` and `concentration_column`
#'   columns.
#' @param metric Name of the metric column (`"gr"` or `"ndr"`).
#' @param concentration_column Column of log10 molar concentrations.
#' @param curve_options Curve options for the AOC dose-response fit, passed to
#'   [fit_curve()].
#' @param gof_floor Lower limit applied to the goodness of fit when adjusting
#'   the score.
#' @return A list with `lgr`, `lgr_adj`, `gof`, `aoc_data`, `predicted` and
#'   `fit` (the [fit_curve()] result).
#' @keywords internal
calculate_lgr <- function(
  data,
  metric,
  concentration_column,
  curve_options = list(
    model = "five_param_sigmoid_log",
    lower_bounds = list(bottom = 0),
    upper_bounds = list(top = 100)
  ),
  gof_floor = 0
) {
  missing_cols <- setdiff(c("x", metric, concentration_column), names(data))
  if (length(missing_cols)) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  keep <- is.finite(data$x) &
    is.finite(data[[metric]]) &
    is.finite(data[[concentration_column]])
  data <- data[keep, , drop = FALSE]

  if (!nrow(data)) {
    stop(sprintf(
      "No finite '%s' values to calculate an LGR score from.",
      metric
    ))
  }

  rng <- range(data[[metric]])
  if (rng[1] < -1 || rng[2] > 1) {
    warning(
      sprintf(
        "'%s' values outside [-1, 1] (observed %.3g to %.3g); AOC will fall outside 0-100%%.",
        metric,
        rng[1],
        rng[2]
      ),
      call. = FALSE
    )
  }

  concs <- sort(unique(data[[concentration_column]]))

  # One bounding box for every dose, so scores are comparable across them
  total_time <- diff(range(data$x))
  if (total_time <= 0) {
    stop("At least two distinct time points are needed to calculate LGR.")
  }
  bounding_area <- 2 * total_time

  auc <- vapply(
    concs,
    function(conc) {
      d <- data[data[[concentration_column]] == conc, , drop = FALSE]
      d <- d[order(d$x), , drop = FALSE]
      auc_trapezoid(d$x, d[[metric]] + 1)
    },
    numeric(1)
  )

  aoc_data <- data.frame(
    concentration = concs,
    auc = auc / bounding_area * 100,
    aoc = (bounding_area - auc) / bounding_area * 100
  )

  fitted_model <- fit_curve(
    data = aoc_data,
    x_var = "concentration",
    # Change to aoc for original implementation
    y_var = "auc",
    curve_options = curve_options
  )

  predicted <- predict_data(
    nls_model = fitted_model$curve$nls_model,
    lower_x = min(fitted_model$data$x),
    upper_x = max(fitted_model$data$x)
  )

  conc_span <- diff(range(aoc_data$concentration))
  lgr <- auc_trapezoid(predicted$x, predicted$y) / conc_span

  # Compare against the rows that were actually fitted
  gof <- calculate_gof(fitted_model$curve$nls_model, fitted_model$data$y)

  list(
    lgr = lgr,
    lgr_adj = lgr * max(gof, gof_floor, na.rm = TRUE),
    gof = gof,
    aoc_data = aoc_data,
    predicted = predicted,
    fit = fitted_model
  )
}

#' Area under a curve by the trapezoidal rule
#'
#' @param x,y Numeric vectors of equal length. Points are sorted by `x`
#'   first.
#' @return The area, or `0` for fewer than two points.
#' @keywords internal
auc_trapezoid <- function(x, y) {
  if (length(x) != length(y)) {
    stop("`x` and `y` must have equal lengths.")
  }
  if (length(x) < 2L) {
    return(0)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum((y[-length(y)] + y[-1L]) / 2 * diff(x))
}

#' Goodness of fit (R squared) of a fitted model
#'
#' @param nls_model The fitted model.
#' @param observed The observed response values the model was fitted to.
#' @return `1 - SS_res / SS_tot`, or `NA` if `observed` has no variance.
#' @keywords internal
calculate_gof <- function(nls_model, observed) {
  ss_res <- sum(stats::residuals(nls_model)^2)
  ss_tot <- sum((observed - mean(observed))^2)

  if (!is.finite(ss_tot) || ss_tot == 0) {
    return(NA_real_)
  }

  1 - (ss_res / ss_tot)
}

#' Time-zero value of a control group
#'
#' @param time_zero_values Per-row time-zero values.
#' @param mask Logical vector selecting the control rows.
#' @param label Group name used in the error message.
#' @return The single time-zero value of the control.
#' @keywords internal
control_time_zero <- function(time_zero_values, mask, label) {
  value <- unique(time_zero_values[mask %in% TRUE])
  if (length(value) != 1) {
    stop(sprintf("%s must have a unique time zero value.", label))
  }

  value
}

#' Mean normalised control value at each x
#'
#' @param normalised Per-row normalised values.
#' @param x Per-row x (time) values.
#' @param mask Logical vector selecting the control rows.
#' @return A vector the length of `x` holding the control mean at each row's
#'   x value (`NA` where the control has no value at that x).
#' @keywords internal
control_norm_by_time <- function(normalised, x, mask) {
  mask <- mask %in% TRUE
  if (!any(mask)) {
    stop("No control rows found; cannot normalise.")
  }

  # Match on the x values themselves: converting them to names and back with
  # as.numeric() does not round-trip every double exactly.
  control_x <- unique(x[mask])
  means <- tapply(
    normalised[mask],
    match(x[mask], control_x),
    mean,
    na.rm = TRUE
  )

  as.numeric(means)[match(x, control_x)]
}

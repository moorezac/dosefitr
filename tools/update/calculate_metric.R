#' Calculate dose-response summary metrics (GR/NDR, LGR score)
#' Fits dose-response curves for growth rate (GR) and, if a positive control
#' is supplied, normalized drug response (NDR), returning IC50-style
#' potency values and an area-under-curve-derived LGR score for each metric.
#' @param model_list List of fitted per-well/per-condition models, as
#'   produced upstream (e.g. by the growth curve fitting pipeline).
#' @param x_var Character. Name of the x-axis variable (e.g. time) used
#'   when fitting the underlying growth models.
#' @param concentration_column Character. Name of the column containing
#'   treatment concentrations.
#' @param treatment_column Character. Name of the column identifying
#'   treatment vs. control wells/conditions.
#' @param negative_control_name Value in `treatment_column` identifying
#'   negative control wells.
#' @param positive_control_name Value in `treatment_column` identifying
#'   positive control wells. If `NA` (default), only the GR metric is
#'   calculated and NDR is skipped.
#' @return A named list with, for each metric (`gr` and optionally `ndr`):
#'   \item{<metric>}{Estimated concentration (as a formatted string) at
#'     which the metric crosses 0.5.}
#'   \item{lgr_<metric>}{GOF-adjusted LGR score summarizing the
#'     dose-response curve.}
#' @export
calc_all_metrics <- function(
  model_list,
  x_var,
  concentration_column,
  treatment_column,
  negative_control_name,
  positive_control_name = NA
) {
  if (!is.na(positive_control_name)) {
    metrics <- c("gr", "ndr")
  } else {
    metrics <- c("gr")
  }

  metric_list <- list()

  for (metric in metrics) {
    data <- calc_inhibition_metrics(
      model_list = model_list,
      x_var = x_var,
      y_var = metric,
      concentration_column = concentration_column,
      treatment_column = treatment_column,
      negative_control_name = negative_control_name,
      positive_control_name = positive_control_name
    )

    data <- assign_condition(
      data,
      treatment_column,
      negative_control_name,
      positive_control_name
    )

    data <- data[is.finite(data[[metric]]), ]

    concentration_info <- extract_concentration_info(
      data,
      concentration_column,
      positive_control_name
    )

    data <- data[is.finite(data$gr), ]
    data <- data[data$condition == "treatment", ]

    fitted_model <- fit_curve(
      data = data,
      x_var = concentration_column,
      y_var = metric,
      curve_opts = list(
        model = "five_param_sigmoid_log",
        lower_bounds = list(bottom = -1),
        upper_bounds = list(top = 1)
      )
    )

    predicted <- predict_data(
      obj = fitted_model$fit$obj,
      lower_x = min(fitted_model$data$x),
      upper_x = max(fitted_model$data$x)
    )

    metric_50 <- try(
      find_x_for_y(
        obj = fitted_model$fit$obj,
        x_values = predicted$x,
        target = 0.5
      ),
      silent = TRUE
    )

    metric_list[[metric]] <- log_m_to_str(metric_50$root)

    auc_values <- lapply(unique(data[[concentration_column]]), function(a) {
      data_filt <- data[data[[concentration_column]] == a, ]
      auc_trapezoid(data_filt$x, data_filt[[metric]] + 1)
    })

    auc_data <- data.frame(
      concentration = unique(data[[concentration_column]]),
      auc = unlist(auc_values)
    )

    total_time <- max(data$x, na.rm = TRUE) - min(data$x, na.rm = TRUE)

    total_bounding_area <- 2 * total_time
    aoc_values <- total_bounding_area - unlist(auc_values)
    aoc_values_norm <- aoc_values / total_bounding_area * 100

    auc_data <- data.frame(
      aoc = aoc_values_norm,
      concentration = unique(data[[concentration_column]])
    )

    fitted_model <- fit_curve(
      data = auc_data,
      x_var = "concentration",
      y_var = "aoc",
      curve_opts = list(
        model = "five_param_sigmoid_log"
      )
    )

    predicted <- predict_data(
      obj = fitted_model$fit$obj,
      lower_x = min(fitted_model$data$x),
      upper_x = max(fitted_model$data$x)
    )

    auc <- auc_trapezoid(predicted$x, predicted$y)

    lgr_score <- auc /
      (max(auc_data$concentration) - min(auc_data$concentration))

    # GOF - this is what dprl uses
    ss_res <- sum(residuals(fitted_model$fit$obj)^2)
    ss_total <- sum((auc_data$aoc - mean(auc_data$aoc))^2)
    gof <- 1 - (ss_res / ss_total)

    lgr_score_adjusted <- lgr_score * gof

    lgr_str <- paste("lgr", metric, sep = "_")
    metric_list[[lgr_str]] <- lgr_score_adjusted
  }

  return(metric_list)
}

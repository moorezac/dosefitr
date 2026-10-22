#' Find Identical Columns
#' @description
#' Find identical columns within a dataframe.
#' @param data A dataframe.
#' @param target A character for which to match for.
#' @param coerc_to_char Logical; convert all columns to characters before match?
#' @param coerce_to_char TODO: description.
#' @keywords internal
find_identical_to_column <- function(data, target, coerce_to_char = TRUE) {
  # Data to match for
  data_target <- data[[target]]
  if (is.character(target) && coerce_to_char) {
    data <- data.frame(
      lapply(data, as.character),
      stringsAsFactors = FALSE
    )
  }
  identical_cols <- vapply(
    colnames(data),
    function(x) {
      if (x != target) {
        # Exclude the target column itself
        identical(data[[x]], data_target)
      } else {
        FALSE
      }
    },
    FUN.VALUE = logical(1)
  )
  colnames(data)[identical_cols]
}

#' Plot Fitted Models
#' @description
#' Create publication-quality plots showing fitted model curves overlaid on
#'    experimental data points, with support for grouped models and outliers.
#' @param list A list as produced by `fit_curve()`.
#' @param return_data Logical indicating whether to return plot data instead of
#'   the plot object (default: FALSE).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @return
#' A ggplot object (default) or a list containing original and predicted
#'   data for external plotting (if return_data = TRUE).
#' @importFrom ggplot2 geom_line geom_point ggplot labs aes guides guide_legend
#' @examples
#' \dontrun{
#' # Plot single model fit
#' model_plot <- plot_model(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "cell_count",
#'   obj = fitted_model
#' )
#' print(model_plot)
#' # Return data for custom plotting
#' plot_data <- plot_model(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "cell_count",
#'   fit = fitted_model,
#'   return_data = TRUE
#' )
#' }
#' @export
plot_model <- function(list, x_var, y_var, return_data = FALSE) {
  data <- list[["data"]]
  fit <- list[["fit"]]

  if (!all(c("x", "y", "x_original", "y_original") %in% colnames(data))) {
    data <- prep_data(data, x_var, y_var)
  }

  # Check for grouping
  formula_vars <- all.vars(formula(fit$obj))
  if ("group" %in% formula_vars) {
    group_vec <- unique(data$group)
    predicted <- predict_data(
      fit$obj,
      min(data$x),
      max(data$x),
      group = group_vec
    )
    use_colour <- TRUE

    original_group <- find_identical_to_column(
      data = data,
      target = "group",
      coerce_to_char = TRUE
    )
  } else {
    predicted <- predict_data(fit$obj, min(data$x), max(data$x))
    use_colour <- FALSE
  }

  # Check for outliers
  outlier_column <- paste("outlier", y_var, sep = "_")

  if (outlier_column %in% colnames(data)) {
    use_shape <- TRUE
  } else {
    use_shape <- FALSE
  }

  # Aesthetics
  map_point <- map_line <- ggplot2::aes(x = x, y = y)
  if (use_colour) {
    map_point$colour <- as.name("group")
    map_line$colour <- as.name("group")
  }
  if (use_shape) {
    map_point$shape <- as.name(outlier_column)
  }

  # If wanted to plot elsewhere
  if (return_data) {
    return(list(data = data, predicted = predicted))
  }

  # Plot
  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data, mapping = map_point, alpha = 0.25) +
    ggplot2::geom_line(data = predicted, mapping = map_line) +
    ggplot2::labs(x = x_var, y = y_var) +
    {
      if ("group" %in% formula_vars) {
        ggplot2::guides(colour = ggplot2::guide_legend(original_group))
      }
    } +
    {
      if (outlier_column %in% colnames(data)) {
        ggplot2::guides(shape = ggplot2::guide_legend("Outlier"))
      }
    } +
    theme_dosefitr()

  return(gg)
}


#' Assign Condition Column Based on Controls
#' @param data Data frame containing treatment information.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @return
#' The input data frame with a new \code{condition} column.
#' @keywords internal
assign_condition <- function(
  data,
  treatment_column,
  negative_control_name,
  positive_control_name = NA
) {
  if (is.na(positive_control_name)) {
    data$condition <- ifelse(
      data[[treatment_column]] == negative_control_name,
      "negative_control",
      "treatment"
    )
  } else {
    data$condition <- ifelse(
      data[[treatment_column]] == negative_control_name,
      "negative_control",
      ifelse(
        data[[treatment_column]] == positive_control_name,
        "positive_control",
        "treatment"
      )
    )
  }
  return(data)
}


#' Generate Concentration Colour Palette
#' @param unique_concs Numeric vector of unique concentrations.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @param start_colour Hex colour string for the highest concentration. Lower
#'   concentrations are displayed as progressively lighter tints. Default is
#'   \code{"#0085ca"}.
#' @param exclude_concs Numeric vector of concentrations to mark as excluded.
#' @return
#' A named character vector of colours.
#' @importFrom tinter tinter
#' @keywords internal
make_concentration_palette <- function(
  unique_concs,
  positive_control_name = NA,
  start_colour = "#0085ca",
  exclude_concs = NA
) {
  total_concs <- length(unique_concs)

  # Generate tints
  if (total_concs > 1) {
    steps <- floor(total_concs / 2) + 1
  } else {
    steps <- 2
  }

  tints <- tinter::tinter(
    x = start_colour,
    steps = steps,
    adjust = -0.1
  )
  tints <- rev(tints[seq_len(total_concs)])
  names(tints) <- paste(unique_concs, "treatment", sep = "_")

  # Mark excluded concentrations
  if (!any(is.na(exclude_concs))) {
    exclude_names <- paste(exclude_concs, "treatment", sep = "_")
    tints[names(tints) %in% exclude_names] <- "#ff0000"
  }

  # Add controls
  tints <- c("negative_control" = "#000000", tints)
  if (!is.na(positive_control_name)) {
    tints <- c(tints, "positive_control" = "#ff0000")
  }

  return(tints)
}


#' Generate Concentration Labels for Plots
#' @param negative_control_conc Concentration value for negative control.
#'   Default is NA.
#' @param positive_control_conc Concentration value for positive control.
#'   Default is NA.
#' @param unique_concs Numeric vector of unique concentrations.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @return
#' Character vector of formatted labels.
#' @keywords internal
make_concentration_labels <- function(
  unique_concs,
  negative_control_name,
  negative_control_conc = NA,
  positive_control_name = NA,
  positive_control_conc = NA
) {
  labels <- log_m_to_str(sort(unique_concs))

  # Add negative control label at front
  if (is.na(negative_control_conc) || negative_control_conc == 0) {
    neg_label <- negative_control_name
  } else {
    neg_label <- paste(
      negative_control_name,
      log_m_to_str(negative_control_conc)
    )
  }
  labels <- c(neg_label, labels)

  # Add positive control label at end
  if (!is.na(positive_control_name)) {
    if (is.na(positive_control_conc) || positive_control_conc == 0) {
      pos_label <- positive_control_name
    } else {
      pos_label <- paste(
        positive_control_name,
        log_m_to_str(positive_control_conc)
      )
    }
    labels <- c(labels, pos_label)
  }

  return(labels)
}


#' Assign Colour Vector for Concentration Plotting
#' @param data A data frame with a \code{condition} column.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @return
#' The input data frame with a new \code{colour_vec} column.
#' @keywords internal
assign_colour_vec <- function(data, concentration_column) {
  data$colour_vec <- ifelse(
    data$condition == "treatment",
    paste(data[[concentration_column]], "treatment", sep = "_"),
    data$condition
  )
  return(data)
}

#' Helper function to extract unique concentrations and control concentrations
#' from a processed data frame.
#' @param data A data frame with a \code{condition} column.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @return A list containing:
#'   \itemize{
#'     \item `unique_concs`: Sorted numeric vector of treatment concentrations.
#'     \item `negative_control_conc`: Concentration value for negative control.
#'     \item `positive_control_conc`: Concentration value for positive control,
#'       or NA if not applicable.
#'     \item `exclude_concs`: Concentration values marked as excluded, or NA.
#'   }
#' @keywords internal
extract_concentration_info <- function(
  data,
  concentration_column,
  positive_control_name = NA
) {
  unique_concs <- sort(unique(
    data[[concentration_column]][data$condition == "treatment"]
  ))

  negative_control_conc <- unique(
    data[[concentration_column]][data$condition == "negative_control"]
  )
  if (!length(negative_control_conc)) {
    negative_control_conc <- NA
  }

  if (!is.na(positive_control_name)) {
    positive_control_conc <-
      unique(data[[concentration_column]][data$condition == "positive_control"])
    if (!length(positive_control_conc)) {
      positive_control_conc <- NA
    }
  } else {
    positive_control_conc <- NA
  }

  exclude_concs <- NA
  if ("exclude" %in% colnames(data)) {
    exclude_concs <- unique(data[[concentration_column]][data$exclude])
  }

  return(list(
    unique_concs = unique_concs,
    negative_control_conc = negative_control_conc,
    positive_control_conc = positive_control_conc,
    exclude_concs = exclude_concs
  ))
}


#' Plot Fitted Curves Across Concentrations
#' Creates a publication-quality plot showing fitted growth curves for each
#' treatment-concentration combination, with concentration-based colour scaling.
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
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
#'   NDR values are calculated. Default is NULL.
#' @param start_colour Hex colour string for the highest concentration. Lower
#'   concentrations are displayed as progressively lighter tints. Default is
#'   \code{"#0085ca"}.
#' @return
#' A \code{ggplot2} object showing observed data as semi-transparent
#'   points and fitted curves as solid lines, coloured by concentration.
#' @seealso \code{\link{interpolate_curve_concentration}},
#'   \code{\link{calc_inhibition_metrics}}, \code{\link{theme_dosefitr}}
#' @examples
#' \dontrun{
#' p <- plot_curve_concentration(
#'   model_list = result$model_list,
#'   x_var = "time",
#'   y_var = "confluence",
#'   treatment_column = "treatment",
#'   concentration_column = "dose_uM",
#'   negative_control_name = "DMSO"
#' )
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_colour_manual scale_fill_manual labs guides guide_legend
#' @export
plot_curve_concentration <- function(
  model_list,
  x_var,
  y_var,
  treatment_column,
  concentration_column,
  negative_control_name,
  positive_control_name = NA,
  start_colour = "#0085ca"
) {
  obj_list <- lapply(model_list, function(x) {
    x$fit$obj
  })
  data_list <- lapply(model_list, function(x) {
    x$data
  })
  all_cols <- Reduce(union, lapply(data_list, names))
  data_list <- lapply(data_list, function(d) {
    missing <- setdiff(all_cols, names(d))
    for (m in missing) {
      d[[m]] <- NA
    }
    d[all_cols]
  })
  data <- do.call(rbind, data_list)
  rownames(data) <- NULL

  # Generate predictions
  prediction_list <- mapply(
    obj_list,
    data_list,
    FUN = function(obj, data) {
      if (all(is.na(obj))) {
        return(NA)
      }

      pred <- predict_data(
        obj = obj,
        lower_x = min(data$x),
        upper_x = max(data$x)
      )
      pred[[concentration_column]] <- unique(data[[concentration_column]])
      pred[[treatment_column]] <- unique(data[[treatment_column]])

      return(pred)
    },
    SIMPLIFY = FALSE
  )
  predicted <- do.call(rbind, prediction_list)
  rownames(predicted) <- NULL

  # Assign conditions
  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )
  predicted <- assign_condition(
    predicted,
    treatment_column,
    negative_control_name,
    positive_control_name
  )

  concentration_info <- extract_concentration_info(
    data,
    concentration_column,
    positive_control_name
  )

  # Generate palette and labels
  tints <- make_concentration_palette(
    unique_concs = concentration_info$unique_concs,
    positive_control_name = positive_control_name,
    start_colour = start_colour,
    exclude_concs = concentration_info$exclude_concs
  )
  labels <- make_concentration_labels(
    unique_concs = concentration_info$unique_concs,
    negative_control_name = negative_control_name,
    negative_control_conc = concentration_info$negative_control_conc,
    positive_control_name = positive_control_name,
    positive_control_conc = concentration_info$positive_control_conc
  )

  # Create colour vectors
  data <- assign_colour_vec(data, concentration_column)
  predicted <- assign_colour_vec(predicted, concentration_column)

  # Remove excluded data points
  if ("exclude" %in% colnames(data)) {
    data <- data[!data$exclude, ]
  }
  predicted <- predicted[!is.na(predicted$x), ]
  predicted <- predicted[!is.na(predicted$y), ]

  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data,
      mapping = ggplot2::aes(x, y, fill = colour_vec, colour = colour_vec),
      alpha = 0.25,
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE
    ) +
    ggplot2::geom_line(
      predicted,
      mapping = ggplot2::aes(x, y, colour = colour_vec, group = colour_vec),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE,
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = tints,
      limits = names(tints),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(concentration_column)
    ) +
    ggplot2::scale_fill_manual(
      values = tints,
      limits = names(tints),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(concentration_column)
    ) +
    ggplot2::labs(x = x_var, y = y_var) +
    theme_dosefitr()

  return(gg)
}


#' Plot Growth Rate Inhibition Metric Over Time
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param metric Character string specifying the response metric to use,
#'   either \code{"gr"} for growth rate inhibition or \code{"ndr"} for
#'   normalised drug response. Default is \code{"gr"}.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @param start_colour Hex colour string for the highest concentration. Lower
#'   concentrations are displayed as progressively lighter tints. Default is
#'   \code{"#0085ca"}.
#' @description
#' Visualises GR or NDR values over time for each concentration, showing the
#' temporal dynamics of drug response.
#' @return
#' A list containing:
#'  \itemize{
#'    \item `data`: The processed data frame with calculated metrics.
#'    \item `plot`: A \code{ggplot2} object showing metric values over time,
#'      with a dashed line at y = 0 (cytostasis).
#'  }
#' @seealso \code{\link{calc_inhibition_metrics}},
#'   \code{\link{interpolate_curve_concentration}}
#' @importFrom ggplot2 ggplot aes geom_hline geom_line scale_colour_manual scale_fill_manual ylim labs guides guide_legend
#' @export
plot_curve_concentration_metric_time <- function(
  model_list,
  x_var,
  metric = "gr",
  concentration_column,
  treatment_column,
  negative_control_name,
  positive_control_name = NA,
  start_colour = "#0085ca"
) {
  data <- calc_inhibition_metrics(
    model_list = model_list,
    x_var = x_var,
    y_var = metric,
    concentration_column = concentration_column,
    treatment_column = treatment_column,
    negative_control_name = negative_control_name,
    positive_control_name = positive_control_name
  )

  data$y <- data[[metric]]

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

  # Generate palette and labels
  tints <- make_concentration_palette(
    unique_concs = concentration_info$unique_concs,
    positive_control_name = positive_control_name,
    start_colour = start_colour,
    exclude_concs = concentration_info$exclude_concs
  )
  labels <- make_concentration_labels(
    unique_concs = concentration_info$unique_concs,
    negative_control_name = negative_control_name,
    negative_control_conc = concentration_info$negative_control_conc,
    positive_control_name = positive_control_name,
    positive_control_conc = concentration_info$positive_control_conc
  )

  # Create colour vector
  data <- assign_colour_vec(data, concentration_column)

  gg <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(
      data,
      mapping = ggplot2::aes(x, y, colour = colour_vec),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE,
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = tints,
      limits = names(tints),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(title = concentration_column)
    ) +
    ggplot2::scale_fill_manual(
      values = tints,
      limits = names(tints),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(title = concentration_column)
    ) +
    ggplot2::ylim(c(-1, 1)) +
    ggplot2::labs(x = x_var, y = metric) +
    theme_dosefitr()

  return(list(data = data, plot = gg))
}


#' Plot Dose-Response Curve for Growth Rate Metric
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param metric Character string specifying the response metric to use,
#'   either \code{"gr"} for growth rate inhibition or \code{"ndr"} for
#'   normalised drug response. Default is \code{"gr"}.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @param curve_opts A named list of curve fitting options. Elements include:
#'   \itemize{
#'     \item `model`: Character string specifying a built-in model from \code{dosefitr_models} or "loess".
#'     \item `model_func`: A function that describes a curve in terms of x (used if \code{model} is NA).
#'     \item `start_func`: A function that generates a named list of starting values for `model_func` (used if \code{model} is NA).
#'     \item `start_values`: A named list of starting values for `model_func`.
#'     \item `lower_bounds`: A named list that specifies the lower bounds for parameters in `model_func`.
#'     \item `upper_bounds`:A named list that specifies the upper bounds for parameters in `model_func`.
#'   }
#' @param start_colour Hex colour string for the highest concentration. Lower
#'   concentrations are displayed as progressively lighter tints. Default is
#'   \code{"#0085ca"}.
#' @param lower_bounds TODO: description.
#' @description
#' Fits and visualises a dose-response curve showing the relationship between
#' drug concentration and growth rate inhibition metric, with time indicated
#' by point colour.
#' @return
#' A list containing:
#'   \describe{
#'     \item{data}{The processed data frame with calculated metrics}
#'     \item{plot}{A \code{ggplot2} object showing points coloured by timepoint,
#'       fitted dose-response curve, vertical line at IC50/GR50 (if calculable),
#'       and dashed horizontal line at y = 0 (cytostasis). The x-axis is
#'       formatted for log-molar concentrations.}
#'   }
#' @seealso \code{\link{calc_inhibition_metrics}}, \code{\link{fit_curve}},
#'   \code{\link{theme_dosefitr}}
#' @examples
#' \dontrun{
#' result <- plot_curve_concentration_metric_dose(
#'   model_list = fitted_models$model_list,
#'   x_var = "time",
#'   metric = "gr",
#'   concentration_column = "dose_uM",
#'   treatment_column = "treatment",
#'   negative_control_name = "DMSO"
#' )
#' result$plot
#' }
#' @importFrom ggplot2 ggplot aes geom_hline geom_point geom_line scale_x_continuous scale_colour_viridis_c ylim labs guides
#' @importFrom geomtextpath geom_textvline
#' @export
plot_curve_concentration_metric_dose <- function(
  model_list,
  x_var,
  metric = "gr",
  concentration_column,
  treatment_column,
  negative_control_name,
  positive_control_name = NA,
  curve_opts = list(
    model = "five_param_sigmoid_log",
    lower_bounds = list(bottom = -1),
    upper_bounds = list(top = 1)
  ),
  start_colour = "#0085ca"
) {
  data <- calc_inhibition_metrics(
    model_list = model_list,
    x_var = x_var,
    y_var = metric,
    concentration_column = concentration_column,
    treatment_column = treatment_column,
    negative_control_name = negative_control_name,
    positive_control_name = positive_control_name
  )

  data$time <- data$x
  data$x <- data[[concentration_column]]
  data$y <- data[[metric]]

  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )

  data <- data[is.finite(data[[metric]]), ]
  data <- data[data$condition == "treatment", ]

  fitted_model <- fit_curve(
    data = data,
    x_var = concentration_column,
    y_var = metric,
    curve_opts = curve_opts
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

  gg <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(
      data,
      mapping = ggplot2::aes(x, y, colour = time),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE
    ) +
    ggplot2::geom_line(
      predicted,
      mapping = ggplot2::aes(x, y),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) log_m_to_str(x),
      minor_breaks = log10(rep(1:9, 21) * (10^rep(-10:10, each = 9)))
    ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::ylim(c(-1, 1)) +
    ggplot2::labs(x = concentration_column, y = metric) +
    ggplot2::guides(colour = "none") +
    theme_dosefitr()

  if (!inherits(metric_50, "try-error")) {
    gg <- gg +
      geomtextpath::geom_textvline(
        xintercept = metric_50$root,
        label = log_m_to_str(metric_50$root),
        linewidth = 0.5,
        hjust = 0.1,
        size = 10 / ggplot2::.pt
      )
  }

  return(list(data = data, plot = gg))
}


#' Plot Longitudinal Growth Rate Score with Dose-Response Curve
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param metric Character string specifying the response metric to use,
#'   either \code{"gr"} for growth rate inhibition or \code{"ndr"} for
#'   normalised drug response. Default is \code{"gr"}.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param curve_opts A named list of curve fitting options. Elements include:
#'   \itemize{
#'     \item `model`: Character string specifying a built-in model from \code{dosefitr_models} or "loess".
#'     \item `model_func`: A function that describes a curve in terms of x (used if \code{model} is NA).
#'     \item `start_func`: A function that generates a named list of starting values for `model_func` (used if \code{model} is NA).
#'     \item `start_values`: A named list of starting values for `model_func`.
#'     \item `lower_bounds`: A named list that specifies the lower bounds for parameters in `model_func`.
#'     \item `upper_bounds`:A named list that specifies the upper bounds for parameters in `model_func`.
#'   }
#' @param lower_bounds TODO: description.
#' @description
#' Calculates and visualises the Longitudinal Growth Rate (LGR) score
#' from time-course drug response data. Fits a dose-response curve to
#' normalised area-over-curve values and reports both raw and
#' goodness-of-fit-adjusted LGR scores.
#' @return
#' A \code{ggplot2} object displaying:
#'   \itemize{
#'     \item Points showing normalised AOC values at each concentration
#'     \item Fitted dose-response curve
#'     \item Title with LGR score and adjusted LGR score (LGR × R²)
#'   }
#'   The x-axis is formatted for log-molar concentrations.
#' @details
#' Higher LGR scores indicate greater overall drug efficacy. The adjusted
#' score penalises poor model fits.
#' @seealso \code{\link{calc_inhibition_metrics}}, \code{\link{fit_curve}},
#'   \code{\link{interpolate_curve_concentration}}, \code{\link{theme_dosefitr}}
#' @examples
#' \dontrun{
#' # Calculate and plot LGR score
#' lgr_plot <- plot_lgr_score(
#'   model_list = fitted_models$model_list,
#'   x_var = "time",
#'   metric = "gr",
#'   treatment_column = "treatment",
#'   negative_control_name = "DMSO",
#'   positive_control_name = NA,
#'   concentration_column = "dose_uM"
#' )
#' # With custom dose-response model
#' lgr_plot <- plot_lgr_score(
#'   model_list = fitted_models$model_list,
#'   x_var = "time",
#'   metric = "gr",
#'   treatment_column = "treatment",
#'   negative_control_name = "DMSO",
#'   positive_control_name = "Staurosporine",
#'   concentration_column = "dose_uM",
#'   curve_opts = list(model = "four_param_sigmoid_log")
#' )
#' }
#' @importFrom ggplot2 ggplot geom_point geom_line aes scale_x_continuous labs guides ggtitle
#' @export
plot_lgr_score <- function(
  model_list,
  x_var,
  metric = "gr",
  treatment_column,
  negative_control_name,
  positive_control_name = NA,
  concentration_column,
  curve_opts = list(
    model = "five_param_sigmoid_log",
    lower_bounds = list(bottom = 0),
    upper_bounds = list(top = 100)
  )
) {
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
    curve_opts = curve_opts
  )

  predicted <- predict_data(
    obj = fitted_model$fit$obj,
    lower_x = min(fitted_model$data$x),
    upper_x = max(fitted_model$data$x)
  )

  auc <- auc_trapezoid(predicted$x, predicted$y)

  lgr_score <- auc / (max(auc_data$concentration) - min(auc_data$concentration))

  # GOF - this is what dprl uses
  ss_res <- sum(residuals(fitted_model$fit$obj)^2)
  ss_total <- sum((auc_data$aoc - mean(auc_data$aoc))^2)
  gof <- 1 - (ss_res / ss_total)

  lgr_score_adjusted <- lgr_score * gof

  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(
      auc_data,
      mapping = ggplot2::aes(
        concentration,
        aoc,
      ),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE,
    ) +
    ggplot2::geom_line(
      predicted,
      mapping = ggplot2::aes(
        x,
        y
      ),
      key_glyph = ggplot2::draw_key_rect,
      show.legend = TRUE
    ) +
    ggplot2::scale_x_continuous(
      labels = \(x) log_m_to_str(x),
      minor_breaks = log10(rep(1:9, 21) * (10^rep(-10:10, each = 9)))
    ) +
    ggplot2::labs(x = concentration_column, y = "Normalised AOC") +
    ggplot2::guides(colour = NULL) +
    ggplot2::ggtitle(
      paste(
        paste("LGR:", format(lgr_score, digits = 4), sep = " "),
        "\n",
        paste("LGR Adj.:", format(lgr_score_adjusted, digits = 4), sep = " "),
        sep = ""
      )
    ) +
    theme_dosefitr()

  return(gg)
}

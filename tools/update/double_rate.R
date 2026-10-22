#' Calculate Doubling Time from Image-Based Growth Data
#' @description
#' Calculate cell doubling time from data by identifying
#' the exponential growth phase and fitting a linear model to log-transformed data.
#' @param data A data frame containing time-series growth measurements.
#' @param cell_column Optional character string specifying a column for grouping multiple known cell populations.
#' @param dimension_factor Numeric multiplier to account for dimensional readout (default = 1 for linear measurements, use 1.5 for area measurements).
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @return
#' A list containing:
#'  \itemize{
#'    \item `double_time`: Numeric value of calculated doubling time.
#'    \item `plot`: \code{ggplot} object showing linear and log-scale growth curves with fitted lines.
#'  }
#' @importFrom geomtextpath geom_textabline
#' @importFrom inflection ese
#' @importFrom stats coef lm
#' @export
calc_double_rate_data <- function(
  data,
  x_var,
  y_var,
  cell_column = NULL,
  dimension_factor = 1
) {
  data <- prep_data(data, x_var, y_var)

  if (!is.null(cell_column)) {
    data$cells <- data[[cell_column]]

    # Checks
    if (length(unique(data$cells)) == 1) {
      stop(sprintf("Only 1 unique value found for %s", cell_column))
    }
    if (!is.numeric(data$cells)) {
      stop(sprintf("Provided cell_column `%s` is not numeric"))
    }

    is_multiple_per_x <- length(unique(data$y)) > length(unique(data$x))

    if (is_multiple_per_x) {
      # Group by and summarise
      data_mean <- do.call(
        rbind,
        lapply(split(data, list(data$x, data$cells)), function(df) {
          data.frame(
            x = df$x[1],
            cells = df$cells[1],
            y = mean(df$y)
          )
        })
      )

      # Filter by cell
      data_filt_list <- lapply(unique(data_mean$cells), function(a) {
        df <- data_mean[data_mean$cells == a, ]

        # Find where section is linear, this is exponential growth
        inflection_points <- inflection::ese(df$x, log(df$y), 0)
        linear_indexes <- inflection_points[1]:inflection_points[2]

        df[linear_indexes, ]
      })
      data_filt <- do.call(rbind, data_filt_list)

      # Add n0 and log_n0 by group
      data_filt <- do.call(
        rbind,
        lapply(split(data_filt, data_filt$cells), function(df) {
          df$n0 <- df$y[1]
          df$log_n0 <- log(df$n0)
          df
        })
      )

      # Left join
      a <- merge(
        data_filt[, c("x", "cells", "n0", "log_n0")],
        data[, c("x", "y", "cells")],
        by = c("x", "cells"),
        all.x = TRUE
      )

      model <- stats::lm(
        log(y) ~ x + offset(log_n0),
        a
      )
      model_double <- stats::lm(
        dimension_factor * log(y) ~ x + offset(log_n0),
        data_filt
      )
    } else {
      # Filter by cell
      data_filt_list <- lapply(unique(data$cells), function(a) {
        df <- data[data$cells == a, ]

        # Find where section is linear, this is exponential growth
        inflection_points <- inflection::ese(df$x, log(df$y), 0)
        linear_indexes <- inflection_points[1]:inflection_points[2]

        df[linear_indexes, ]
      })
      data_filt <- do.call(rbind, data_filt_list)

      # Add n0 and log_n0 by group
      data_filt <- do.call(
        rbind,
        lapply(split(data_filt, data_filt$cells), function(df) {
          df$n0 <- df$y[1]
          df$log_n0 <- log(df$n0)
          df
        })
      )

      model <- stats::lm(
        log(y) ~ x + offset(log_n0),
        data_filt
      )
      model_double <- stats::lm(
        dimension_factor * log(y) ~ x + offset(log_n0),
        data_filt
      )
    }
  } else {
    # Single sample
    is_multiple_per_x <- length(unique(data$y)) > length(unique(data$x))

    if (is_multiple_per_x) {
      # Group by and summarise
      data_mean <- do.call(
        rbind,
        lapply(split(data, data$x), function(df) {
          data.frame(x = df$x[1], y = mean(df$y))
        })
      )

      # Find where section is linear, this is exponential growth
      inflection_points <- inflection::ese(data_mean$x, log(data_mean$y), 0)
      linear_indexes <- inflection_points[1]:inflection_points[2]

      data_filt <- data_mean[linear_indexes, ]

      # Left join
      a <- merge(
        data_filt[, "x", drop = FALSE],
        data[, c("x", "y")],
        by = "x",
        all.x = TRUE
      )

      # Fit model on single
      model <- stats::lm(log(y) ~ x, data_filt)
      model_double <- stats::lm(
        dimension_factor * log(y) ~ x,
        a
      )
    } else {
      # Find where section is linear, this is exponential growth
      inflection_points <- inflection::ese(data$x, log(data$y), 0)
      linear_indexes <- inflection_points[1]:inflection_points[2]

      data_filt <- data[linear_indexes, ]

      # Fit model on single
      model <- lm(log(y) ~ x, data_filt)
      model_double <- lm(
        dimension_factor * log(y) ~ x,
        data_filt
      )
    }
  }

  # Double calculation
  double_time <- log(2) / stats::coef(model_double)[["x"]]

  ## Plots
  # Add in linear section
  if (!is.null(cell_column)) {
    data_list <- lapply(unique(data_filt$cells), function(a) {
      data_filt_cell <- data_filt[data_filt$cells == a, ]
      data_cell <- data[data$cells == a, ]
      data_cell$linear <- data_cell$x %in% data_filt_cell$x
      data_cell
    })
    data <- do.call(rbind, data_list)
  } else {
    data$linear <- data$x %in% data_filt$x
  }

  # Colours
  if (!is.null(cell_column)) {
    n_colors <- length(unique(data$cells))
    hues <- seq(15, 375, length.out = n_colors + 1)
    colours <- grDevices::hcl(h = hues, c = 100, l = 65)[1:n_colors]
    names(colours) <- unique(data$cells)
  }

  # Create aes()
  if (!is.null(cell_column)) {
    aes_lin <- ggplot2::aes(x, y, colour = factor(cells))
    aes_log <- ggplot2::aes(x, log(y), colour = factor(cells))
  } else {
    aes_lin <- ggplot2::aes(x, y)
    aes_log <- ggplot2::aes(x, log(y))
  }

  # Linear
  if (is_multiple_per_x) {
    gg_lin <- ggplot2::ggplot(data, aes_lin) +
      ggplot2::geom_point() +
      ggplot2::guides(colour = guide_legend(title = cell_column))
  } else {
    gg_lin <- ggplot2::ggplot(data, aes_lin) +
      ggplot2::geom_line() +
      ggplot2::guides(colour = guide_legend(title = cell_column))
  }

  gg_lin <- gg_lin +
    ggplot2::labs(x = x_var, y = y_var) +
    theme_dosefitr()

  # Add vertical lines
  if (!is.null(cell_column)) {
    for (a in unique(data$cells)) {
      data_filt_cell <- data_filt[data_filt$cells == a, ]
      gg_lin <- gg_lin +
        ggplot2::geom_vline(
          xintercept = min(data_filt_cell$x),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        ) +
        ggplot2::geom_vline(
          xintercept = max(data_filt_cell$x),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        )
    }
  } else {
    gg_lin <- gg_lin +
      ggplot2::geom_vline(
        xintercept = min(data_filt$x),
        linetype = "dashed"
      ) +
      ggplot2::geom_vline(
        xintercept = max(data_filt$x),
        linetype = "dashed"
      )
  }

  # Logarthmic
  if (is_multiple_per_x) {
    gg_log <- ggplot2::ggplot(data, aes_log) +
      ggplot2::geom_point() +
      ggplot2::guides(colour = ggplot2::guide_legend(title = cell_column))
  } else {
    gg_log <- ggplot2::ggplot(data, aes_log) +
      ggplot2::geom_line() +
      ggplot2::guides(colour = ggplot2::guide_legend(title = cell_column))
  }

  # Intercept(s)
  if (!is.null(cell_column)) {
    intercept <- sapply(unique(data$cells), function(a) {
      data_filt_cell <- data_filt[data_filt$cells == a, ]
      intercept <- data_filt_cell$log_n0[1] + stats::coef(model)[1]
      intercept
    })
  } else {
    intercept <- stats::coef(model)[[1]]
  }

  gg_log <- gg_log +
    ggplot2::labs(x = x_var, y = paste0("log(", y_var, ")")) +
    theme_dosefitr()

  # Add vertical lines
  if (!is.null(cell_column)) {
    for (a in unique(data$cells)) {
      data_filt_cell <- data_filt[data_filt$cells == a, ]
      gg_log <- gg_log +
        ggplot2::geom_vline(
          xintercept = min(data_filt_cell$x),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        ) +
        ggplot2::geom_vline(
          xintercept = max(data_filt_cell$x),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        )
    }
  } else {
    gg_log <- gg_log +
      ggplot2::geom_vline(
        xintercept = min(data_filt$x),
        linetype = "dashed"
      ) +
      ggplot2::geom_vline(
        xintercept = max(data_filt$x),
        linetype = "dashed"
      )
  }

  # Add regression lines
  for (i in seq_along(intercept)) {
    slope <- stats::coef(model)[[2]]
    label_text <- paste0(
      "y = ",
      signif(slope, digits = 3),
      "\u00d7x + ",
      signif(intercept[i], digits = 3)
    )

    gg_log <- gg_log +
      geomtextpath::geom_textabline(
        intercept = intercept[i],
        slope = stats::coef(model)[[2]],
        label = label_text,
        size = 6 / ggplot2::.pt,
        vjust = 1.5
      )
  }

  return(list(double_time = double_time, linear = gg_lin, log = gg_log))
}

#' Calculate Doubling Time from Passage Data
#' @description
#' Calculate cell doubling time from known cell counts at two different time points,
#' typically from passage data or manual cell counting.
#' @param start_date Character string or Date object representing the first measurement date in format `YYYY-MM-DD`.
#' @param end_date Character string or Date object representing the second measurement date in format `YYYY-MM-DD`.
#' @param start_n Numeric value for the number of cells at the first measurement.
#' @param end_n Numeric value for the number of cells at the second measurement.
#' @return
#' Numeric value representing the doubling time in hours.
#' @details
#' Uses the exponential growth formula to calculate doubling time:
#' t_double = (t2 - t1) * ln(2) / ln(N2/N1)
#' @export
calc_double_rate_passage <- function(start_date, end_date, start_n, end_n) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  days_between <- as.numeric(difftime(end_date, start_date))

  return(days_between * log(2) / (log(end_n / start_n)) * 24)
}

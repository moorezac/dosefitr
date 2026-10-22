#' Calculate Doubling Time from Image-Based Growth Data
#'
#' @description
#' Calculate cell doubling time from data by identifying the exponential
#' growth phase and fitting a linear model to log-transformed data.
#'
#' @param data A data frame containing time-series growth measurements.
#' @param cell_column Optional character string specifying a column for
#'   grouping multiple known cell populations.
#' @param doubling_readout_factor Numeric value greater than 1 describing
#'   the fold-change in the measured readout corresponding to one biological
#'   doubling. For example, use 2 for cell number, where a two-fold increase
#'   represents one doubling, or 1.5 for spheroid cross-sectional area, where
#'   a 1.5-fold increase is interpreted as one doubling.
#' @param x_var Character string specifying the column name for the
#'   independent variable (typically time).
#' @param y_var Character string specifying the column name for the
#'   dependent variable.
#'
#' @return
#' A list containing:
#'   \itemize{
#'     \item `double_time`: Numeric value of calculated doubling time, in
#'       the units of `x_var`.
#'     \item `model`: Fitted linear model.
#'     \item `window`: Data used to identify the exponential window.
#'     \item `linear`: ggplot object showing the linear-scale growth curve.
#'     \item `log`: ggplot object showing the log-scale growth curve.
#'   }
#'
#' @export
calc_double_rate_data <- function(
  data,
  x_var,
  y_var,
  cell_column = NULL,
  doubling_readout_factor = 2
) {
  data <- prepare_growth_data(
    data = data,
    x_var = x_var,
    y_var = y_var,
    cell_column = cell_column
  )

  window_data <- find_exponential_windows(data)

  fit_data <- prepare_doubling_fit_data(
    data = data,
    window_data = window_data
  )

  fit <- fit_doubling_model(
    data = fit_data,
    doubling_readout_factor = doubling_readout_factor
  )

  plots <- plot_doubling_data(
    data = data,
    window_data = window_data,
    model = fit$model,
    x_lab = x_var,
    y_lab = y_var,
    legend_title = cell_column
  )

  list(
    double_time = fit$double_time,
    model = fit$model,
    window = window_data,
    linear = plots$linear,
    log = plots$log
  )
}


#' Prepare growth data
#'
#' @param data A data frame.
#' @param x_var,y_var Column names of time and the growth readout.
#' @param cell_column Optional grouping column, or `NULL`.
#' @return Prepared data with a `cells` column (`"all"` when ungrouped).
#' @keywords internal
prepare_growth_data <- function(
  data,
  x_var,
  y_var,
  cell_column = NULL
) {
  data <- prepare_data(data, x_var, y_var)

  if (is.null(cell_column)) {
    data$cells <- "all"
  } else {
    if (!cell_column %in% colnames(data)) {
      stop(sprintf("Column not found in data: %s", cell_column))
    }

    data$cells <- data[[cell_column]]

    if (length(unique(data$cells)) == 1L) {
      stop(sprintf(
        "Only 1 unique value found for '%s'.",
        cell_column
      ))
    }
  }

  data
}


#' Identify exponential growth windows
#'
#' Replicates are collapsed only to identify the exponential window.
#' The full replicate-level data are retained for model fitting.
#'
#' @param data Data from [prepare_growth_data()].
#' @return The rows inside each group's exponential window, with `n0` and
#'   `log_n0` columns.
#' @keywords internal
find_exponential_windows <- function(data) {
  has_replicates <- anyDuplicated(
    data[, c("x", "cells")]
  ) >
    0L

  window_source <- if (has_replicates) {
    collapse_replicates(data)
  } else {
    data
  }

  window_data <- lapply(
    split(window_source, window_source$cells),
    function(df) {
      add_n0(exponential_window(df))
    }
  )

  window_data <- do.call(rbind, window_data)

  rownames(window_data) <- NULL

  window_data
}


#' Prepare data for doubling-time model
#'
#' @param data Data from [prepare_growth_data()].
#' @param window_data Windows from [find_exponential_windows()].
#' @return Replicate-level observations inside each window, with `n0` and
#'   `log_n0`.
#' @keywords internal
prepare_doubling_fit_data <- function(data, window_data) {
  merge(
    window_data[, c("x", "cells", "n0", "log_n0")],
    data[, c("x", "y", "cells")],
    by = c("x", "cells"),
    all.x = TRUE
  )
}


#' Fit doubling-time model
#'
#' @param data Data from [prepare_doubling_fit_data()].
#' @param doubling_readout_factor Fold-change in the readout per doubling.
#' @return A list with the fitted `model` and `double_time`.
#' @keywords internal
fit_doubling_model <- function(
  data,
  doubling_readout_factor = 2
) {
  if (
    !is.numeric(doubling_readout_factor) ||
      length(doubling_readout_factor) != 1L ||
      is.na(doubling_readout_factor) ||
      doubling_readout_factor <= 1
  ) {
    stop(
      "`doubling_readout_factor` must be a single numeric value greater than 1."
    )
  }

  model <- stats::lm(
    log(y) ~ x + offset(log_n0),
    data = data
  )

  slope <- stats::coef(model)[["x"]]

  if (!is.finite(slope) || slope <= 0) {
    warning(
      "No growth detected in the exponential window (slope <= 0); the doubling time is not meaningful.",
      call. = FALSE
    )
  }

  double_time <- log(doubling_readout_factor) / slope

  list(
    model = model,
    double_time = double_time
  )
}


#' Plot doubling-time data
#'
#' @param data Data from [prepare_growth_data()].
#' @param window_data Windows from [find_exponential_windows()].
#' @param model The fitted model from [fit_doubling_model()].
#' @param x_lab,y_lab Axis labels.
#' @param legend_title Legend title for grouped data.
#' @return A list with the `linear` and `log` plots.
#' @keywords internal
plot_doubling_data <- function(
  data,
  window_data,
  model,
  x_lab,
  y_lab,
  legend_title = NULL
) {
  has_replicates <- anyDuplicated(
    data[, c("x", "cells")]
  ) >
    0L

  grouped <- length(unique(data$cells)) > 1L

  colours <- growth_colours(
    cells = data$cells,
    grouped = grouped
  )

  linear_plot <- growth_panel(
    data = data,
    y_expr = quote(y),
    y_lab = y_lab,
    x_lab = x_lab,
    points = has_replicates,
    colours = colours,
    legend_title = legend_title
  ) +
    window_vlines(window_data, colours)

  log_plot <- growth_panel(
    data = data,
    y_expr = quote(log(y)),
    y_lab = paste0("log(", y_lab, ")"),
    x_lab = x_lab,
    points = has_replicates,
    colours = colours,
    legend_title = legend_title
  ) +
    window_vlines(window_data, colours)

  log_plot <- add_growth_fit_lines(
    plot = log_plot,
    model = model,
    window_data = window_data
  )

  list(
    linear = linear_plot,
    log = log_plot
  )
}

#' Generate colours for growth groups
#'
#' @param cells The `cells` column.
#' @param grouped Whether the data are grouped.
#' @return A named character vector of colours, or `NULL` if ungrouped.
#' @keywords internal
growth_colours <- function(cells, grouped = TRUE) {
  if (!grouped) {
    return(NULL)
  }

  cell_levels <- sort(unique(cells))

  hues <- seq(
    15,
    375,
    length.out = length(cell_levels) + 1
  )

  colours <- grDevices::hcl(
    h = hues,
    c = 100,
    l = 65
  )[seq_along(cell_levels)]

  names(colours) <- cell_levels

  colours
}


#' Add fitted growth lines to log plot
#'
#' @param plot The log-scale plot.
#' @param model The fitted model from [fit_doubling_model()].
#' @param window_data Windows from [find_exponential_windows()].
#' @return `plot` with one labelled fitted line per group.
#' @keywords internal
add_growth_fit_lines <- function(
  plot,
  model,
  window_data
) {
  cell_levels <- sort(unique(window_data$cells))

  slope <- stats::coef(model)[["x"]]

  intercepts <- vapply(
    cell_levels,
    function(cell) {
      log_n0 <- window_data$log_n0[
        window_data$cells == cell
      ][1]

      log_n0 + stats::coef(model)[[1]]
    },
    numeric(1)
  )

  for (intercept in intercepts) {
    plot <- plot +
      geomtextpath::geom_textabline(
        intercept = intercept,
        slope = slope,
        label = sprintf(
          "y = %s\u00d7x + %s",
          signif(slope, 3),
          signif(intercept, 3)
        ),
        size = 8 / ggplot2::.pt,
        vjust = 1.5,
        colour = "black"
      )
  }

  plot
}


#' Layers marking the exponential growth window
#'
#' @description
#' Dashed vertical lines at the first and last timepoint of each group's
#' exponential window, coloured to match the group.
#'
#' @param window_data The windowed data, with `x` and `cells` columns.
#' @param colours Named character vector of group colours, or `NULL` for
#'   a single ungrouped black line.
#' @return A list of ggplot2 layers.
#' @keywords internal
window_vlines <- function(window_data, colours = NULL) {
  bounds <- stats::aggregate(
    x ~ cells,
    data = window_data,
    FUN = function(v) c(min = min(v), max = max(v))
  )

  bounds <- data.frame(
    cells = bounds$cells,
    min = bounds$x[, "min"],
    max = bounds$x[, "max"]
  )

  unlist(
    lapply(seq_len(nrow(bounds)), function(i) {
      colour <- if (is.null(colours)) {
        "black"
      } else {
        unname(colours[as.character(bounds$cells[i])])
      }

      list(
        ggplot2::geom_vline(
          xintercept = bounds$min[i],
          colour = colour,
          linetype = "dashed"
        ),
        ggplot2::geom_vline(
          xintercept = bounds$max[i],
          colour = colour,
          linetype = "dashed"
        )
      )
    }),
    recursive = FALSE
  )
}


#' Identify exponential growth window
#'
#' @param df One group's data, with `x` and positive `y`.
#' @return The rows between the extreme points found by
#'   \code{inflection::ese()} on `log(y)`.
#' @keywords internal
exponential_window <- function(df) {
  df <- df[order(df$x), , drop = FALSE]

  if (any(df$y <= 0)) {
    stop(
      "Growth values must be positive to identify an exponential window."
    )
  }

  points <- inflection::ese(
    df$x,
    log(df$y),
    0
  )

  df[points[1]:points[2], , drop = FALSE]
}


#' Add initial value to window data
#'
#' @param df One group's windowed data.
#' @return `df` with `n0` (the first `y` in the window) and `log_n0`.
#' @keywords internal
add_n0 <- function(df) {
  df$n0 <- df$y[1]
  df$log_n0 <- log(df$n0)

  df
}


#' Collapse replicates
#'
#' @param data Data with `y` and the `by` columns.
#' @param by Columns identifying replicates.
#' @return The mean `y` for each combination of `by`.
#' @keywords internal
collapse_replicates <- function(
  data,
  by = c("x", "cells")
) {
  out <- stats::aggregate(
    data[, "y", drop = FALSE],
    by = data[, by, drop = FALSE],
    FUN = mean,
    na.rm = TRUE
  )

  out[order(out$cells, out$x), , drop = FALSE]
}


#' Plot growth data
#'
#' @param data Data from [prepare_growth_data()].
#' @param y_expr Quoted expression for the y aesthetic, e.g. `quote(log(y))`.
#' @param y_lab,x_lab Axis labels.
#' @param points Draw points (for replicate data) instead of lines.
#' @param colours Named group colours, or `NULL` if ungrouped.
#' @param legend_title Legend title for grouped data.
#' @return A `ggplot` object.
#' @keywords internal
growth_panel <- function(
  data,
  y_expr,
  y_lab,
  x_lab,
  points,
  colours = NULL,
  legend_title = NULL
) {
  mapping <- if (is.null(colours)) {
    ggplot2::aes(x = x, y = !!y_expr)
  } else {
    ggplot2::aes(x = x, y = !!y_expr, colour = factor(cells))
  }

  gg <- ggplot2::ggplot(data, mapping) +
    (if (points) {
      ggplot2::geom_point()
    } else {
      ggplot2::geom_line()
    }) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_dosefitr()

  if (!is.null(colours)) {
    gg <- gg +
      ggplot2::scale_colour_manual(
        values = colours
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          title = legend_title
        )
      )
  }

  gg
}


#' Calculate Doubling Time from Passage Data
#'
#' @description
#' Calculate cell doubling time from known cell counts at two different
#' time points, typically from passage data or manual cell counting.
#'
#' @param start_date Character string or Date object representing the first
#'   measurement date in format `YYYY-MM-DD`.
#' @param end_date Character string or Date object representing the second
#'   measurement date in format `YYYY-MM-DD`.
#' @param start_n Numeric value for the number of cells at the first
#'   measurement.
#' @param end_n Numeric value for the number of cells at the second
#'   measurement.
#'
#' @return
#' Numeric value representing the doubling time in hours.
#'
#' @details
#' Uses the exponential growth formula to calculate doubling time:
#' double_time = (end_date - start_date) * ln(2) / ln(end_n/start_n)
#'
#' @export
calc_double_rate_passage <- function(
  start_date,
  end_date,
  start_n,
  end_n
) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  days_between <- as.numeric(
    difftime(
      end_date,
      start_date,
      units = "days"
    )
  )

  days_between * log(2) / log(end_n / start_n) * 24
}

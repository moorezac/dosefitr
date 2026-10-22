#' Calculate Area Under Curve for Fitted Model via Integration
#' @description
#' Calculate the area under the curve (AUC) for a fitted model between specified bounds.
#' @param lower_x The lower bound for integration.
#' @param upper_x The upper bound for integration.
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#' @return
#' A list of class `integrate` containing:
#'  \itemize{
#'    \item{value}{The calculated area under the curve}
#'    \item{abs.error}{Estimate of absolute error}
#'    \item{subdivisions}{Number of subdivisions used}
#'    \item{message}{Status message}
#'  }
#' @examples
#' \dontrun{
#' # Calculate total growth over time
#' auc_result <- calculate_auc(
#'   obj = growth_fit,
#'   lower_x = 0,
#'   upper_x = 72  # 72 hours
#' )
#' total_growth <- auc_result$value
#' }
#' @export
calculate_auc <- function(obj, lower_x, upper_x) {
  integrate(
    f = function(x) {
      predict(obj, data.frame(x = x))
    },
    lower = lower_x,
    upper = upper_x
  )
}


#' Generate Predictions from Fitted Model Object
#' @description
#' Generate smooth prediction curves from fitted models over a specified range
#'    of x values. Supports both single and grouped models.
#' @param lower_x The minimum x value for predictions.
#' @param upper_x The maximum x value for predictions.
#' @param group Optional character vector of group values for shared parameter models.
#' @param num_points Number of prediction points between lower_x and upper_x (default: 1000).
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#' @return
#' A tibble containing:
#'  \itemize{
#'    \item{x}{Prediction x values}
#'    \item{y}{Predicted y values}
#'    \item{group}{Group identifier (if applicable)}
#'  }
#' @examples
#' \dontrun{
#' # Generate smooth prediction curve
#' predictions <- predict_data(
#'   obj = growth_fit,
#'   lower_x = 0,
#'   upper_x = 100,
#'   num_points = 500
#' )
#' # For grouped model
#' predictions <- predict_data(
#'   obj = grouped_fit,
#'   lower_x = 0,
#'   upper_x = 100,
#'   group = c("cell_line_A", "cell_line_B")
#' )
#' }
#' @export
predict_data <- function(
  obj,
  lower_x,
  upper_x,
  group = NA,
  num_points = 1e3
) {
  x_values <- seq(lower_x, upper_x, length.out = num_points)

  if (!any(is.na(group))) {
    # For each group, predict and store results
    res_list <- lapply(group, function(i) {
      y <- predict(obj, newdata = data.frame(x = x_values, group = i))
      data.frame(x = x_values, y = y, group = i, stringsAsFactors = FALSE)
    })
    # Combine all into one data.frame
    do.call(rbind, res_list)
  } else {
    y <- predict(obj, newdata = data.frame(x = x_values))
    data.frame(x = x_values, y = y, stringsAsFactors = FALSE)
  }
}


#' Find X Value for Target Y
#' Uses root-finding to determine the x value at which the fitted model
#' predicts a specified target y value.
#' @param x_values Numeric vector of x values defining the search interval.
#' @param target Numeric value specifying the target y value.
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#' @return
#' A list as returned by \code{\link[stats]{uniroot}}, containing:
#'  \itemize{
#'    \item `root` T:he x value where predicted y equals target.
#'    \item `f.root` :The value of the function at the root (should be near 0).
#'    \item `iter` :Number of iterations.
#'    \item `estim.prec` :Estimated precision.
#'  }
#' @seealso \code{\link{fit_curve}}, \code{\link[stats]{uniroot}}
#' @examples
#' \dontrun{
#' # Find IC50 (concentration at which response = 0.5)
#' ic50 <- find_x_for_y(
#'   obj = dose_response_fit,
#'   x_values = predicted$x,
#'   target = 0.5
#' )
#' }
#' @export
find_x_for_y <- function(obj, x_values, target) {
  uniroot(
    f = function(x) {
      predict(obj, data.frame(x = x)) - target
    },
    interval = c(min(x_values), max(x_values))
  )
}


#' Calculate Area Under Curve Using Trapezoidal Rule
#' @param x Numeric vector of x values (must be sorted).
#' @param y Numeric vector of y values.
#' @return
#' Numeric scalar representing the area.
#' @keywords internal
auc_trapezoid <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have equal lengths")
  }
  if (length(x) < 2) {
    return(0)
  }
  sum((y[-length(y)] + y[-1]) / 2 * diff(x))
}

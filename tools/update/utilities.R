#' Prepare Data
#' @param data A dataframe.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @description
#' Create `x` and `y` columns from specified variable names.
#' @keywords internal
prep_data <- function(data, x_var, y_var) {
  # Extract columns by name
  x <- data[[x_var]]
  y <- data[[y_var]]

  # If already done
  if (all(c("x", "y", "x_original", "y_original") %in% colnames(data))) {
    return(data)
  }

  # Handle existing columns
  cols_to_rename <- c("x", "y")[c("x", "y") %in% colnames(data)]
  if (length(cols_to_rename) > 0) {
    names(data)[names(data) %in% cols_to_rename] <- paste(
      cols_to_rename,
      "_original",
      sep = ""
    )
  }

  # Subset data and coerce to numeric
  keep <- !(is.na(x) | is.na(y))

  data <- data[keep, , drop = FALSE]
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  # Reorder columns to have x, y first
  result <- cbind(data.frame(x = x, y = y), data)
  rownames(result) <- NULL

  return(result)
}


#' Convert Log Molar Concentrations to Human-Readable Strings
#' @description
#' Convert log10 molar concentration values to strings with appropriate
#'    unit prefixes (M, mM, µM, nM, pM, fM) for display in plots and tables.
#' @param log_m Numeric vector of log10 molar concentration values.
#' @param digits Number of decimal places for formatting (default: 2).
#' @return
#' Character vector of formatted concentration strings.
#' @examples
#' # Convert log concentrations to readable format
#' log_concs <- c(-9, -6, -3, 0)
#' formatted <- log_m_to_str(log_concs)
#' # Returns: c("1.00 nM", "1.00 µM", "1.00 mM", "1.00 M")
#' @export
log_m_to_str <- function(log_m, digits = 2) {
  sapply(
    log_m,
    function(x) {
      if (is.na(x)) {
        return(NA_character_)
      }
      if (x == 0) {
        return("Negative Control")
      }
      units <- c("M", "mM", "µM", "nM", "pM", "fM")
      scale <- c(0, -3, -6, -9, -12, -15)
      idx <- min(which(x >= scale))
      value <- 10^x / 10^scale[idx]
      sprintf(paste0("%.", digits, "f %s"), value, units[idx])
    },
    USE.NAMES = FALSE
  )
}

stop_unable_to_fit <- function(
  message,
  errors = character(),
  call = sys.call(-1)
) {
  stop(structure(
    class = c("dosefitr_fit_error", "error", "condition"),
    list(
      message = message,
      call = call,
      errors = errors
    )
  ))
}

check_fit_data <- function(data, x_var, y_var, n_params) {
  missing_cols <- setdiff(c(x_var, y_var), names(data))
  if (length(missing_cols)) {
    stop_unable_to_fit(sprintf(
      "Column%s not found in data: %s",
      if (length(missing_cols) > 1) "s" else "",
      paste(missing_cols, collapse = ", ")
    ))
  }

  x <- data[[x_var]]
  y <- data[[y_var]]
  keep <- is.finite(x) & is.finite(y)

  if (!any(keep)) {
    stop_unable_to_fit(
      "No finite observations: every row has a missing or infinite value."
    )
  }
  if (sum(keep) < n_params) {
    stop_unable_to_fit(sprintf(
      "Not enough data to fit: %d usable observation%s for a model with %d parameters.",
      sum(keep),
      if (sum(keep) == 1) "" else "s",
      n_params
    ))
  }
  if (length(unique(y[keep])) == 1) {
    stop_unable_to_fit(sprintf(
      "'%s' is constant (%g) across all observations, so no curve can be fitted.",
      y_var,
      y[keep][1]
    ))
  }
  invisible(TRUE)
}

#' Translate a fitting error message into an actionable hint
#'
#' @param msg Character vector of error messages from a failed fit.
#' @return Character vector of hints, one per matched message. Messages with
#'   no matching pattern are dropped, so the result may be shorter than `msg`
#'   or empty.
#' @keywords internal
translate_fit_error <- function(msg) {
  patterns <- list(
    c(
      "singular gradient",
      "The model could not be distinguished from a simpler one with this data (singular gradient). Try a model with fewer parameters, or supply starting values."
    ),
    c(
      "NA/NaN/Inf|Missing value or an infinity",
      "The model produced non-finite values during fitting. Check for zero or negative values where the model expects positive ones (e.g. log or power terms)."
    ),
    c(
      "number of iterations exceeded|maximum number of iterations",
      "Fitting did not converge within the iteration limit. Try increasing 'maxiter', or supply starting values closer to the expected result."
    ),
    c(
      "parameters without starting value|'start'",
      "One or more model parameters had no starting value."
    ),
    c(
      "could not find function|object '.*' not found",
      "The model formula referred to something that doesn't exist. This is likely a problem with the model definition rather than your data."
    )
  )

  hints <- vapply(
    msg,
    function(m) {
      for (p in patterns) {
        if (grepl(p[1], m, ignore.case = TRUE)) return(p[2])
      }
      NA_character_
    },
    character(1),
    USE.NAMES = FALSE
  )

  unique(hints[!is.na(hints)])
}

#' Build a diagnostic message for a failed fit
#'
#' @param errors Character vector of underlying error messages.
#' @param prefix Character; the leading sentence, without trailing punctuation.
#' @return A single character string: `prefix`, followed by any actionable
#'   hints on their own lines, or a pointer to the condition's `errors`
#'   element when no hint applies.
#' @keywords internal
# no_fit_message <- function(
#   errors,
#   prefix = "Unable to fit a model to this data"
# ) {
#   hints <- translate_fit_error(errors)

#   if (!length(hints)) {
#     return(paste0(
#       prefix,
#       ". Use conditionCall() or the 'errors' element of this condition for the underlying message."
#     ))
#   }

#   paste0(prefix, ":\n", paste(hints, collapse = "\n"))
# }

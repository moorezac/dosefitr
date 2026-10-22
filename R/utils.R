#' Standardise data to `x` and `y` columns
#'
#' Copies `x_var` and `y_var` into numeric `x` and `y` columns at the front of
#' the data, dropping rows where either is missing or non-finite. Any existing
#' `x`/`y` columns are renamed to `x_original`/`y_original`. Data that already
#' has all four columns is returned unchanged.
#'
#' @param data A data frame.
#' @param x_var,y_var Column names of the independent and dependent variables.
#' @return `data` with `x` and `y` prepended and row names reset.
#' @keywords internal
prepare_data <- function(data, x_var, y_var) {
    if (all(c("x", "y", "x_original", "y_original") %in% colnames(data))) {
        return(data)
    }

    x <- as.numeric(data[[x_var]])
    y <- as.numeric(data[[y_var]])

    cols_to_rename <- c("x", "y")[c("x", "y") %in% colnames(data)]
    if (length(cols_to_rename) > 0) {
        names(data)[names(data) %in% cols_to_rename] <- paste0(
            cols_to_rename,
            "_original"
        )
    }

    # Matches check_fit_data(): Inf/-Inf rows are dropped along with NA/NaN
    keep <- is.finite(x) & is.finite(y)

    data <- data[keep, , drop = FALSE]
    x <- x[keep]
    y <- y[keep]

    result <- cbind(data.frame(x = x, y = y), data)
    rownames(result) <- NULL

    result
}

#' Format log10 molar concentrations as readable strings
#'
#' @param log_m Numeric vector of log10 molar concentrations.
#' @param digits Number of decimal places.
#' @return Character vector such as `"1.00 µM"`. `0` is labelled
#'   `"Negative Control"` and `NA` stays `NA`.
#' @keywords internal
log_m_to_str <- function(log_m, digits = 2) {
    vapply(
        log_m,
        function(x) {
            if (is.na(x)) {
                return(NA_character_)
            }
            if (x == 0) {
                return("Negative Control")
            }
            units <- c("M", "mM", "\u03bcM", "nM", "pM", "fM")
            scale <- c(0, -3, -6, -9, -12, -15)
            idx <- min(which(x >= scale), length(scale))
            value <- 10^x / 10^scale[idx]
            sprintf(paste0("%.", digits, "f %s"), value, units[idx])
        },
        character(1),
        USE.NAMES = FALSE
    )
}

#' Assign Condition Column Based on Controls
#' @param data Data frame containing treatment information.
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. Default is
#'   \code{NA}.
#' @return
#' The input data frame with a new \code{condition} column.
#' @keywords internal
assign_condition <- function(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name = NA
) {
    data$condition <- "treatment"
    data$condition[data[[treatment_column]] == negative_control_name] <-
        "negative_control"

    if (!is.na(positive_control_name)) {
        data$condition[data[[treatment_column]] == positive_control_name] <-
            "positive_control"
    }

    data
}

check_fit_data <- function(data, x_var, y_var, n_params) {
    missing_cols <- setdiff(c(x_var, y_var), names(data))
    if (length(missing_cols)) {
        stop_unable_to_fit(sprintf(
            "Column%s not found in data: %s",
            if (length(missing_cols) > 1) "s" else "",
            paste(missing_cols, collapse = ", ")
        ))
    }

    # Coerce as prepare_data() does; it will emit any coercion warning
    x <- suppressWarnings(as.numeric(data[[x_var]]))
    y <- suppressWarnings(as.numeric(data[[y_var]]))
    keep <- is.finite(x) & is.finite(y)

    if (!any(keep)) {
        stop_unable_to_fit(
            "No finite observations: every row has a missing or infinite value."
        )
    }
    if (!is.na(n_params) && sum(keep) < n_params) {
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

relocate_after <- function(data, col, after = "y") {
    if (!all(c(col, after) %in% names(data)) || identical(col, after)) {
        return(data)
    }

    target <- match(after, names(data))
    moving <- match(col, names(data))

    order <- seq_along(data)
    order <- order[order != moving]
    order <- append(order, moving, after = match(target, order))

    data[, order, drop = FALSE]
}

#' Copy the grouping variable into a standard `group` column
#'
#' @param data Prepared data.
#' @param shared_group Name of the column identifying groups.
#' @param ignore_group Skip the collision check, for data that already went
#'   through this once (a refit).
#' @return `data` with a character `group` column placed after `y`.
#' @keywords internal
add_group_column <- function(data, shared_group, ignore_group = FALSE) {
    if (!shared_group %in% colnames(data)) {
        stop_unable_to_fit(sprintf(
            "Grouping column not found in data: %s",
            shared_group
        ))
    }

    if (!ignore_group && "group" %in% colnames(data)) {
        message("Detected 'group' column in data, renaming to 'group_original'")
        data$group_original <- data$group
    }

    data$group <- as.character(data[[shared_group]])

    relocate_after(data, "group")
}

#' Build the group-expanded version of a curve function
#'
#' Parameters named in `shared_parameters` stay single; every other parameter
#' is given one copy per group, and the body dispatches on `group`.
#'
#' @param curve_function The single-group curve function.
#' @param group_vec Unique group labels, from the `group` column.
#' @param shared_parameters Names of the parameters held in common.
#' @return A function whose formals are the expanded parameter set.
#' @keywords internal
build_shared_function <- function(
    curve_function,
    group_vec,
    shared_parameters
) {
    shared_function <- function() {}

    formals(shared_function) <- as.pairlist(make_shared_formals(
        curve_function,
        group_vec,
        shared_parameters
    ))
    body(shared_function) <- make_shared_body(
        curve_function,
        group_vec,
        shared_parameters
    )

    shared_function
}

#' Build the argument list for `nlsLM()`
#'
#' @param curve_function The curve function referenced by the formula.
#' @param start_values Named list of starting values.
#' @param lower_bounds,upper_bounds Named lists of parameter bounds, or `NA`.
#' @return A list with `formula`, `start` and, if any bounds were given,
#'   `lower` and `upper` ordered to match `start`.
#' @keywords internal
create_function_arguments <- function(
    curve_function,
    start_values,
    lower_bounds = NA,
    upper_bounds = NA
) {
    final_arguments <- list(
        formula = make_formula(curve_function),
        start = unlist(start_values)
    )

    if (!all(is.na(lower_bounds)) || !all(is.na(upper_bounds))) {
        bounds <- make_bounds(curve_function, lower_bounds, upper_bounds)

        parameter_order <- names(final_arguments$start)
        if (!is.null(bounds$lower)) {
            bounds$lower <- bounds$lower[parameter_order]
        }
        if (!is.null(bounds$upper)) {
            bounds$upper <- bounds$upper[parameter_order]
        }

        final_arguments <- append(final_arguments, bounds)
    }

    final_arguments[!vapply(final_arguments, is.null, logical(1))]
}

#' Expand user-supplied bounds to every parameter of a curve function
#'
#' Unbounded parameters get `-Inf`/`Inf`. A bound named after a parameter
#' applies to it directly; for a group-expanded function (see
#' [build_shared_function()]) a bound named after a base parameter, e.g.
#' `bottom`, also applies to every per-group copy (`bottom_A`, `bottom_B`, ...).
#' Exact names take precedence over this expansion.
#'
#' @param curve_function The (possibly group-expanded) curve function.
#' @param lower_bounds,upper_bounds Named lists of single numeric bounds, or
#'   `NA` for none.
#' @return A list with named numeric vectors `lower` and `upper`, one element
#'   per parameter.
#' @keywords internal
make_bounds <- function(curve_function, lower_bounds = NA, upper_bounds = NA) {
    params <- setdiff(names(formals(curve_function)), c("x", "group"))

    expand <- function(bounds, default, label) {
        out <- stats::setNames(rep(default, length(params)), params)

        if (!is.list(bounds) && all(is.na(bounds))) {
            return(out)
        }

        bound_names <- names(bounds)
        if (is.null(bound_names) || any(bound_names == "")) {
            stop(sprintf("`%s_bounds` must be a fully named list.", label))
        }
        if (!all(lengths(bounds) == 1L)) {
            stop(sprintf(
                "Each element of `%s_bounds` must be a single number.",
                label
            ))
        }

        exact <- bound_names %in% params
        prefix_targets <- lapply(bound_names[!exact], function(nm) {
            params[startsWith(params, paste0(nm, "_"))]
        })

        # Every name must match; a single typo is an error, not ignored
        unknown <- bound_names[!exact][lengths(prefix_targets) == 0L]
        if (length(unknown)) {
            stop(sprintf(
                "Parameter bound(s) '%s' not found in model function. Available parameters: %s",
                paste(unknown, collapse = ", "),
                paste(params, collapse = ", ")
            ))
        }

        # Expanded (per-group) bounds first, so exact names win
        for (i in seq_along(prefix_targets)) {
            out[prefix_targets[[i]]] <- as.numeric(bounds[[
                bound_names[!exact][i]
            ]])
        }
        for (nm in bound_names[exact]) {
            out[nm] <- as.numeric(bounds[[nm]])
        }

        out
    }

    list(
        lower = expand(lower_bounds, -Inf, "lower"),
        upper = expand(upper_bounds, Inf, "upper")
    )
}


#' Default option lists
#'
#' @description
#' The single source of truth for option defaults. [resolve_model_options()] merges
#' user-supplied options over these.
#' @keywords internal
dosefitr_defaults <- list(
    curve = list(
        model = NA,
        curve_function = NA,
        start_function = NA,
        start_values = NA,
        lower_bounds = NA,
        upper_bounds = NA,
        loess_span = 0.3
    ),
    shared = list(
        shared_group = NA,
        shared_parameters = NA,
        return_func = FALSE
    ),
    outlier = list(
        huber = FALSE,
        huber_iter = 100L,
        huber_k = 1.345,
        huber_tol = 1e-6,
        rout = FALSE,
        rout_q = 1e-3,
        rout_scale = "mad"
    )
)


make_formula <- function(func) {
    if (!is.function(func)) {
        stop("`func` must be a function; got ", class(func)[1], ".")
    }

    function_out <- utils::capture.output(args(func))
    function_out <- function_out[function_out != "NULL"]
    function_out[-1] <- trimws(function_out[-1])
    function_out <- paste(function_out, collapse = "")
    function_out <- sub("^function ", "f", function_out)
    function_out <- paste("y ~", function_out)

    str2lang(function_out)
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
    nls_model,
    lower_x,
    upper_x,
    group = NA,
    num_points = 1e3
) {
    x_values <- seq(lower_x, upper_x, length.out = num_points)

    if (!any(is.na(group))) {
        # For each group, predict and store results
        res_list <- lapply(group, function(i) {
            y <- stats::predict(
                nls_model,
                newdata = data.frame(x = x_values, group = i)
            )
            data.frame(x = x_values, y = y, group = i, stringsAsFactors = FALSE)
        })
        # Combine all into one data.frame
        do.call(rbind, res_list)
    } else {
        y <- stats::predict(nls_model, newdata = data.frame(x = x_values))
        data.frame(x = x_values, y = y, stringsAsFactors = FALSE)
    }
}

find_x_for_y <- function(nls_model, x_values, target) {
    stats::uniroot(
        f = function(x) {
            stats::predict(nls_model, data.frame(x = x)) - target
        },
        interval = c(min(x_values), max(x_values))
    )
}


#' Predict a curve for one treatment-concentration record
#'
#' @param nls_model The fitted model, or `NULL`/`NA` if the record was not
#'   fitted.
#' @param data The record's data.
#' @param concentration_column,treatment_column Column names copied onto the
#'   predictions.
#' @param range Length-2 x range to predict over; defaults to the range of
#'   `data$x`.
#' @param on_missing For excluded or unfitted records: `"na"` returns `NA`,
#'   `"zero"` returns a flat prediction of 0 marked `exclude = TRUE`.
#' @param extra Named list of extra columns to add.
#' @return A data frame of predictions, or `NA`.
#' @keywords internal
predict_for_group <- function(
    nls_model,
    data,
    concentration_column,
    treatment_column,
    range = NULL,
    on_missing = c("na", "zero"),
    extra = list()
) {
    on_missing <- match.arg(on_missing)

    if (is.null(range)) {
        range <- range(data$x, na.rm = TRUE)
    }

    conc <- unique(data[[concentration_column]])
    treat <- unique(data[[treatment_column]])

    is_excluded <- "exclude" %in% names(data) && all(data$exclude %in% TRUE)
    is_unfitted <- is.null(nls_model) ||
        (length(nls_model) == 1L && all(is.na(nls_model)))

    if (is_excluded || is_unfitted) {
        if (on_missing == "na") {
            return(NA)
        }
        predicted <- data.frame(
            x = seq(range[1], range[2], length.out = 1e3),
            y = 0
        )
        predicted$exclude <- TRUE
    } else {
        predicted <- predict_data(
            nls_model = nls_model,
            lower_x = range[1],
            upper_x = range[2]
        )
        if ("exclude" %in% names(data)) {
            predicted$exclude <- FALSE
        }
    }

    predicted[[concentration_column]] <- conc
    predicted[[treatment_column]] <- treat

    for (nm in names(extra)) {
        predicted[[nm]] <- extra[[nm]]
    }

    predicted
}

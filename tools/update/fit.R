#' Attempt Fit
#' @description
#' Attempt to fit a nonlinear model using `minpack.lm::nlsLM()`.
#' @param data A data frame with columns `x` and `y`.
#' @param args A named list of arguments to pass to `minpack.lm::nlsLM()`.
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @importFrom minpack.lm nlsLM
#' @keywords internal
attempt_fit <- function(args, data, func) {
  warning_list <- character()
  error_list <- character()

  call_obj <- as.call(c(
    list(quote(minpack.lm::nlsLM)),
    list(data = quote(data)),
    args
  ))

  fit_env <- new.env(parent = parent.frame())
  fit_env$data <- data
  fit_env$f <- func

  result <- withCallingHandlers(
    tryCatch(
      list(obj = eval(call_obj, envir = fit_env)),
      error = function(e) {
        error_list <<- c(error_list, conditionMessage(e))
        list(obj = NA)
      }
    ),
    warning = function(w) {
      warning_list <<- c(warning_list, conditionMessage(w))
      tryInvokeRestart("muffleWarning")
    }
  )

  result$ok <- !length(error_list)
  result$warnings <- if (length(warning_list)) warning_list else NA_character_
  result$errors <- if (length(error_list)) error_list else NA_character_

  result
}


#' Make Formula
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @description
#' Generate a formula object in the form `y ~ f(arguments)` from a function.
#' @return
#' A language object representing the formula `y ~ f(parameters)`.
#' @keywords internal
make_formula <- function(func) {
  # Capture printed args output
  function_out <- capture.output(args(func))
  function_out <- function_out[function_out != "NULL"]

  # Trim white space from all but the first line
  function_out[-1] <- trimws(function_out[-1])

  # Collapse into a single string
  function_out <- paste(function_out, collapse = "")

  # Replace 'function ' with 'f'
  function_out <- sub("^function ", "f", function_out)

  # Prepend 'y ~'
  function_out <- paste("y ~", function_out)

  # Turn into a formula object
  return(str2lang(function_out))
}


#' Make Bounds
#' @description
#' Generate named lists of upper and lower bounds for `minpack.lm::nlsLM()`,
#'    setting specified bounds and defaults for unspecified parameters.
#' @param lower_x A named list specifying lower bounds for parameters.
#'   Unspecified parameters default to `-Inf`.
#' @param upper_x A named list specifying upper bounds for parameters.
#'   Unspecified parameters default to `Inf`.
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @return
#' A named list containing `lower` and `upper` elements with bounds
#'   for all parameters in `func`.
#' @details
#' This ensures all parameters in the model function have defined bounds,
#'    either from user specification or from defaults. Names must match
#'    those in the model function arguments.
#' @keywords internal
make_bounds <- function(func, lower_x = NA, upper_x = NA) {
  formals <- names(formals(func))
  formals <- formals[!formals %in% c("x", "group")]

  # If only one supplied - fix: use all(is.na()) for lists
  final_list <- list(lower = lower_x, upper = upper_x)
  final_list <- final_list[!sapply(final_list, function(x) all(is.na(x)))]

  # Check bounded args exist
  for (i in names(final_list)) {
    x <- final_list[[i]]
    if (!any(names(x) %in% formals)) {
      missing <- names(x)[!names(x) %in% formals]
      stop(sprintf(
        "Parameter bound(s) '%s' not found in model function. Available parameters: %s",
        paste(missing, collapse = ", "),
        paste(formals, collapse = ", ")
      ))
    }
  }

  # Bounds need to be in order
  result <- list()
  for (i in names(final_list)) {
    x <- final_list[[i]]

    # Create full vector with defaults
    if (i == "lower") {
      full <- rep(-Inf, length(formals))
    } else if (i == "upper") {
      full <- rep(Inf, length(formals))
    }
    names(full) <- formals

    # Replace with provided values - FIX: use name matching
    replace <- unlist(x)
    full[names(replace)] <- replace

    result[[i]] <- full
  }
  return(result)
}

#' Make Shared Start Values
#' @description
#' Generate appropriate starting parameter values when fitting models with
#'    shared parameters across experimental groups.
#' @param data A data frame containing data with group information.
#' @param model_func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param start_func A function for generating a maed list of starting parameter values.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to be shared.
#' @return
#' A named list of starting values.
#' @details
#' This function creates starting values by:
#' \itemize{
#'   \item Calculating global starting values.
#'   \item Calculating group-specific starting values for non-shared parameters.
#'   \item Properly naming group-specific parameters with group suffixes.
#'   \item Combining shared and group-specific values into a single list.
#' }
#' @keywords internal
make_shared_start_values <- function(
  data,
  model_func,
  start_func,
  group_vec,
  share_params
) {
  # Global starting values for shared parameters
  start_values_global <- start_func(
    x = data$x,
    y = data$y
  )

  # Group-specific starting values
  start_values_group <- lapply(group_vec, function(i) {
    group_data <- data[data$group == i, ]
    start_func(
      x = group_data$x,
      y = group_data$y
    )
  })
  names(start_values_group) <- group_vec

  # Modify names for group-specific parameters
  for (i in seq_along(start_values_group)) {
    group_name <- names(start_values_group)[i]
    x <- start_values_group[[i]]

    # Rename parameters that are not shared
    new_names <- ifelse(
      !names(x) %in% share_params,
      paste(names(x), "_", group_name, sep = ""),
      names(x)
    )
    names(x) <- new_names

    # Keep only group-specific parameters
    start_values_group[[i]] <- x[!names(x) %in% names(start_values_global)]
  }

  # Combine shared and group-specific values
  start_values_global <- start_values_global[
    names(start_values_global) %in% share_params
  ]

  # Flatten and combine
  group_values_flat <- unlist(unname(start_values_group))

  return(append(start_values_global, group_values_flat))
}

#' Make Shared Formals
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to be shared.
#' @description
#' Modify function formal arguments to support shared parameters across groups
#'    by creating group-specific parameter names for non-shared parameters.
#' @return
#' A pairlist object suitable for creating a new function with
#'   group-specific parameters.
#' @keywords internal
make_shared_formals <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  # Create group-specific parameter names
  append_arguments <- lapply(
    unique_arguments,
    function(arg) paste0(arg, "_", group_vec)
  )

  # Validate shared parameters exist
  missing_params <- share_params[!share_params %in% formal_arguments]
  if (length(missing_params) > 0) {
    stop(sprintf(
      "Shared parameter(s) '%s' not found in model function. Available parameters: %s",
      paste(missing_params, collapse = ", "),
      paste(formal_arguments, collapse = ", ")
    ))
  }

  # Combine all arguments
  final_arguments <- c(
    "x",
    formal_arguments[formal_arguments %in% share_params],
    unlist(append_arguments),
    "group"
  )

  # Create formals list
  final_formals <- rep(list(quote(expr = )), length(final_arguments))
  names(final_formals) <- final_arguments

  return(as.pairlist(final_formals))
}

#' Make Shared Body
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to be shared.
#' @description
#' Generate a function body that implements shared parameter logic using
#'    nested `ifelse()` functions.
#' @return
#' A call object representing the modified function body.
#' @keywords internal
make_shared_bod <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  # Create parameter mapping expressions
  append_arguments <- lapply(
    unique_arguments,
    function(arg) paste0(arg, "_", group_vec)
  )

  # Generate ifelse expressions for each parameter
  expression_list <- list()
  for (i in seq_along(unique_arguments)) {
    arg <- unique_arguments[i]
    group_params <- append_arguments[[i]]

    # Build nested ifelse expression
    if (length(group_vec) == 1) {
      # Simple case: only one group
      expr_string <- sprintf("%s <- %s", arg, group_params[1])
    } else {
      # Build nested ifelse chain from inside out
      expr <- group_params[length(group_params)]

      # Work backwards through the groups
      for (j in rev(seq_len(length(group_vec) - 1))) {
        expr <- sprintf(
          "ifelse(group == '%s', %s, %s)",
          group_vec[j],
          group_params[j],
          expr
        )
      }

      expr_string <- sprintf("%s <- %s", arg, expr)
    }

    expression_list[[i]] <- str2lang(expr_string)
  }

  # Combine with original function body
  return(as.call(c(as.name("{"), expression_list, body(func))))
}

#' Prepared Shared Parameters
#' @description
#' Execute the complete pipeline for creating shared parameter models including
#'    data preparation, function modification, and starting value generation.
#' @param data Input data frame.
#' @param share_params Character vector of parameter names to be shared.
#' @param share_group Character vector of unique group identifiers.
#' @param func Function defining the mathematical model.
#' @param start_func Function for generating starting parameter values.
#' @param model_func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param ignore_group TODO: description.
#' @return
#' A list containing:
#'  \itemize{
#'    \item Modified data frame with standardized group column
#'    \item Modified model function with shared parameter support
#'    \item Starting parameter values for the shared parameter model
#'  }
#' @details
#' This function orchestrates the complete shared parameter workflow.
#' @keywords internal
prep_shared_params <- function(
  data,
  share_params,
  share_group,
  model_func,
  start_func,
  ignore_group = FALSE
) {
  # Handle existing 'group' column
  if (!ignore_group) {
    if ("group" %in% colnames(data)) {
      message("Detected 'group' column in data, renaming to 'group_original'")
      data$group_original <- data$group
    }
  }
  # Add group column
  data$group <- as.character(data[[share_group]])

  # Relocate group column after y
  y_index <- which(names(data) == "y")
  if (length(y_index) > 0) {
    col_order <- seq_len(ncol(data))
    group_index <- which(names(data) == "group")
    col_order <- col_order[col_order != group_index]
    col_order <- append(col_order, group_index, after = y_index)
    data <- data[, col_order]
  }

  group_vec <- unique(data$group)

  start_values <- make_shared_start_values(
    data,
    model_func,
    start_func,
    group_vec,
    share_params
  )
  shared_arguments <- make_shared_formals(
    model_func,
    group_vec,
    share_params
  )
  shared_body <- make_shared_bod(model_func, group_vec, share_params)

  # Create modified function
  model_func <- function() {}
  formals(model_func) <- as.pairlist(shared_arguments)
  body(model_func) <- shared_body

  return(
    list(data = data, model_func = model_func, start_values = start_values)
  )
}

#' Prepare Arguments
#' @description
#' Assemble all necessary arguments for `minpack.lm::nlsLM()` including formula,
#'    starting values, and optional parameter bounds.
#' @param lower_bounds Optional named list of lower parameter bounds.
#' @param upper_bounds Optional named list of upper parameter bounds.
#' @param model_func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#' @param start_values A named list of starting parameter values.
#' @return
#' A named list to be injected into `minpack.lm::nlsLM()` containing:
#'  \describe{
#'    \item{formula}{Model formula object}
#'    \item{start}{Starting parameter values}
#'    \item{lower}{Lower parameter bounds (if specified)}
#'    \item{upper}{Upper parameter bounds (if specified)}
#'  }
#' @keywords internal
prep_args <- function(
  model_func,
  start_values,
  lower_bounds = NA,
  upper_bounds = NA
) {
  # Formula
  func_formula <- make_formula(model_func)

  final_args <- list(
    formula = func_formula,
    start = unlist(start_values)
  )

  # Bounds
  if (!all(is.na(lower_bounds)) | !all(is.na(upper_bounds))) {
    bounds <- make_bounds(model_func, lower_bounds, upper_bounds)

    # CRITICAL: Reorder bounds to match start parameter order
    param_order <- names(final_args$start)
    if (!is.null(bounds$lower)) {
      bounds$lower <- bounds$lower[param_order]
    }
    if (!is.null(bounds$upper)) {
      bounds$upper <- bounds$upper[param_order]
    }

    final_args <- append(final_args, bounds)
  }

  # print("Final args being passed to nlsLM:")
  # print(final_args)
  return(final_args[!sapply(final_args, is.null)])
}


#' Set Default Options
#' @param curve_opts A named list of curve fitting options. Elements include:
#'   \itemize{
#'     \item `model`: Character string specifying a built-in model from \code{dosefitr_models} or "loess".
#'     \item `model_func`: A function that describes a curve in terms of x (used if \code{model} is NA).
#'     \item `start_func`: A function that generates a named list of starting values for `model_func` (used if \code{model} is NA).
#'     \item `start_values`: A named list of starting values for `model_func`.
#'     \item `lower_bounds`: A named list that specifies the lower bounds for parameters in `model_func`.
#'     \item `upper_bounds`:A named list that specifies the upper bounds for parameters in `model_func`.
#'   }
#' @param shared_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `share_group`: Character string specifying the column name for the grouping variable (default: NA).
#'     \item `share_params`: Character string specifying which parameters in the model are to be shared (default: NA).
#'     \item `return_func`: Logical; return fitted function? (default: FALSE).
#'   }
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
#' Set Default Options for Curve Fitting. Merges user-supplied options with
#'    default values for curve fitting, shared parameters, and outlier handling.
#' @keywords internal
set_default_opts <- function(
  curve_opts = NA,
  shared_opts = NA,
  outlier_opts = NA
) {
  # Define default lists
  curve_defaults <- list(
    model = NA,
    model_func = NA,
    start_func = NA,
    start_values = NA,
    lower_bounds = NA,
    upper_bounds = NA,
    loess_span = 0.3
  )

  shared_defaults <- list(
    share_group = NA,
    share_params = NA,
    return_func = FALSE
  )

  outlier_defaults <- list(
    huber = FALSE,
    huber_iter = 100,
    huber_k = 1.345,
    huber_tol = 1e-6,
    rout = FALSE,
    rout_q = 1e-3,
    rout_scale = "mad"
  )

  # Ensure NAs become empty lists
  curve_opts <- if (all(is.na(curve_opts))) list() else curve_opts
  shared_opts <- if (all(is.na(shared_opts))) list() else shared_opts
  outlier_opts <- if (all(is.na(outlier_opts))) list() else outlier_opts

  # Merge defaults
  curve_opts <- modifyList(curve_defaults, curve_opts)
  shared_opts <- modifyList(shared_defaults, shared_opts)
  outlier_opts <- modifyList(outlier_defaults, outlier_opts)

  # Return all together as a named list
  return(
    list(
      curve_opts = curve_opts,
      shared_opts = shared_opts,
      outlier_opts = outlier_opts
    )
  )
}


#' Fit a Curve to Data
#' @description
#' Fits a non-linear model to x-y data with support for built-in models,
#'    custom model functions, shared parameters across groups, and robust
#'    outlier handling via Huber regression and/or ROUT detection.
#' @param data A data frame containing the variables to fit.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param curve_opts A named list of curve fitting options. Elements include:
#'   \itemize{
#'     \item `model`: Character string specifying a built-in model from \code{dosefitr_models} or "loess".
#'     \item `model_func`: A function that describes a curve in terms of x (used if \code{model} is NA).
#'     \item `start_func`: A function that generates a named list of starting values for `model_func` (used if \code{model} is NA).
#'     \item `start_values`: A named list of starting values for `model_func`.
#'     \item `lower_bounds`: A named list that specifies the lower bounds for parameters in `model_func`.
#'     \item `upper_bounds`:A named list that specifies the upper bounds for parameters in `model_func`.
#'   }
#' @param shared_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `share_group`: Character string specifying the column name for the grouping variable (default: NA).
#'     \item `share_params`: Character string specifying which parameters in the model are to be shared (default: NA).
#'     \item `return_func`: Logical; return fitted function? (default: FALSE).
#'   }
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
#' A list containing:
#'  \itemize{
#'    \item `data`: The prepared data frame, with an \code{outlier_<y_var>} column added if ROUT detection was used.
#'    \item `fit:` A named list with the fitted object `obj`, in addition to `warnings` and `errors` that occured during the fit.
#'    \item `func`: The model function (only if \code{shared_opts$return_func = TRUE}).
#'  }
#' @examples
#' \dontrun{
#' # Fit a built-in model
#' result <- fit_curve(
#'   data = my_data,
#'   x_var = "time",
#'   y_var = "response",
#'   curve_opts = list(model = "exponential")
#' )
#' # Fit with ROUT outlier detection
#' result <- fit_curve(
#'   data = my_data,
#'   x_var = "time",
#'   y_var = "response",
#'   curve_opts = list(model = "logistic"),
#'   outlier_opts = list(rout = TRUE, rout_q = 0.01)
#' )
#' }
#' @export
fit_curve <- function(
  data,
  x_var,
  y_var,
  curve_opts = NA,
  shared_opts = NA,
  outlier_opts = NA
) {
  default_opts <- set_default_opts(
    curve_opts = curve_opts,
    shared_opts = shared_opts,
    outlier_opts = outlier_opts
  )
  curve_opts <- default_opts$curve_opts
  shared_opts <- default_opts$shared_opts
  outlier_opts <- default_opts$outlier_opts

  is_shared <- !is.na(shared_opts$share_group) &&
    !any(is.na(shared_opts$share_params))
  is_outlier <- isTRUE(outlier_opts$huber) || isTRUE(outlier_opts$rout)

  # Built-in models
  if (!is.na(curve_opts$model)) {
    if (curve_opts$model == "loess") {
      model_func <- NA
      start_func <- NA
    } else {
      if (!curve_opts$model %in% names(dosefitr_models)) {
        stop(sprintf("Not a built-in model: %s", curve_opts$model))
      }
      model_func <- dosefitr_models[[curve_opts$model]]$model_func
      start_func <- dosefitr_models[[curve_opts$model]]$start_func
    }
  } else {
    model_func <- curve_opts$model_func
    start_func <- curve_opts$start_func
  }

  check_fit_data(
    data = data,
    x_var = x_var,
    y_var = y_var,
    n_params = length(formals(model_func)) - 1
  )

  # Prepare
  data <- prep_data(
    data = data,
    x_var = x_var,
    y_var = y_var
  )

  if (is.na(curve_opts$start_values)) {
    start_values <- start_func(x = data$x, y = data$y)
  } else {
    start_values <- curve_opts$start_values
  }

  if (is_shared) {
    params <- prep_shared_params(
      data = data,
      share_params = shared_opts$share_params,
      share_group = shared_opts$share_group,
      model_func = model_func,
      start_func = start_func
    )
  } else {
    params <- list(
      data = data,
      model_func = model_func,
      start_values = start_values
    )
  }

  # Create arguments
  args <- prep_args(
    model_func = params$model_func,
    start_values = params$start_values,
    lower_bounds = curve_opts$lower_bounds,
    upper_bounds = curve_opts$upper_bounds
  )

  fit <- attempt_fit(
    args = args,
    data = params$data,
    func = params$model_func
  )

  if (!fit$ok) {
    errors <- unique(unlist(fit$errors))
    hints <- unique(Filter(Negate(is.null), translate_fit_error(errors)))

    msg <- if (length(hints)) {
      paste0(
        "Unable to fit a model to this data:\n",
        paste0(hints, collapse = "\n")
      )
    } else {
      "Unable to fit a model to this data. Use conditionCall() or the 'errors' element of this condition for the underlying message."
    }

    stop_unable_to_fit(msg, errors = errors)
  }

  # if (!is.na(fit$errors)) {
  #   stop(paste("Error in fit:", fit$error, sep = " "))
  # }

  if (is_outlier) {
    if (outlier_opts$huber) {
      fit$obj <- tryCatch(
        huber_irwls(
          obj = fit$obj,
          args = args,
          data = params$data,
          model_func = params$model_func,
          outlier_opts = outlier_opts
        ),
        error = function(e) {
          warning(
            "Huber refinement failed; using the unweighted fit.",
            call. = FALSE
          )
          fit$obj
        }
      )
    }
    if (outlier_opts$rout) {
      outlier_indices <- find_rout_indices(
        obj = fit$obj,
        outlier_opts = outlier_opts
      )

      outlier_column <- paste("outlier_", y_var, sep = "")
      params$data[[outlier_column]] <- seq_len(nrow(data)) %in% outlier_indices

      # Relocate outlier column after y
      y_index <- which(names(params$data) == "y")
      if (length(y_index) > 0) {
        col_order <- seq_len(ncol(params$data))
        outlier_index <- which(names(params$data) == outlier_column)
        col_order <- col_order[col_order != outlier_index]
        col_order <- append(col_order, outlier_index, after = y_index)
        params$data <- params$data[, col_order]
      }

      # If no outliers, move on
      if (!all(params$data[[outlier_column]] == FALSE)) {
        # Refit without outliers
        data_filtered <- params$data[params$data[[outlier_column]] == FALSE, ]

        # Generate new starting values
        if (is_shared) {
          params_filtered <- prep_shared_params(
            data = data_filtered,
            share_params = shared_opts$share_params,
            share_group = shared_opts$share_group,
            model_func = model_func,
            start_func = start_func,
            ignore_group = TRUE
          )
        } else {
          params_filtered <- list(
            data = data_filtered,
            model_func = model_func,
            start_values = start_values
          )
        }

        args_filtered <- prep_args(
          model_func = params_filtered$model_func,
          start_values = params_filtered$start_values,
          lower_bounds = curve_opts$lower_bounds,
          upper_bounds = curve_opts$upper_bounds
        )

        fit_filtered <- attempt_fit(
          args = args_filtered,
          data = params_filtered$data,
          func = params_filtered$model_func
        )

        if (!fit_filtered$ok) {
          message(
            "Could not refit with removal of outliers; returning the original fit."
          )
          message(paste(fit_filtered$errors, collapse = "; "))
        } else {
          fit$obj <- fit_filtered$obj
        }
      }
    }
  }

  if (shared_opts$return_func) {
    return(list(data = params$data, fit = fit, func = params$model_func))
  } else {
    return(list(data = params$data, fit = fit))
  }
}

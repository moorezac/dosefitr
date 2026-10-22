is_scalar_str <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

#' Merge user-supplied options over the package defaults
#'
#' Unrecognised names are an error rather than being silently carried through,
#' so a mistyped option fails loudly instead of leaving the default in place.
#'
#' @param curve_options,shared_options,outlier_options Lists of user options,
#'   or `NA` to use the defaults.
#' @return A list with `curve`, `shared` and `outlier` option lists.
#' @keywords internal
resolve_model_options <- function(
  curve_options = NA,
  shared_options = NA,
  outlier_options = NA
) {
  resolve <- function(user_options, default_options, label) {
    if (!is.list(user_options)) {
      return(default_options)
    }

    unknown <- setdiff(names(user_options), names(default_options))
    if (length(unknown)) {
      stop(sprintf(
        "Unknown %s option%s: %s",
        label,
        if (length(unknown) > 1) "s" else "",
        paste(unknown, collapse = ", ")
      ))
    }

    utils::modifyList(default_options, user_options)
  }

  list(
    curve = resolve(curve_options, dosefitr_defaults$curve, "curve"),
    shared = resolve(shared_options, dosefitr_defaults$shared, "shared"),
    outlier = resolve(outlier_options, dosefitr_defaults$outlier, "outlier")
  )
}

#' Resolve curve and starting-value functions from curve options
#'
#' @param curve_options Resolved curve options from [resolve_model_options()].
#' @return A list with `curve_function` and `start_function`. Both are `NULL`
#'   for `loess`, which is fitted without a parametric curve.
#' @keywords internal
resolve_model_functions <- function(curve_options) {
  model_name <- curve_options$model

  if (!is_scalar_str(model_name)) {
    if (!is.function(curve_options$curve_function)) {
      stop(
        "Supply either `model` (a built-in model name or \"loess\") or ",
        "`curve_function` in `curve_options`."
      )
    }
    return(curve_options[c("curve_function", "start_function")])
  }

  if (identical(model_name, "loess")) {
    return(list(
      curve_function = NULL,
      start_function = NULL
    ))
  }

  if (!model_name %in% names(dosefitr_models)) {
    stop(sprintf("Not a built-in model: %s", model_name))
  }

  dosefitr_models[[model_name]][c("curve_function", "start_function")]
}

#' Should parameters be shared across groups?
#'
#' Sharing needs both a grouping column and at least one parameter to share.
#' `NULL` and `NA` both count as "not supplied".
#'
#' @param shared_options Resolved shared options from [resolve_model_options()].
#' @return `TRUE` or `FALSE`.
#' @keywords internal
use_shared_parameters <- function(shared_options) {
  is_scalar_str(shared_options$shared_group) &&
    length(shared_options$shared_parameters) > 0 &&
    !any(is.na(shared_options$shared_parameters))
}

make_shared_formals <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  append_arguments <- lapply(
    unique_arguments,
    function(arg) paste0(arg, "_", group_vec)
  )

  missing_params <- share_params[!share_params %in% formal_arguments]
  if (length(missing_params) > 0) {
    stop(sprintf(
      "Shared parameter(s) '%s' not found in model function. Available parameters: %s",
      paste(missing_params, collapse = ", "),
      paste(formal_arguments, collapse = ", ")
    ))
  }

  final_arguments <- c(
    "x",
    formal_arguments[formal_arguments %in% share_params],
    unlist(append_arguments),
    "group"
  )
  final_formals <- rep(list(quote(expr = )), length(final_arguments))
  names(final_formals) <- final_arguments

  as.pairlist(final_formals)
}

#' Assemble everything needed for a single `nlsLM()` call
#'
#' When parameters are shared, the data gains a `group` column, the curve is
#' replaced by its group-expanded version, and starting values are derived per
#' group -- so `start_values` passed in is not used on that path.
#'
#' @param data Prepared data with `x` and `y` columns.
#' @param model_functions Resolved functions from [resolve_model_functions()].
#' @param start_values Starting values, or `NULL` for a model without them.
#' @param model_options Resolved options from [resolve_model_options()].
#' @param use_shared_params Whether to share parameters across groups.
#' @param ignore_group Skip the `group` column collision check (used on refit,
#'   where the column has already been renamed).
#' @param regenerate_start Recompute starting values from `data` rather than
#'   reusing the ones passed in.
#' @return A list with `data`, `curve_function`, `start_values` and
#'   `function_arguments`, ready to hand to [fit_with_nls()].
#' @keywords internal
prepare_nls <- function(
  data,
  model_functions,
  start_values,
  model_options,
  use_shared_params,
  ignore_group = FALSE,
  regenerate_start = FALSE
) {
  if (use_shared_params) {
    shared_parameters <- model_options$shared$shared_parameters
    data <- add_group_column(
      data = data,
      shared_group = model_options$shared$shared_group,
      ignore_group = ignore_group
    )
    group_vec <- unique(data$group)

    curve_function <- build_shared_function(
      curve_function = model_functions$curve_function,
      group_vec = group_vec,
      shared_parameters = shared_parameters
    )
    start_values <- make_shared_start_values(
      data,
      model_functions$curve_function,
      model_functions$start_function,
      group_vec,
      shared_parameters
    )
  } else {
    curve_function <- model_functions$curve_function

    if (regenerate_start && is.function(model_functions$start_function)) {
      start_values <- model_functions$start_function(
        x = data$x,
        y = data$y
      )
    }
  }

  list(
    data = data,
    curve_function = curve_function,
    start_values = start_values,
    function_arguments = create_function_arguments(
      curve_function = curve_function,
      start_values = start_values,
      lower_bounds = model_options$curve$lower_bounds,
      upper_bounds = model_options$curve$upper_bounds
    )
  )
}

#' Attempt a single `nlsLM()` fit, capturing warnings and errors
#'
#' @param function_arguments Argument list from [create_function_arguments()].
#' @param data Data to fit on.
#' @param curve_function The curve function referenced by the formula.
#' @return A list with `nls_model`, `ok`, `warnings` and `errors`.
#' @keywords internal
fit_with_nls <- function(function_arguments, data, curve_function) {
  warning_list <- character()
  error_list <- character()

  final_call <- as.call(c(
    list(quote(minpack.lm::nlsLM)),
    list(data = quote(data)),
    function_arguments
  ))

  fit_env <- new.env(parent = parent.frame())
  fit_env$data <- data
  fit_env$f <- curve_function

  result <- withCallingHandlers(
    tryCatch(
      list(nls_model = eval(final_call, envir = fit_env)),
      error = function(e) {
        error_list <<- c(error_list, conditionMessage(e))
        list(nls_model = NA)
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

#' Refine a nls_model with Huber IRWLS weighting
#'
#' @param nls_model The fitted model object.
#' @param function_arguments,data,curve_function,outlier_options Passed to
#'   [huber_irwls()].
#' @return The refined model, or the original one if refinement failed.
#' @keywords internal
refine_with_huber <- function(
  nls_model,
  function_arguments,
  data,
  curve_function,
  outlier_options
) {
  tryCatch(
    huber_irwls(
      nls_model = nls_model,
      function_arguments = function_arguments,
      data = data,
      curve_function = curve_function,
      outlier_options = outlier_options
    ),
    error = function(e) {
      warning(
        "Huber refinement failed; using the unweighted fit.",
        call. = FALSE
      )
      nls_model
    }
  )
}

make_shared_body <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  # For each per-group parameter, build
  #   arg <- ifelse(group == g1, arg_g1, ifelse(group == g2, arg_g2, arg_gN))
  # directly as language objects, so group labels containing spaces, quotes
  # or other non-syntactic characters cannot break parsing.
  expression_list <- lapply(unique_arguments, function(arg) {
    group_params <- paste0(arg, "_", group_vec)
    n_groups <- length(group_params)

    expr <- as.name(group_params[n_groups])
    for (j in rev(seq_len(n_groups - 1L))) {
      expr <- call(
        "ifelse",
        call("==", as.name("group"), group_vec[j]),
        as.name(group_params[j]),
        expr
      )
    }

    call("<-", as.name(arg), expr)
  })

  # Combine with original function body
  as.call(c(as.name("{"), expression_list, body(func)))
}

#' Flag ROUT outliers and refit without them
#'
#' @param model The fit returned by [fit_with_nls()].
#' @param data Prepared data to flag and refit on.
#' @param model_functions,start_values,model_options,use_shared_params Passed
#'   to [prepare_nls()] for the refit.
#' @param has_user_start Whether the user supplied starting values.
#' @param y_var Name of the response variable, used for the outlier column.
#' @return A list with the updated `model` and the flagged `data`.
#' @keywords internal
refine_with_rout <- function(
  model,
  data,
  model_functions,
  start_values,
  model_options,
  use_shared_params,
  has_user_start,
  y_var
) {
  indices <- find_rout_indices(
    nls_model = model$nls_model,
    outlier_options = model_options$outlier
  )
  column <- paste0("outlier_", y_var)

  data[[column]] <- seq_len(nrow(data)) %in% indices
  data <- relocate_after(data, column)

  if (!any(data[[column]])) {
    return(list(model = model, data = data))
  }

  filtered <- data[!data[[column]], , drop = FALSE]

  if ("group" %in% names(data) && !setequal(filtered$group, data$group)) {
    message(
      "Outlier removal would drop an entire group; returning the original fit."
    )
    return(list(model = model, data = data))
  }

  remodel_inputs <- prepare_nls(
    data = filtered,
    model_functions = model_functions,
    start_values = start_values,
    model_options = model_options,
    use_shared_params = use_shared_params,
    ignore_group = TRUE,
    regenerate_start = !has_user_start
  )

  remodel <- fit_with_nls(
    function_arguments = remodel_inputs$function_arguments,
    data = remodel_inputs$data,
    curve_function = remodel_inputs$curve_function
  )

  if (remodel$ok) {
    model <- remodel
  } else {
    message(
      "Could not refit with removal of outliers; returning the original fit."
    )
    message(paste(remodel$errors, collapse = "; "))
  }

  list(model = model, data = data)
}

make_shared_start_values <- function(
  data,
  curve_function,
  start_function,
  group_vec,
  share_parameters
) {
  start_values_global <- start_function(
    x = data$x,
    y = data$y
  )
  start_values_group <- lapply(group_vec, function(i) {
    group_data <- data[data$group == i, ]
    start_function(
      x = group_data$x,
      y = group_data$y
    )
  })
  names(start_values_group) <- group_vec

  for (i in seq_along(start_values_group)) {
    group_name <- names(start_values_group)[i]
    x <- start_values_group[[i]]
    new_names <- ifelse(
      !names(x) %in% share_parameters,
      paste(names(x), "_", group_name, sep = ""),
      names(x)
    )
    names(x) <- new_names

    start_values_group[[i]] <- x[!names(x) %in% names(start_values_global)]
  }
  start_values_global <- start_values_global[
    names(start_values_global) %in% share_parameters
  ]
  group_values_flat <- unlist(unname(start_values_group))

  append(start_values_global, group_values_flat)
}

#' Fit a loess smoother in place of a parametric curve
#'
#' @param data Prepared data with `x` and `y` columns.
#' @param model_options Resolved options from [resolve_model_options()].
#' @return A list shaped like the output of [fit_curve()], with the `loess`
#'   object stored as `curve$nls_model`.
#' @keywords internal
fit_loess_curve <- function(data, model_options) {
  if (
    use_shared_parameters(model_options$shared) ||
      isTRUE(model_options$outlier$huber) ||
      isTRUE(model_options$outlier$rout)
  ) {
    stop(
      "Shared parameters and outlier handling are not supported with ",
      "model = \"loess\"."
    )
  }

  warning_list <- character()

  loess_model <- withCallingHandlers(
    tryCatch(
      stats::loess(
        y ~ x,
        data = data,
        span = model_options$curve$loess_span
      ),
      error = function(e) {
        stop_unable_to_fit(
          no_fit_message(conditionMessage(e)),
          errors = conditionMessage(e)
        )
      }
    ),
    warning = function(w) {
      warning_list <<- c(warning_list, conditionMessage(w))
      tryInvokeRestart("muffleWarning")
    }
  )

  list(
    data = data,
    curve = list(
      nls_model = loess_model,
      ok = TRUE,
      warnings = if (length(warning_list)) warning_list else NA_character_,
      errors = NA_character_
    )
  )
}

#' Fit a Curve to Data
#' @description
#' Fits a built-in or user-supplied curve to one set of observations with
#' \code{minpack.lm::nlsLM()}, optionally sharing parameters across groups and
#' refining the fit with Huber weighting and/or ROUT outlier removal.
#' @param data A data frame.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param curve_options A named list of curve fitting options, or \code{NA} for
#'   the defaults. Elements include:
#'   \itemize{
#'     \item `model`: Character string naming a built-in model from
#'       \code{dosefitr_models}, or \code{"loess"}.
#'     \item `curve_function`: A function describing the curve in terms of `x`
#'       (used if `model` is \code{NA}).
#'     \item `start_function`: A function of `x` and `y` returning a named list
#'       of starting values for `curve_function` (used if `model` is
#'       \code{NA}).
#'     \item `start_values`: A named list of starting values.
#'     \item `lower_bounds`: A named list of lower bounds for parameters.
#'     \item `upper_bounds`: A named list of upper bounds for parameters.
#'     \item `loess_span`: Span used when `model` is \code{"loess"}.
#'   }
#' @param shared_options A named list of options for sharing parameters across
#'   groups, or \code{NA} for the defaults. Elements include:
#'   \itemize{
#'     \item `shared_group`: Column name identifying the groups.
#'     \item `shared_parameters`: Character vector of parameters held in
#'       common across groups; all others are fitted per group.
#'     \item `return_func`: Logical; if \code{TRUE}, the group-expanded curve
#'       function is returned as `func`.
#'   }
#' @param outlier_options A named list of outlier-handling options, or
#'   \code{NA} for the defaults. Elements include `huber`, `huber_iter`,
#'   `huber_k` and `huber_tol` (Huber IRWLS refinement), and `rout`, `rout_q`
#'   and `rout_scale` (ROUT outlier removal).
#' @return
#' A list containing:
#'  \itemize{
#'    \item `data`: The prepared data, with `x` and `y` columns, a `group`
#'      column when parameters are shared, and an `outlier_<y_var>` column
#'      when ROUT is used.
#'    \item `curve`: A list with the fitted model (`nls_model`), `ok`, and any
#'      captured `warnings` and `errors`.
#'    \item `func`: The group-expanded curve function, only if
#'      `return_func = TRUE`.
#'  }
#' @details
#' Unrecognised option names are an error. If no model can be fitted, an error
#' of class \code{dosefitr_fit_error} is raised; its `errors` element holds the
#' underlying messages.
#' @seealso \code{\link{plot_model}}, \code{\link{fit_concentrations}}
#' @export
fit_curve <- function(
  data,
  x_var,
  y_var,
  curve_options = NA,
  shared_options = NA,
  outlier_options = NA
) {
  model_options <- resolve_model_options(
    curve_options = curve_options,
    shared_options = shared_options,
    outlier_options = outlier_options
  )
  model_functions <- resolve_model_functions(model_options$curve)

  check_fit_data(
    data = data,
    x_var = x_var,
    y_var = y_var,
    n_params = if (is.function(model_functions$curve_function)) {
      length(formals(model_functions$curve_function)) - 1
    } else {
      NA
    }
  )
  data <- prepare_data(data, x_var, y_var)

  if (identical(model_options$curve$model, "loess")) {
    return(fit_loess_curve(data, model_options))
  }

  has_user_start <- !all(is.na(model_options$curve$start_values))
  if (has_user_start) {
    start_values <- model_options$curve$start_values
  } else if (is.function(model_functions$start_function)) {
    start_values <- model_functions$start_function(data$x, data$y)
  } else {
    start_values <- NULL
  }

  use_shared_params <- use_shared_parameters(model_options$shared)

  model_inputs <- prepare_nls(
    data = data,
    model_functions = model_functions,
    start_values = start_values,
    model_options = model_options,
    use_shared_params = use_shared_params
  )

  model <- fit_with_nls(
    function_arguments = model_inputs$function_arguments,
    data = model_inputs$data,
    curve_function = model_inputs$curve_function
  )

  if (!model$ok) {
    errors <- model$errors
    stop_unable_to_fit(no_fit_message(errors), errors = errors)
  }

  fit_data <- model_inputs$data

  if (isTRUE(model_options$outlier$huber)) {
    model$nls_model <- refine_with_huber(
      nls_model = model$nls_model,
      function_arguments = model_inputs$function_arguments,
      data = fit_data,
      curve_function = model_inputs$curve_function,
      outlier_options = model_options$outlier
    )
  }
  if (isTRUE(model_options$outlier$rout)) {
    routed <- refine_with_rout(
      model = model,
      data = fit_data,
      model_functions = model_functions,
      start_values = start_values,
      model_options = model_options,
      use_shared_params = use_shared_params,
      has_user_start = has_user_start,
      y_var = y_var
    )

    model <- routed$model
    fit_data <- routed$data
  }

  result <- list(data = fit_data, curve = model)
  if (isTRUE(model_options$shared$return_func)) {
    result$func <- model_inputs$curve_function
  }

  result
}

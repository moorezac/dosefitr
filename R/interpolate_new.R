#' Fit every candidate model to one combination's data
#'
#' Each model is attempted independently; a model that cannot be fitted is
#' recorded in `failures` rather than aborting the rest.
#'
#' @param data,x_var,y_var Passed to [fit_curve()].
#' @param models Names of the models to attempt.
#' @param outlier_options Passed to [fit_curve()].
#' @return A list with `fits` (named, successful only) and `failures` (named
#'   error messages).
#' @keywords internal
fit_candidate_models <- function(
    data,
    x_var,
    y_var,
    models,
    outlier_options = NA
) {
    fits <- list()
    failures <- list()

    for (model_name in models) {
        curve_result <- tryCatch(
            suppressMessages(suppressWarnings(fit_curve(
                data = data,
                x_var = x_var,
                y_var = y_var,
                curve_options = list(model = model_name),
                outlier_options = outlier_options
            ))),
            dosefitr_fit_error = function(e) {
                # `errors` defaults to character(0), not NULL, when the
                # failure came from a data check rather than nlsLM()
                failures[[model_name]] <<- if (length(e$errors)) {
                    e$errors
                } else {
                    conditionMessage(e)
                }
                NULL
            }
        )

        if (!is.null(curve_result)) {
            fits[[model_name]] <- curve_result
        }
    }

    list(fits = fits, failures = failures)
}

calculate_bic <- function(nls_model) {
    residuals <- stats::residuals(nls_model)
    residuals <- residuals[is.finite(residuals)]

    n <- length(residuals)
    k <- length(stats::coef(nls_model)) + 1L

    rss <- sum(residuals^2)

    n * log(rss / n) + k * log(n)
}

#' Pick the candidate with the lowest BIC
#'
#' @param fits Named list of fits from [fit_candidate_models()].
#' @return A list with `selected` (`NULL` if nothing could be scored),
#'   `chosen_model` and `bic_values`.
#' @keywords internal
select_best_model <- function(fits) {
    bic_values <- if (length(fits)) {
        vapply(
            fits,
            function(x) {
                tryCatch(
                    calculate_bic(x$curve$nls_model),
                    error = function(e) NA_real_
                )
            },
            numeric(1)
        )
    } else {
        NULL
    }

    chosen <- if (length(bic_values) && any(is.finite(bic_values))) {
        names(bic_values)[which.min(bic_values)]
    } else {
        NA_character_
    }

    list(
        selected = if (is.na(chosen)) NULL else fits[[chosen]],
        chosen_model = chosen,
        bic_values = bic_values
    )
}

observed_combinations <- function(
    data,
    treatment_column,
    concentration_column
) {
    combos <- unique(data[,
        c(treatment_column, concentration_column),
        drop = FALSE
    ])
    combos <- combos[
        order(combos[[treatment_column]], combos[[concentration_column]]),
        ,
        drop = FALSE
    ]

    data.frame(
        treatment = as.character(combos[[treatment_column]]),
        concentration = combos[[concentration_column]],
        stringsAsFactors = FALSE
    )
}

#' Warn or stop when no model could be fitted for a combination
#'
#' @param data The combination's data, returned with `exclude` set.
#' @param failures Per-model error messages from [fit_candidate_models()].
#' @param conc,treat Identify the combination in the message.
#' @param is_negative_control Whether this combination is the control, in
#'   which case the failure is fatal rather than a warning.
#' @return An excluded record from [new_curve_record()].
#' @keywords internal
handle_fit_failure <- function(
    data,
    failures,
    conc,
    treat,
    is_negative_control
) {
    raw <- unique(unlist(failures))
    msg <- no_fit_message(
        raw,
        prefix = sprintf(
            "No model could be fitted for '%s' at concentration %s",
            treat,
            conc
        )
    )

    # The control anchors every GR/NDR/LGR calculation
    if (is_negative_control) {
        stop_unable_to_fit(
            paste0(
                msg,
                "\n'",
                treat,
                "' is the negative control, so no metrics can be calculated."
            ),
            errors = raw
        )
    }

    warning(msg, call. = FALSE)
    data$exclude <- TRUE

    new_curve_record(data = data, excluded = TRUE)
}

#' Predict from every candidate model, labelling the selected one
#'
#' The rejected models are kept so the diagnostic plot can show them behind
#' the selected curve.
#'
#' @param all_fits Named list of fits from [fit_candidate_models()].
#' @param chosen Name of the selected model.
#' @param treatment_column,concentration_column Passed to
#'   [predict_for_group()].
#' @return One data frame of predictions, or `NULL` if none could be made.
#' @keywords internal
predict_all_candidates <- function(
    all_fits,
    chosen,
    treatment_column,
    concentration_column
) {
    predictions <- Map(
        function(f, name) {
            predict_for_group(
                nls_model = f$curve$nls_model,
                data = f$data,
                concentration_column = concentration_column,
                treatment_column = treatment_column,
                extra = list(
                    model = name,
                    label = if (identical(name, chosen)) {
                        chosen
                    } else {
                        NA_character_
                    }
                )
            )
        },
        all_fits,
        names(all_fits)
    )

    do.call(rbind, Filter(is.data.frame, predictions))
}

#' Fit and select a model for one treatment-concentration combination
#'
#' @param data The combination's prepared data.
#' @param conc,treat Identify the combination in messages.
#' @param x_var,y_var,treatment_column,concentration_column Column names.
#' @param is_negative_control Whether this combination is the control.
#' @param models Names of the candidate models.
#' @param outlier_options Passed to [fit_curve()] when refitting the winner.
#' @param loess,loess_span Fit a loess smoother instead of the candidates.
#' @return A record from [new_curve_record()].
#' @keywords internal
fit_one_combination <- function(
    data,
    conc,
    treat,
    x_var,
    y_var,
    treatment_column,
    concentration_column,
    is_negative_control,
    models,
    outlier_options,
    loess = FALSE,
    loess_span = 0.3
) {
    if ("exclude" %in% names(data) && any(data$exclude %in% TRUE)) {
        return(new_curve_record(data = data, excluded = TRUE))
    }

    if (loess) {
        chosen <- "loess"
        all_fits <- list(
            loess = list(
                curve = list(
                    nls_model = stats::loess(
                        y ~ x,
                        data = data,
                        span = loess_span
                    )
                ),
                data = data
            )
        )
    } else {
        candidates <- fit_candidate_models(
            data = data,
            x_var = x_var,
            y_var = y_var,
            models = models,
            # Candidates are compared without outlier handling; only the
            # winner is refitted with `outlier_options` below
            outlier_options = NA
        )
        selection <- select_best_model(candidates$fits)

        if (is.null(selection$selected)) {
            return(handle_fit_failure(
                data = data,
                failures = candidates$failures,
                conc = conc,
                treat = treat,
                is_negative_control = is_negative_control
            ))
        }

        chosen <- selection$chosen_model
        all_fits <- candidates$fits

        # refit the winner with outlier detection enabled
        if (!identical(outlier_options, NA)) {
            refined <- tryCatch(
                suppressMessages(suppressWarnings(fit_curve(
                    data = data,
                    x_var = x_var,
                    y_var = y_var,
                    curve_options = list(model = chosen),
                    outlier_options = outlier_options
                ))),
                dosefitr_fit_error = function(e) NULL
            )
            if (!is.null(refined)) {
                all_fits[[chosen]] <- refined
            }
        }
    }

    new_curve_record(
        curve = all_fits[[chosen]]$curve,
        data = all_fits[[chosen]]$data,
        prediction_data = predict_all_candidates(
            all_fits = all_fits,
            chosen = chosen,
            treatment_column = treatment_column,
            concentration_column = concentration_column
        )
    )
}

#' Fit Growth Curves for Every Treatment and Concentration
#' @description
#' Fits every candidate model to each treatment-concentration combination,
#' keeps the one with the lowest BIC, and optionally refits the winner with
#' outlier handling. Can alternatively fit a loess smoother to each
#' combination, and flag or exclude concentrations whose control-adjusted
#' response is non-monotonic.
#' @param data A data frame containing time-course measurements.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}. If no
#'   model can be fitted to the negative control, an error is raised.
#' @param positive_control_name Character string identifying the positive
#'   control in \code{treatment_column}. Default is \code{NA}.
#' @param loess Logical; if \code{TRUE}, fit a loess smoother instead of the
#'   candidate models.
#' @param loess_span Span passed to \code{stats::loess()} when
#'   \code{loess = TRUE}.
#' @param outlier_options Outlier-handling options (see
#'   \code{\link{fit_curve}}), applied when refitting the selected model.
#'   Default is \code{NA} (none).
#' @param models Character vector of built-in model names to compare. Default
#'   is all of \code{dosefitr_models}.
#' @param monotonic_check Logical; if \code{TRUE}, check each treatment for
#'   concentrations whose response is non-monotonic relative to the negative
#'   control, and report them.
#' @param monotonic_exclude Logical; if \code{TRUE} (and
#'   \code{monotonic_check = TRUE}), exclude the flagged combinations.
#' @return
#' A list containing:
#'  \itemize{
#'    \item `model_list`: A named list with one element per
#'      treatment-concentration combination, each containing `curve` (the
#'      selected fit, or \code{NULL} if excluded), `data` (the subset of data)
#'      and `excluded` (logical).
#'    \item `plot`: A \code{ggplot2} object comparing every candidate fit per
#'      combination, with the selected model highlighted.
#'  }
#' @details
#' Combinations for which no model can be fitted are excluded with a warning
#' and their data are marked \code{exclude = TRUE}.
#' @seealso \code{\link{fit_curve}}, \code{\link{plot_curve_concentration}},
#'   \code{\link{calc_inhibition_metrics}}
#' @export
fit_concentrations <- function(
    data,
    x_var,
    y_var,
    treatment_column,
    concentration_column,
    negative_control_name,
    positive_control_name = NA,
    loess = FALSE,
    loess_span = 0.3,
    outlier_options = NA,
    models = names(dosefitr_models)[
        names(
            dosefitr_models
        ) !=
            "beta_growth_decay"
    ],
    monotonic_check = FALSE,
    monotonic_exclude = FALSE
) {
    data <- prepare_data(data, x_var, y_var)
    data <- assign_condition(
        data = data,
        treatment_column = treatment_column,
        negative_control_name = negative_control_name,
        positive_control_name = positive_control_name
    )

    combos <- observed_combinations(
        data = data,
        treatment_column = treatment_column,
        concentration_column = concentration_column
    )

    model_list <- Map(
        function(conc, treat) {
            data_filt <- data[
                data[[concentration_column]] == conc &
                    as.character(data[[treatment_column]]) == treat,
                ,
                drop = FALSE
            ]
            fit_one_combination(
                data = data_filt,
                conc = conc,
                treat = treat,
                x_var = x_var,
                y_var = y_var,
                treatment_column = treatment_column,
                concentration_column = concentration_column,
                is_negative_control = identical(
                    treat,
                    as.character(negative_control_name)
                ),
                models = models,
                outlier_options = outlier_options,
                loess = loess,
                loess_span = loess_span
            )
        },
        combos$concentration,
        combos$treatment
    )

    names(model_list) <- paste(
        combos$treatment,
        combos$concentration,
        sep = "_"
    )

    if (monotonic_check) {
        model_list <- apply_monotonic_check(
            model_list,
            treatment_column,
            concentration_column,
            negative_control_name,
            exclude = monotonic_exclude
        )
    }

    list(
        model_list = lapply(model_list, function(x) {
            x[c("curve", "data", "excluded")]
        }),
        plot = plot_fit_comparison(
            model_list,
            data,
            treatment_column,
            concentration_column,
            x_var,
            y_var
        )
    )
}


plot_fit_comparison <- function(
    model_list,
    data,
    treatment_column,
    concentration_column,
    x_var,
    y_var
) {
    facet_key <- function(d) {
        paste(d[[treatment_column]], d[[concentration_column]], sep = "_")
    }

    data$facet <- facet_key(data)

    gg <- ggplot2::ggplot(
        mapping = ggplot2::aes(x = x, y = y, colour = model)
    ) +
        ggplot2::geom_point(
            data = data,
            mapping = ggplot2::aes(x = x, y = y),
            inherit.aes = FALSE,
            colour = "black",
            alpha = 0.25,
            show.legend = FALSE
        )

    predictions <- collect_predictions(model_list)

    if (!is.null(predictions) && nrow(predictions)) {
        predictions$facet <- facet_key(predictions)
        selected <- predictions[!is.na(predictions$label), , drop = FALSE]

        gg <- gg +
            ggplot2::geom_line(data = predictions, alpha = 0.5) +
            ggplot2::geom_line(
                data = selected,
                linewidth = 2,
                alpha = 0.75,
                show.legend = FALSE
            ) +
            ggplot2::guides(
                colour = ggplot2::guide_legend(
                    title = "Model",
                    override.aes = list(alpha = 1)
                )
            )
    }

    gg +
        ggplot2::labs(x = x_var, y = y_var) +
        theme_dosefitr() +
        ggplot2::facet_wrap(~facet, scales = "free_y")
}

#' Create a per-combination result record
#'
#' @param curve The selected fit, or `NULL` if none.
#' @param data The combination's data.
#' @param prediction_data Predictions from [predict_all_candidates()].
#' @param excluded Whether the combination was excluded.
#' @return A list with `curve`, `data`, `prediction_data` and `excluded`.
#' @keywords internal
new_curve_record <- function(
    curve = NULL,
    data,
    prediction_data = NULL,
    excluded = FALSE
) {
    list(
        curve = curve,
        data = data,
        prediction_data = prediction_data,
        excluded = excluded
    )
}

#' Find concentrations whose control-adjusted response is non-monotonic
#'
#' Predictions for the selected model of each combination are checked for
#' violations separately for each treatment, with violations scaled by the
#' fitted negative control at the same x.
#'
#' @param model_list Records from [fit_one_combination()].
#' @param treatment_column,concentration_column Column names in the data.
#' @param negative_control_name Treatment value identifying the control.
#' @return A data frame with `treatment` and `concentration` columns, one row
#'   per offending combination, or `NULL` if there are none or the check could
#'   not be run.
#' @keywords internal
find_monotonic_violations <- function(
    model_list,
    treatment_column,
    concentration_column,
    negative_control_name
) {
    selected <- collect_predictions(model_list, selected_only = TRUE)
    if (is.null(selected) || !nrow(selected)) {
        return(NULL)
    }

    treatment <- as.character(selected[[treatment_column]])
    is_control <- treatment %in% as.character(negative_control_name)

    control <- selected[is_control, c("x", "y"), drop = FALSE]
    control <- control[is.finite(control$y), , drop = FALSE]
    if (!nrow(control)) {
        warning(
            "Monotonicity check skipped: no fitted negative control to compare against.",
            call. = FALSE
        )
        return(NULL)
    }

    # One control value per x, so the merge below cannot duplicate rows
    control <- stats::aggregate(y ~ x, data = control, FUN = stats::median)
    names(control)[2] <- "control_y"

    treated <- selected[!is_control, , drop = FALSE]
    if (!nrow(treated)) {
        return(NULL)
    }

    # Each treatment is checked on its own, so drugs sharing a concentration
    # series are never pooled together
    per_treatment <- lapply(
        split(treated, treatment[!is_control]),
        function(d) {
            treat <- as.character(d[[treatment_column]][1])

            # The control is only used to scale violations. Subtracting it is
            # unnecessary: it is constant at each x, so it cancels out of the
            # comparison between concentrations.
            check <- merge(d, control, by = "x", all.x = TRUE)

            result <- tryCatch(
                check_monotonic_violations(
                    check,
                    conc = concentration_column,
                    x = "x",
                    y = "y",
                    reference = "control_y"
                ),
                error = function(e) {
                    warning(
                        sprintf(
                            "Monotonicity check skipped for '%s': %s",
                            treat,
                            conditionMessage(e)
                        ),
                        call. = FALSE
                    )
                    NULL
                }
            )

            if (is.null(result) || !any(result$flagged)) {
                return(NULL)
            }

            data.frame(
                treatment = treat,
                concentration = result$concentration[result$flagged],
                stringsAsFactors = FALSE
            )
        }
    )

    flagged <- do.call(rbind, Filter(Negate(is.null), per_treatment))
    if (is.null(flagged) || !nrow(flagged)) {
        return(NULL)
    }
    rownames(flagged) <- NULL

    flagged
}

apply_monotonic_check <- function(
    model_list,
    treatment_column,
    concentration_column,
    negative_control_name,
    exclude = FALSE
) {
    flagged <- find_monotonic_violations(
        model_list = model_list,
        treatment_column = treatment_column,
        concentration_column = concentration_column,
        negative_control_name = negative_control_name
    )

    if (is.null(flagged)) {
        return(model_list)
    }

    flag_labels <- paste(
        sprintf("%s at %s", flagged$treatment, flagged$concentration),
        collapse = ", "
    )

    if (!exclude) {
        message(
            "These treatment/concentration combination(s) show potential ",
            "non-monotonic responses. Consider reviewing the data before ",
            "deciding whether to exclude them:\n",
            flag_labels
        )
        return(model_list)
    }

    message(
        "The following treatment/concentration combination(s) have been ",
        "excluded for non-monotonic responses:\n",
        flag_labels
    )

    lapply(model_list, function(record) {
        if (is.null(record$data) || !nrow(record$data)) {
            return(record)
        }

        record_treatment <- unique(as.character(
            record$data[[treatment_column]]
        ))
        record_concentration <- unique(record$data[[concentration_column]])

        flagged_concentration <- flagged$concentration[
            flagged$treatment %in% record_treatment
        ]
        is_flagged <- any(record_concentration %in% flagged_concentration)

        record$data$exclude <- is_flagged

        if (is_flagged) {
            new_curve_record(
                data = record$data,
                excluded = TRUE
            )
        } else {
            record
        }
    })
}


collect_predictions <- function(model_list, selected_only = FALSE) {
    preds <- lapply(model_list, function(x) x$prediction_data)
    preds <- Filter(is.data.frame, preds)

    if (!length(preds)) {
        return(NULL)
    }

    out <- do.call(rbind, preds)
    rownames(out) <- NULL

    if (selected_only) {
        out <- out[!is.na(out$label), , drop = FALSE]
        rownames(out) <- NULL
    }

    out
}

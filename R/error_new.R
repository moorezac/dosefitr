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

no_fit_message <- function(
    errors,
    prefix = "Unable to fit a model to this data"
) {
    hints <- translate_fit_error(errors)

    if (!length(hints)) {
        return(paste0(
            prefix,
            ". See the 'errors' element of this condition for the underlying message(s)."
        ))
    }

    paste0(prefix, ":\n", paste(hints, collapse = "\n"))
}

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
                if (grepl(p[1], m, ignore.case = TRUE)) {
                    return(p[2])
                }
            }
            NA_character_
        },
        character(1),
        USE.NAMES = FALSE
    )

    unique(hints[!is.na(hints)])
}

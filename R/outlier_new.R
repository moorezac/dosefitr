#' Refit a model with Huber weights from iteratively reweighted least squares
#'
#' Huber weights are estimated by IRWLS starting from the unweighted fit, then
#' the model is refitted with `nlsLM()` using the final weights.
#'
#' @param nls_model The unweighted fitted model.
#' @param function_arguments Argument list from [create_function_arguments()].
#' @param data The data the model was fitted on.
#' @param curve_function The curve function referenced by the formula.
#' @param outlier_options Outlier options; `huber_iter`, `huber_k` and
#'   `huber_tol` are used. Missing entries take their defaults.
#' @return The weighted model, or `nls_model` if the refit failed.
#' @keywords internal
huber_irwls <- function(
    nls_model,
    function_arguments,
    data,
    curve_function,
    outlier_options = NULL
) {
    outlier_options <- utils::modifyList(
        dosefitr_defaults$outlier,
        if (is.list(outlier_options)) outlier_options else list()
    )

    iter_max <- outlier_options$huber_iter
    k <- outlier_options$huber_k
    tol <- outlier_options$huber_tol

    coef_curr <- stats::coef(nls_model)
    coef_names <- names(coef_curr)

    m <- nls_model$m
    env <- m$getEnv()
    grad_fun <- m$gradient
    x <- get("x", env)
    y <- get("y", env)

    if ("group" %in% names(env)) {
        group <- unique(get("group", env))
    } else {
        group <- NULL
    }

    huber_psi <- function(r, k) {
        s <- stats::median(abs(r - stats::median(r))) * 1.4826
        s <- max(s, 1e-6)
        u <- r / s
        w <- ifelse(abs(u) <= k, 1, k / abs(u))
        w[!is.finite(w)] <- 1
        w <- pmax(w, 0.01)
        w <- pmin(w, 100)

        w
    }

    damp <- 0.5
    w <- rep(1, length(y))

    for (i in seq_len(iter_max)) {
        y_hat <- numeric(length(y))
        if (!is.null(group)) {
            for (j in group) {
                idx <- which(env$group == j)
                y_hat[idx] <- do.call(
                    curve_function,
                    append(as.list(coef_curr), list(x = x[idx], group = j))
                )
            }
        } else {
            y_hat <- do.call(
                curve_function,
                append(as.list(coef_curr), list(x = x))
            )
        }

        if (any(!is.finite(y_hat))) {
            break
        }

        res <- y - y_hat
        w <- huber_psi(res, k)
        X <- grad_fun()
        Xw <- X * sqrt(w)
        rw <- res * sqrt(w)

        XtX <- crossprod(Xw)
        Xtr <- crossprod(Xw, rw)

        step <- tryCatch(
            solve(XtX + diag(ncol(XtX)), Xtr),
            error = function(e) qr.solve(XtX + diag(ncol(XtX)), Xtr)
        )

        coef_new <- as.vector(coef_curr + (damp * step))
        names(coef_new) <- coef_names

        if (sqrt(sum((coef_new - coef_curr)^2)) < tol) {
            break
        }

        coef_curr <- coef_new
    }

    weighted_function_arguments <- function_arguments
    weighted_function_arguments$weights <- w

    new_fit <- fit_with_nls(
        function_arguments = weighted_function_arguments,
        data = data,
        curve_function = curve_function
    )

    if (!new_fit$ok) {
        warning(
            "Huber refinement failed; using the unweighted fit.",
            call. = FALSE
        )
        return(nls_model)
    }

    new_fit$nls_model
}


#' Identify outliers with the ROUT method
#'
#' Residuals are scaled (by MAD, or by the robust standard deviation of the
#' residuals as in Motulsky & Brown, 2006) and tested with a false discovery
#' rate of `rout_q`.
#'
#' @param nls_model The fitted model.
#' @param outlier_options Outlier options; `rout_q` and `rout_scale` are used.
#'   Missing entries take their defaults.
#' @return Row indices of the outliers, or `NULL` if there are none.
#' @keywords internal
find_rout_indices <- function(
    nls_model,
    outlier_options = NULL
) {
    outlier_options <- utils::modifyList(
        dosefitr_defaults$outlier,
        if (is.list(outlier_options)) outlier_options else list()
    )

    q <- outlier_options$rout_q
    scale_method <- outlier_options$rout_scale

    residuals_vec <- stats::residuals(nls_model)
    abs_residuals <- abs(residuals_vec)
    n <- length(residuals_vec)
    df_resid <- n - length(stats::coef(nls_model))

    scale_factor <- if (scale_method == "mad") {
        stats::mad(residuals_vec, center = 0)
    } else if (scale_method == "quantile") {
        unname(stats::quantile(
            abs_residuals,
            probs = stats::pnorm(1) - stats::pnorm(-1)
        )) *
            (n / df_resid)
    } else {
        stop("`rout_scale` must be either 'mad' or 'quantile'")
    }

    if (!is.finite(scale_factor) || scale_factor <= 0) {
        warning(
            "Could not estimate the residual scale; no outliers flagged.",
            call. = FALSE
        )
        return(NULL)
    }

    sorted <- sort(abs_residuals, index.return = TRUE)

    alphas <- q * seq(from = n, to = 1, by = -1) / n
    p_values <- 2 *
        stats::pt(
            q = sorted$x / scale_factor,
            df = df_resid,
            lower.tail = FALSE
        )

    flagged <- which(p_values < alphas)
    if (!length(flagged)) {
        return(NULL)
    }

    sorted$ix[seq.int(min(flagged), n)]
}

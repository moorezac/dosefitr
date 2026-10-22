#' Excess of each element over the running minimum of all lower-indexed elements
#'
#' Given a vector ordered by ascending concentration, returns how far each value
#' sits *above* the best (lowest) value achieved by any weaker dose, scaled by
#' `scale`. Under a perfectly monotone-decreasing relationship every element is
#' <= 0. Positive values are monotonicity violations.
#'
#' Using the running minimum rather than the immediately preceding dose matters:
#' it means a violation cannot be masked by an adjacent dose that is itself
#' already compromised.
#'
#' @param x Numeric vector, ordered by ascending concentration.
#' @param scale Positive scalar used to normalise the excess.
#' @return Numeric vector the same length as `x`. First element is 0 by
#'   definition (nothing to compare it against).
#' @keywords internal
monotone_excess <- function(x, scale = 1) {
    n <- length(x)
    if (n < 2L) {
        return(rep(0, n))
    }
    ceiling_lower <- c(Inf, cummin(x)[-n])
    out <- (x - ceiling_lower) / scale
    out[1] <- 0
    out
}


#' Flag doses that violate the monotone dose-response assumption
#'
#' For each x value (typically time), responses are compared across
#' concentrations in ascending order, and a concentration is counted as
#' violating when its response exceeds the running minimum of all lower
#' concentrations by more than `tol` times the reference response at that x
#' (see [monotone_excess()]). A concentration is flagged when it violates at a
#' sufficient fraction of x values.
#'
#' @param data Long-format data frame with one row per observation.
#' @param conc,x,y Column names, as strings, for concentration, the
#'   independent variable (typically time) and the response.
#' @param reference Column name, as a string, of the reference response used
#'   for scaling (typically the negative control at the same x). It should be
#'   constant within each x value; replicates are collapsed with `collapse`.
#' @param tol Excess above which a concentration counts as violating at a
#'   given x value, as a fraction of the absolute reference response at that
#'   x. The default, 0.05, is 5% of the reference.
#' @param min_frac The fraction of x values that must violate before a
#'   concentration is flagged. Real artefacts persist; noise does not.
#' @param collapse Function used to collapse replicates at each
#'   (concentration, x). Defaults to `median`.
#'
#' @return A data frame, one row per concentration, with:
#'   \item{concentration}{the concentration}
#'   \item{statistic}{reserved; currently always `NA`}
#'   \item{max_excess}{largest normalised violation observed}
#'   \item{frac_violating}{fraction of x values violating}
#'   \item{flagged}{whether this concentration violated the assumption}
#'
#' @details Requires at least 3 concentrations and 2 x values, otherwise an
#'   error is raised. x values where the reference is zero or missing are
#'   ignored. Tests for a monotone *decreasing* relationship; for a
#'   relationship expected to be monotone *increasing*, negate `y` before
#'   calling.
#'
#'   Scaling by the reference, rather than by the response at the lowest
#'   concentration, keeps the denominator away from zero: when responses are
#'   expressed relative to a control, the lowest dose is often close to zero.
#' @keywords internal
check_monotonic_violations <- function(
    data,
    conc = "concentration",
    x = "x",
    y = "y",
    reference = "reference",
    tol = 0.05,
    min_frac = 0.5,
    collapse = stats::median
) {
    stopifnot(all(c(conc, x, y, reference) %in% names(data)))

    d <- data.frame(
        co = as.numeric(data[[conc]]),
        x = as.numeric(data[[x]]),
        y = as.numeric(data[[y]]),
        ref = as.numeric(data[[reference]])
    )

    d <- d[stats::complete.cases(d), , drop = FALSE]

    if (!nrow(d)) {
        stop("no complete observations to check")
    }

    # Collapse replicates at concentration x time
    d <- stats::aggregate(
        cbind(y, ref) ~ co + x,
        data = d,
        FUN = collapse
    )

    # Concentrations and time points
    concs <- sort(unique(d$co))
    xs <- sort(unique(d$x))

    if (length(concs) < 3L) {
        stop("need at least 3 concentrations")
    }

    if (length(xs) < 2L) {
        stop("need at least 2 x/time points")
    }

    m <- matrix(
        NA_real_,
        nrow = length(xs),
        ncol = length(concs),
        dimnames = list(
            as.character(xs),
            as.character(concs)
        )
    )

    m[cbind(
        match(d$x, xs),
        match(d$co, concs)
    )] <- d$y

    # One reference value per x
    reference_by_x <- vapply(
        xs,
        function(value) collapse(d$ref[d$x == value]),
        numeric(1)
    )

    ex <- matrix(
        NA_real_,
        nrow = nrow(m),
        ncol = ncol(m)
    )

    for (i in seq_len(nrow(m))) {
        ex[i, ] <- monotone_excess(
            m[i, ],
            scale = abs(reference_by_x[i])
        )
    }

    # A zero or missing reference gives Inf/NaN; those x values are ignored
    ex[!is.finite(ex)] <- NA_real_

    frac <- colMeans(
        ex > tol,
        na.rm = TRUE
    )

    frac[is.nan(frac)] <- 0

    mx <- apply(
        ex,
        2,
        function(z) {
            if (all(is.na(z))) {
                NA_real_
            } else {
                max(z, na.rm = TRUE)
            }
        }
    )

    out <- data.frame(
        concentration = concs,
        statistic = NA_real_,
        max_excess = as.numeric(mx),
        frac_violating = as.numeric(frac),
        flagged = frac >= min_frac
    )

    out$flagged[is.na(out$flagged)] <- FALSE
    rownames(out) <- NULL

    return(out)
}

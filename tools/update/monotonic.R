#' Excess of each element over the running minimum of all lower-indexed elements
#' Given a vector ordered by ascending concentration, returns how far each value
#' sits *above* the best (lowest) value achieved by any weaker dose, scaled by
#' `scale`. Under a perfectly monotone-decreasing relationship every element is
#' <= 0. Positive values are monotonicity violations.
#' Using the running minimum rather than the immediately preceding dose matters:
#' it means a violation cannot be masked by an adjacent dose that is itself
#' already compromised.
#' @param x Numeric vector, ordered by ascending concentration.
#' @param scale Positive scalar used to normalise the excess.
#' @return Numeric vector the same length as `x`. First element is 0 by
#'   definition (nothing to compare it against).
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
#' @param data Long-format data frame with one row per (well, concentration, time).
#' @param conc,time,response Column names, as strings.
#' @param metric Which quantity to test for monotonicity.
#'   \describe{
#'     \item{"growth"}{Slope of log(response) vs time, one value per
#'       concentration. Collapses the whole trajectory into a single robust
#'       number, so it is insensitive to frame-to-frame jitter. Usually the
#'       better default.}
#'     \item{"area"}{Raw response, tested independently at every timepoint. More
#'       sensitive, but needs `min_frac` to avoid firing on single noisy frames.}
#'   }
#' @param tol Excess above which a dose counts as violating. For `metric =
#'   "area"` this is a fraction of the lowest-dose area at that timepoint
#'   (default 0.05, i.e. 5% of vehicle). For `metric = "growth"` it is a fraction
#'   of the full range of growth rates across the plate (default 0.10). Pass
#'   `NULL` to use the metric-appropriate default.
#' @param min_frac For `metric = "area"` only: the fraction of timepoints that
#'   must violate before a dose is flagged. Real death artefacts persist; noise
#'   does not.
#' @param collapse Function used to collapse replicate wells at each
#'   (concentration, time). Defaults to `median`.
#' @param conc TODO: description.
#' @param x TODO: description.
#' @param y TODO: description.
#' @return A data frame, one row per concentration, with:
#'   \item{statistic}{the tested quantity (growth slope, or NA for "area")}
#'   \item{max_excess}{largest normalised violation observed}
#'   \item{frac_violating}{fraction of timepoints violating ("area" only)}
#'   \item{flagged}{did this dose violate the assumption}
#'   \item{excluded}{is it part of the contiguous flagged block at the top of the
#'     dose range -- this is the column to filter on}
#' @details `excluded` is deliberately stricter than `flagged`. Loss of the
#'   spheroid is a top-of-range phenomenon: it should affect the highest doses
#'   and everything above them. Requiring the flagged set to be contiguous with
#'   the top dose prevents a single noisy mid-range dose from being discarded.
#'   For a relationship expected to be monotone *increasing*, negate `response`
#'   before calling.
check_monotonic_violations <- function(
    data,
    conc = "concentration",
    x = "x",
    y = "y",
    tol = 0.05,
    min_frac = 0.5,
    collapse = stats::median
) {
    stopifnot(all(c(conc, x, y) %in% names(data)))

    d <- data.frame(
        co = as.numeric(data[[conc]]),
        x = as.numeric(data[[x]]),
        y = as.numeric(data[[y]])
    )

    d <- d[stats::complete.cases(d), , drop = FALSE]

    # Collapse replicates at concentration x time
    d <- stats::aggregate(
        y ~ co + x,
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

    ex <- matrix(
        NA_real_,
        nrow = nrow(m),
        ncol = ncol(m)
    )

    for (i in seq_len(nrow(m))) {
        ex[i, ] <- monotone_excess(
            m[i, ],
            scale = abs(m[i, 1L])
        )
    }

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

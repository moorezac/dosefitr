dosefitr_models <- list(
    exponential_growth = list(
        curve_function = function(x, y0, k) {
            y0 * exp(k * x)
        },
        start_function = function(x, y) {
            list(y0 = min(y), k = 0)
        }
    ),
    exponential_plateau = list(
        curve_function = function(x, y0, ym, k) {
            ym - (ym - y0) * exp(-k * x)
        },
        start_function = function(x, y) {
            list(ym = max(y), y0 = min(y), k = 0)
        }
    ),
    four_param_sigmoid = list(
        curve_function = function(x, top, bottom, ec50, slope) {
            bottom + (top - bottom) / (1 + (ec50 / x)^slope)
        },
        start_function = function(x, y) {
            x_at_y_mid <- function(x, y) {
                if (length(x) != length(y)) {
                    stop("x and y must have equal lengths")
                }
                mid <- min(y) + (max(y) - min(y)) / 2
                index <- which.min(abs(y - mid))
                return(x[index])
            }
            list(
                bottom = min(y),
                top = max(y),
                ec50 = x_at_y_mid(x, y),
                slope = sign(
                    mean(y[x %in% max(x)]) -
                        mean(y[
                            x %in%
                                min(x)
                        ])
                )
            )
        }
    ),
    four_param_sigmoid_log = list(
        curve_function = function(x, top, bottom, log_ec50, slope) {
            bottom + (top - bottom) / (1 + 10^((log_ec50 - x) * slope))
        },
        start_function = function(x, y) {
            x_at_y_mid <- function(x, y) {
                if (length(x) != length(y)) {
                    stop("x and y must have equal lengths")
                }
                mid <- min(y) + (max(y) - min(y)) / 2
                index <- which.min(abs(y - mid))
                return(x[index])
            }
            list(
                bottom = min(y),
                top = max(y),
                log_ec50 = x_at_y_mid(x, y),
                slope = sign(y[which.max(x)] - y[which.min(x)])
            )
        }
    ),
    biphasic_four_param_sigmoid = list(
        curve_function = function(
            x,
            top,
            bottom,
            frac,
            ec50_1,
            nh1,
            ec50_2,
            nh2
        ) {
            bottom +
                (top - bottom) * frac / (1 + (ec50_1 / x)^nh1) +
                (top -
                    bottom) *
                    (1 - frac) /
                    (1 + (ec50_2 / x)^nh2)
        },
        start_function = function(x, y) {
            list(
                bottom = min(y),
                top = max(y),
                frac = 0.5,
                ec50_1 = min(x) +
                    0.3 * (max(x) - min(x)),
                nh1 = 1,
                ec50_2 = min(x) +
                    0.7 *
                        (max(x) - min(x)),
                nh2 = 1
            )
        }
    ),
    five_param_sigmoid = list(
        curve_function = function(x, top, bottom, s, ec50, slope) {
            bottom +
                (top - bottom) / (1 + (2^(1 / s) - 1) * ((ec50 / x)^slope))^s
        },
        start_function = function(x, y) {
            x_at_y_mid <- function(x, y) {
                if (length(x) != length(y)) {
                    stop("x and y must have equal lengths")
                }
                mid <- min(y) + (max(y) - min(y)) / 2
                index <- which.min(abs(y - mid))
                return(x[index])
            }
            list(
                bottom = min(y),
                top = max(y),
                ec50 = x_at_y_mid(x, y),
                slope = sign(
                    mean(y[x %in% max(x)]) -
                        mean(y[
                            x %in%
                                min(x)
                        ])
                ),
                s = 0.5
            )
        }
    ),
    five_param_sigmoid_log = list(
        curve_function = function(x, top, bottom, s, log_ec50, slope) {
            bottom +
                (top - bottom) /
                    ((1 +
                        10^(((log_ec50 +
                            (1 / slope) *
                                log10((2^(1 / s)) - 1)) -
                            x) *
                            slope))^s)
        },
        start_function = function(x, y) {
            x_at_y_mid <- function(x, y) {
                if (length(x) != length(y)) {
                    stop("x and y must have equal lengths")
                }
                mid <- min(y) + (max(y) - min(y)) / 2
                index <- which.min(abs(y - mid))
                return(x[index])
            }
            list(
                bottom = min(y),
                top = max(y),
                log_ec50 = x_at_y_mid(x, y),
                slope = sign(
                    mean(y[x %in% max(x)]) -
                        mean(y[
                            x %in%
                                min(x)
                        ])
                ),
                s = 0.5
            )
        }
    ),
    logistic_growth = list(
        curve_function = function(x, ym, y0, k) {
            ym * y0 / ((ym - y0) * exp(-k * x) + y0)
        },
        start_function = function(x, y) {
            list(ym = max(y), y0 = min(y), k = 0)
        }
    ),
    gompertz_growth = list(
        curve_function = function(x, ym, y0, k) {
            ym * (y0 / ym)^(exp(-k * x))
        },
        start_function = function(x, y) {
            list(ym = max(y), y0 = min(y), k = 0)
        }
    ),
    beta_growth_decay = list(
        curve_function = function(x, ym, te, tm) {
            ym * (1 + (te - x) / (te - tm)) * (x / te)^(te / (te - tm))
        },
        start_function = function(x, y) {
            # te (end of growth) must differ from tm (time of maximum growth
            # rate), otherwise te - tm is zero and the curve is undefined
            te <- x[which.max(y)]
            if (!is.finite(te) || te <= 0) {
                te <- max(x)
            }
            list(
                ym = max(y),
                te = te,
                tm = te / 2
            )
        }
    )
)

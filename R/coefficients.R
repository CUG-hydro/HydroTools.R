#' skewness and kurtosis
#'
#' Inherit from package `e1071`
#' @param x a numeric vector containing the values whose skewness is to be
#' computed.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @param type an integer between 1 and 3 selecting one of the algorithms for
#' computing skewness.
#'
#' @examples
#' x <- rnorm(100)
#' coef_kurtosis <- kurtosis(x)
#' coef_skewness <- skewness(x)
#' @keywords internal
#' @export
kurtosis <- function(x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x^4) / (sum(x^2)^2)

    y <- if (type == 1) {
        r - 3
    } else if (type == 2) {
        if (n < 4) {
            warning("Need at least 4 complete observations.")
            return(NA_real_)
        }
        ((n + 1) * (r - 3) + 6) * (n - 1) / ((n - 2) * (n - 3))
    } else {
        r * (1 - 1 / n)^2 - 3
    }
    y
}

#' @rdname kurtosis
#' @export
skewness <- function(x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3) / (sum(x^2)^(3 / 2))

    if (type == 2) {
        if (n < 3) {
            warning("Need at least 3 complete observations.")
            return(NA_real_)
        }
        y <- y * sqrt(n * (n - 1)) / (n - 2)
    } else if (type == 3) {
        y <- y * ((1 - 1 / n))^(3 / 2)
    }
    y
}

#' weighted CV
#' @param x Numeric vector
#' @param w weights of different point
#'
#' @return Named numeric vector, (mean, sd, cv).
#' @examples
#' x <- rnorm(100)
#' coefs <- cv_coef(x)
#' @keywords internal
#' @export
cv_coef <- function(x, w) {
    if (missing(w)) w <- rep(1, length(x))
    if (length(x) == 0) {
        return(c(mean = NA_real_, sd = NA_real_, cv = NA_real_))
    }
    # rm NA_real_
    I <- is.finite(x)
    x <- x[I]
    w <- w[I]

    mean <- sum(x * w) / sum(w)
    sd <- sqrt(sum((x - mean)^2 * w) / sum(w))
    cv <- sd / mean
    c(mean = mean, sd = sd, cv = cv) # quickly return
}


#' Critical value of determined correlation
#'
#' @param n length of observation.
#' @param NumberOfPredictor Number of predictor, including constant.
#' @param alpha significant level.
#'
#' @return `F` statistic and `R2` at significant level.
#'
#' @keywords internal
#' @references
#' Chen Yanguang (2012), Geographical Data analysis with MATLAB.
#' @examples
#' R2_critical <- R2_sign(30, NumberOfPredictor = 2, alpha = 0.05)
#' @export
R2_sign <- function(n, NumberOfPredictor = 2, alpha = 0.05) {
    freedom_r <- NumberOfPredictor - 1 # regression
    freedom_e <- n - NumberOfPredictor # error

    F <- qf(1 - alpha, freedom_r, freedom_e)
    R2 <- 1 - 1 / (1 + F * freedom_r / freedom_e)

    # F = 485.1
    # F = R2/freedom_r/((1-R2)/freedom_e)
    # Rc = sqrt(/(qf(1 - alpha, 1, freedom) + freedom)) %TRUE>% print  # 0.11215
    return(list(F = F, R2 = R2))
}

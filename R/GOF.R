#' @param ... ignored
#' 
#' @rdname GOF
#' @export
NSE <- function(yobs, ysim, w, ...) {
    if (missing(w)) w <- rep(1, length(yobs))

    ind <- valindex(yobs, ysim)
    w <- w[ind]

    y_mean <- sum(yobs[ind] * w) / sum(w)
    # R2: the portion of regression explained variance, also known as
    # coefficient of determination

    # SSR <- sum((ysim - y_mean)^2 * w)
    SST <- sum((yobs[ind] - y_mean)^2 * w)
    # R2     <- SSR / SST

    RE <- ysim[ind] - yobs[ind]
    # Bias <- sum(w * RE) / sum(w) # bias
    # Bias_perc <- Bias / y_mean # bias percentage
    # MAE <- sum(w * abs(RE)) / sum(w) # mean absolute error
    RMSE <- sqrt(sum(w * (RE)^2) / sum(w)) # root mean sqrt error

    NSE <- 1 - sum((RE)^2 * w) / SST # NSE coefficient
    NSE
}

#' @rdname GOF
#' @export
KGE <- function(yobs, ysim, ...) {
  ind <- valindex(yobs, ysim)
  yobs <- yobs[ind]
  ysim <- ysim[ind]
  
  c1 = cor(yobs, ysim)
  c2 = sd(ysim) / sd(yobs)
  c3 = mean(ysim) / mean(yobs)  
  
  1 - sqrt((c1 - 1)^2 + (c2 - 1)^2 + (c3 - 1)^2)
}

#' GOF
#'
#' Good of fitting
#'
#' @param yobs Numeric vector, observations
#' @param ysim Numeric vector, corresponding simulated values
#' @param w Numeric vector, weights of every points. If w included, when
#' calculating mean, Bias, MAE, RMSE and NSE, w will be taken into considered.
#' @param include.cv If true, cv will be included.
#' @param include.r If true, r and R2 will be included.
#' 
#' @return
#' * `RMSE` root mean square error
#' * `NSE` NASH coefficient
#' * `MAE` mean absolute error
#' * `AI` Agreement index (only good points (w == 1)) participate to
#' calculate. See details in Zhang et al., (2015).
#' * `Bias` bias
#' * `Bias_perc` bias percentage
#' * `n_sim` number of valid obs
#' * `cv` Coefficient of variation
#' * `R2` correlation of determination
#' * `R` pearson correlation
#' * `pvalue` pvalue of `R`
#'
#' @references
#' 1. https://en.wikipedia.org/wiki/Coefficient_of_determination
#' 2. https://en.wikipedia.org/wiki/Explained_sum_of_squares
#' 3. https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
#' 4. Zhang Xiaoyang (2015), http://dx.doi.org/10.1016/j.rse.2014.10.012
#'
#' @examples
#' yobs = rnorm(100)
#' ysim = yobs + rnorm(100)/4
#' GOF(yobs, ysim)
#' 
#' @importFrom dplyr tibble
#' @export
GOF <- function(yobs, ysim, w, include.cv = FALSE, include.r = TRUE){
    if (missing(w)) w <- rep(1, length(yobs))

    # remove NA_real_ and Inf values in ysim, yobs and w
    valid <- function(x) !is.na(x) & is.finite(x)

    I <- which(valid(ysim) & valid(yobs) & valid(w))
    # n_obs <- length(yobs)
    n_sim <- length(I)

    ysim <- ysim[I]
    yobs <- yobs[I]
    w     <- w[I]

    R      <- NA_real_
    pvalue <- NA_real_
    
    if (include.cv) {
        CV_obs <- cv_coef(yobs, w)
        CV_sim <- cv_coef(ysim, w)
    }
    if (is_empty(yobs)){
        out <- tibble(
          R, pvalue, R2 = NA_real_,
          NSE = NA_real_, KGE = NA_real_, RMSE = NA_real_, MAE = NA_real_,
          Bias = NA_real_, Bias_perc = NA_real_, AI = NA_real_, n_sim = NA_integer_
        )
        if (include.cv) out <- cbind(out, CV_obs, CV_sim)
        return(out)
    }

    # R2: the portion of regression explained variance, also known as
    # coefficient of determination
    KGE = KGE(ysim, yobs)
    y_mean <- sum(yobs * w) / sum(w)

    SSR    <- sum( (ysim - y_mean)^2 * w)
    SST    <- sum( (yobs - y_mean)^2 * w)
    # R2     <- SSR / SST

    RE     <- ysim - yobs
    Bias   <- sum ( w*RE)     /sum(w)                     # bias
    Bias_perc <- Bias/y_mean                              # bias percentage
    MAE    <- sum ( w*abs(RE))/sum(w)                     # mean absolute error
    RMSE   <- sqrt( sum(w*(RE)^2)/sum(w) )                # root mean sqrt error

    # https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
    NSE  <- 1  - sum( (RE)^2 * w) / SST # NSE coefficient

    # Observations number are not same, so comparing correlation coefficient
    # was meaningless.
    # In the current, I have no idea how to add weights `R`.
    if (include.r){
        
        tryCatch({
            cor.obj <- cor.test(yobs, ysim, use = "complete.obs")
            R       <- cor.obj$estimate[[1]]
            pvalue  <- cor.obj$p.value
        }, error = function(e){
            message(e$message)
        })
        R2 = R^2
    }
    # In Linear regression, R2 = R^2 (R is pearson cor)
    # R2     <- summary(lm(ysim ~ yobs))$r.squared # low efficient

    # AI: Agreement Index (only good values(w==1) calculate AI)
    AI <- NA_real_
    I2 <- which(w == 1)
    if (length(I2) >= 2) {
        yobs = yobs[I2]
        ysim = ysim[I2]
        y_mean = mean(yobs)
        AI = 1 - sum( (ysim - yobs)^2 ) / sum( (abs(ysim - y_mean) + abs(yobs - y_mean))^2 )
    }

    out <- tibble(R, pvalue, R2, NSE, KGE, RMSE, MAE, 
             Bias, Bias_perc, AI = AI, n_sim = n_sim)
    if (include.cv) out <- cbind(out, CV_obs, CV_sim)
    return(out)
}

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

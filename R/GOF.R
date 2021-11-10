#' @importFrom hydroGOF KGE
#' @export
hydroGOF::KGE

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
#' @importFrom hydroGOF KGE
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

    if (include.cv) {
        CV_obs <- cv_coef(yobs, w)
        CV_sim <- cv_coef(ysim, w)
    }
    if (is_empty(yobs)){
        out <- c(RMSE = NA_real_, 
            KGE = NA_real_,
            NSE = NA_real_, MAE = NA_real_, AI = NA_real_,
            Bias = NA_real_, Bias_perc = NA_real_, n_sim = NA_real_)

        if (include.r) out <- c(out, R2 = NA_real_, R = NA_real_, pvalue = NA_real_)
        if (include.cv) out <- c(out, obs = CV_obs, sim = CV_sim)
        return(out)
    }

    # R2: the portion of regression explained variance, also known as
    # coefficient of determination
    KGE = KGE(ysim, yobs)
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    # https://en.wikipedia.org/wiki/Explained_sum_of_squares
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
        R      <- NA_real_
        pvalue <- NA_real_
        
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

    out <- c(RMSE = RMSE, KGE = KGE, NSE = NSE, MAE = MAE, AI = AI,
             Bias = Bias, Bias_perc = Bias_perc, n_sim = n_sim)

    if (include.r) out <- c(out, R2 = R2, R = R, pvalue = pvalue)
    if (include.cv) out <- c(out, obs = CV_obs, sim = CV_sim)
    return(out)
}

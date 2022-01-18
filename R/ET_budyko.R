#' Budyko Curve
#' 
#' @param PET.P PET divide P, 
#' - `PET`: potential evapotranspiration, can be calculated by [ET0_Penman48()].
#' - `P`: Precipitation
#' - `AET.P`: `AET / P`
#' - `PET.P`: `PET / P`
#' @param par scalar value, the parameter of Budyko curve
#' @param method one of `c("budyko", "Fu1981", "Zhang2001", "Yang2008", "Pike1964")`, 
#' or `all`.
#' 
#' @references 
#' 1. Budyko, M. I., Climate and Life, 508 pp., Academic, San Diego, Calif., 1974.
#' 2. Fu, Baw-Puh., 1981. On the Calculation of the Evaporation from Land
#'    Surface. Chinese Journal of Atmospheric Sciences, 5(1), 23-31.
#' 3. Zhang, L., Dawes, W. R., & Walker, G. R. (2001). Response of mean annual
#'    evapotranspiration to vegetation changes at catchment scale. Water
#'    Resources Research, 37(3), 701–708. <doi:10.1029/2000WR900325>
#' 4. Yang, H., Yang, D., Lei, Z., & Sun, F. (2008). New analytical derivation
#'    of the mean annual water-energy balance equation. Water Resources
#'    Research. <doi:10.1029/2007WR006135>
#' 
#' @export
budyko <- function(PET.P, par = 2, method = "Fu1981") {
    # PET.P = PET/P
    # AET.P = ET/P
    w = par
    n = par
    p = par
    if (method == "budyko") {
        # ET = ( (1 - exp(-PET.P)) * P * PET * tanh(1 / PET.P) ) ^ (1/w) # budyko, w = 2
        AET.P = ((1 - exp(-PET.P)) * PET.P * tanh(1 / PET.P)) ^ (1/w)
    } else if (method == "Fu1981") {
        # 傅抱璞 1981
        AET.P = 1 + PET.P - (1 + (PET.P)^w)^(1 / w)                   # Fu1981, Eq. 21
    } else if (method == "Zhang2001") {
        AET.P = (1 + w * PET.P) / ( 1 + w * PET.P + PET.P^-1 )      # Zhang2001, Eq. 6    
    } else if (method == "Yang2008") {
        # ET    = PET * P / (P^n + PET^n) ^ (1/n) # Yang2008, Eq. 25    
        # AET.P = PET / (P^n + PET^n) ^ (1/n)
        AET.P = 1 / (1 + (1 / PET.P)^n)^(1 / n)
    } else if (method == "Pike1964") {
        # p = 2, default
        AET.P = (1 + (PET.P)^(-p))^(-1 / p) # Pike 1964, also known as "turc-pike"
    }
}

#' @rdname budyko
#' @export
budyko_goal <- function(AET.P, PET.P, par = 2, method = "Fu1981", ...) {
    .methods = c("budyko", "Fu1981", "Zhang2001", "Yang2008", "Pike1964") 
    methods <- if (method == "all") .methods else method
    
    map(methods, function(method) {
        ysim = budyko(PET.P, par, method = method)
        GOF(AET.P, ysim) %>% as.list() %>% as.data.table() %>% cbind(method, .)
    }) %>% do.call(rbind, .)    
}

#' @param ... ignored
#' 
#' @rdname budyko
#' @importFrom purrr map
#' @export
budyko_fit <- function(AET.P, PET.P, par = 2, method = "Fu1981", ...) {
    methods = c("budyko", "Fu1981", "Zhang2001", "Yang2008", "Pike1964")

    if (method == "all") {
        r = map(methods, ~ budyko_fit(AET.P, PET.P, par = par, method = .x))
        map(r, "gof") %>% do.call(rbind, .) %>% as.data.table() %>% 
            cbind(method = methods, .)
    } else {
        data = data.frame(AET.P, PET.P)
        l = nls(formula = AET.P ~ budyko(PET.P, par, method = method), data, start = c(par = 2))
        s = summary(l)
        par = s$coefficients[, 1]
        ysim = predict(l)
        gof = GOF(AET.P, ysim)
        listk(method, par, gof)
    }
}

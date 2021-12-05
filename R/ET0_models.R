#' @name ET0_models
#' @title Potential Evapotranspiration models
#' 
#' @description
#' - `ET0_PM93`: Penman equation simplified by Shuttleworth 1993
#' - `ET0_PM98`: Penman-Monteith reference crop evapotranspiration, FAO56
#' 
#' @references
#' 1. https://en.wikipedia.org/w/index.php?title=Penman_equation
#' 2. Allen, R. G., & Luis S. Pereira. (1998). Crop
#'    evapotranspiration-Guidelines for computing crop water requirements-FAO
#'    Irrigation and drainage paper 56. European Journal of Agronomy, 34(3),
#'    144–152. doi:10.1016/j.eja.2010.12.001
NULL

#' @param Rn net radiation (W m-2)
#' @param Tair 2m air temperature (degC)
#' @param D vapor pressure deficit (kPa)
#' @param Pa surface air pressure (kPa)
#' @param wind wind speed at the height of `z.wind`
#' @inheritParams cal_U2
#' 
#' @rdname ET0_models
#' @export
ET0_eq <- function(Rn, Tair, D, Pa = atm, ...) {
    
    lambda <- cal_lambda(Tair) # MJ kg-1
    slope <- cal_slope(Tair) # kPa degC-1
    gamma <- Cp * Pa / (epsilon * lambda) # kPa degC-1
    
    coef_W2mm <- 0.086400 / lambda
    Eeq <- slope / (slope + gamma) * Rn * coef_W2mm

    data.table(lambda, slope, gamma, Eeq)
}

#' @param Rn land surface net radiation, W m-2
#' 
#' @importFrom dplyr mutate
#' @rdname ET0_models
#' @export
ET0_PM93 <- function(Rn, Tair, D, Pa = atm, 
    wind, z.wind = 10) 
{
    dat = ET0_eq(Rn, Tair, D, Pa)
    U2 = cal_U2(wind, z.wind)

    mutate(dat, 
        Evp = gamma / (slope + gamma) * 6.43 * (1 + 0.536 * U2) * D / lambda, 
        ET0_PM93 = Evp + Eeq)
}

#' @rdname ET0_models
#' @export
ET0_PM98 <- function(Rn, Tair, D, Pa = atm, 
    wind, z.wind = 10, tall_crop = FALSE)
{
    if (tall_crop) {
        p1 <- 1600
        p2 <- 0.38
    } else {
        p1 <- 900
        p2 <- 0.34
    }

    dat = ET0_eq(Rn, Tair, D, Pa)
    U2 = cal_U2(wind, z.wind)
    coef_W2mm <- 0.086400 / dat$lambda

    dat %>% mutate(
        Eeq_pm = slope / (slope + (gamma * (1 + p2 * U2))) * Rn * coef_W2mm,
        Evp_pm = gamma * p1 / (Tair + 273.15) * U2 * D / (slope + (gamma * (1 + p2 * U2))), 
        ET0_PM98 = Eeq_pm + Evp_pm
    )
}

#' @rdname ET0_models
#' @export
ET0_PT72 <- function(Rn, Tair, D, Pa = atm, alpha = 1.26, ...) {
    dat = ET0_eq(Rn, Tair, D, Pa)
    mutate(dat, ET0_PM72 = Eeq * alpha)
}

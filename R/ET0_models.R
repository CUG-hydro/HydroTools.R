#' @name ET0_models
#' @title Potential Evapotranspiration models
#'
#' @description
#' - `ET0_eq`   : Equilibrium evaporation, `slope / (slope + gamma) * Rn`
#' - `ET0_Penman48` : Penman 1948 equation simplified by Shuttleworth 1993
#' - `ET0_Monteith65`: Penman-Monteith 1965
#' - `ET0_PT72` : Priestley Taylor 1972, `ET0_PT72 = ET0_eq * 1.26`
#' - `ET0_FAO98` : Penman-Monteith reference crop evapotranspiration, FAO56
#' 
#' @example R/example/ex-ET0.R
#'
#' @references
#' 1. Penman equation. Wikipedia 2021. 
#'    <https://en.wikipedia.org/w/index.php?title=Penman_equation>
#' 2. Allen, R. G., & Luis S. Pereira. (1998). Crop
#'    evapotranspiration-Guidelines for computing crop water requirements-FAO
#'    Irrigation and drainage paper 56. European Journal of Agronomy, 34(3), 
#'    144-152. <doi:10.1016/j.eja.2010.12.001>
#' 3. Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). FAO Irrigation
#'    and drainage paper No. 56. Rome: Food and Agriculture Organization of the
#'    United Nations, 56(97), e156.
#' 4. Chapter 2 - FAO Penman-Monteith equation,
#'    <https://www.fao.org/3/x0490E/x0490e06.htm>
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
ET0_eq <- function(Rn, Tair, Pa = atm, ...) {

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
ET0_Penman48 <- function(Rn, Tair, Pa = atm, D,
    wind, z.wind = 10)
{
    dat = ET0_eq(Rn, Tair, Pa)
    U2 = cal_U2(wind, z.wind)

    mutate(dat,
        Evp = gamma / (slope + gamma) * 6.43 * (1 + 0.536 * U2) * D / lambda,
        ET0 = Evp + Eeq)
}

#' @rdname ET0_models
#' @export
ET0_Monteith65 <- function(Rn, Tair, Pa = atm, D, wind, z.wind = 10, rs = 70, ...) {
    dat = ET0_eq(Rn, Tair, Pa)
    U2 = cal_U2(wind, z.wind)
    rH = cal_rH(U2, h = 0.12)

    coef_W2mm <- 0.086400 / dat$lambda
    rou_a = 3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5, kg m-3

    dat %>% mutate(
        Eeq = slope / (slope + (gamma * (1 + rs / rH))) * Rn * coef_W2mm,
        Evp = (rou_a * Cp * D / rH) / (slope + (gamma * (1 + rs / rH))) * 86400 / lambda,
        ET0 = Eeq + Evp
    )
}

#' @rdname ET0_models
#' @export
ET0_PT72 <- function(Rn, Tair, Pa = atm, alpha = 1.26, ...) {
    dat <- ET0_eq(Rn, Tair, Pa)
    mutate(dat, ET0 = Eeq * alpha)
}

#' @rdname ET0_models
#' @export
ET0_FAO98 <- function(Rn, Tair, Pa = atm, D,
    wind, z.wind = 10, tall_crop = FALSE)
{
    if (tall_crop) {
        p1 <- 1600
        p2 <- 0.38
    } else {
        p1 <- 900
        p2 <- 0.34
    }

    dat = ET0_eq(Rn, Tair, Pa)
    U2 = cal_U2(wind, z.wind)
    coef_W2mm <- 0.086400 / dat$lambda

    dat %>% mutate(
        Eeq = slope / (slope + (gamma * (1 + p2 * U2))) * Rn * coef_W2mm,
        Evp = gamma * p1 / (Tair + 273.15) * U2 * D / (slope + (gamma * (1 + p2 * U2))),
        ET0 = Eeq + Evp
    )
}

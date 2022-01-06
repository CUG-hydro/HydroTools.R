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
#' 2. Shuttleworth, W. J. (1993) Evaporation. In: Handbook of Hydrology (ed. by
#'    D. Maidment). McGraw-Hill, New York.
#' 3. Allen, R. G., & Luis S. Pereira. (1998). Crop
#'    evapotranspiration-Guidelines for computing crop water requirements-FAO
#'    Irrigation and drainage paper 56. European Journal of Agronomy, 34(3),
#'    144-152. <doi:10.1016/j.eja.2010.12.001>
#' 4. Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). FAO Irrigation
#'    and drainage paper No. 56. Rome: Food and Agriculture Organization of the
#'    United Nations, 56(97), e156.
#' 5. Chapter 2 - FAO Penman-Monteith equation,
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
    # rou_a * Cp * dT / rH (MJ m-2 s-1)
    # rou_a ≈ 1.225 kg/m3
    # rou_a * Cp / rH = f(U2)
    # `f(U2) = 2.6 * (1 + 0.54U2)` is equivalent to Shuttleworth1993
    mutate(dat,
        Evp = gamma / (slope + gamma) * 6.43 * (1 + 0.536 * U2) * D / lambda,
        ET0 = Evp + Eeq)
}

#' @rdname ET0_models
#' @export
ET0_Monteith65 <- function(Rn, Tair, Pa = atm, D, wind, z.wind = 10, rs = 70, rH = NULL, ...) {
    dat = ET0_eq(Rn, Tair, Pa)
    U2 = cal_U2(wind, z.wind)
    if (is.null(rH)) rH = cal_rH(U2, h = 0.12)
    
    coef_W2mm <- 0.086400 / dat$lambda
    rou_a = 3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5, kg m-3
    # Cp = 1.013 * 1e-3 # MJ kg-1 degC-1
    # rou_a * Cp * dT * gH (in MJ m-2 s-1)
    # = kg m-3 * MJ kg-1 degC-1 * degC * m s-1
    # = MJ m-2 s-1
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

#' Estimate ET0 by Hargreaves equation.
#'
#' @description Estimate reference evapotranspiration (ET0) from a hypothetical
#' short grass reference surface using the the Hargreaves equation.
#'
#' @param Tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param Tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param Tavg Daily mean air temperature at 2m height `[deg Celsius]`. If not
#' provided it would estimated by averaging the Tmax and tmin.
#' @param Ra  Clear sky incoming shortwave radiation, i. e. extraterrestrial
#' radiation multiply by clear sky transmissivity (i. e. a + b, a and b are
#' coefficients of Angstrom formula. Normally 0.75) `[MJ m-2 day-1]`. If not
#' provided, must provide lat and dates.
#' @param lat Latitude `[degree]`.
#' @param dates A R Date type of a vector of Date type. If not provided, it will Regard
#' the ssd series is begin on the first day of a year.
#'
#' @return Reference evapotranspiration ET0 from a hypothetical grass reference
#' surface `[mm day-1]`.
#'
#' @export
PET_hg <- function(Tmax, Tmin, Tavg = NULL, Ra = NULL, lat = NULL, dates=NULL) {
  if(is.null(Tavg)) Tavg <- (Tmax + Tmin) / 2
  if(is.null(Ra)) Ra <- cal_Ra(lat, dates)
  0.0023 * 0.408 * Ra * sqrt(Tmax - Tmin) * (Tavg + 17.8)
}

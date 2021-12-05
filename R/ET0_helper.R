#' @name ET0_helper
#' @title helper functions for potential evapotranspiration (ET0)
#' 
#' @description
#' - `lambda`: latent heat of vaporization, about `[2.5 MJ kg-1]`.
#' - `slope`: The slope of the saturation vapour pressure curve at certain air
#'   temperature `Tair` ([deg Celsius]), `[kPa degC-1]`.
#' - `gamma`: psychrometric constant (`[kPa degC-1]`), `Cp * Pa / (epsilon * lambda)`.
#' - `U2`: 10m wind speed (m/s). According to wind profile relationship, convert U_z to U_2.
#' - `es`: saturation vapor pressure (kPa)
#' 
#' @inheritParams cal_es
#' 
#' @references 
#' 1. Allen, R. G., & Luis S. Pereira. (1998). Crop
#'    evapotranspiration-Guidelines for computing crop water requirements-FAO
#'    Irrigation and drainage paper 56. European Journal of Agronomy, 34(3),
#'    144–152. doi:10.1016/j.eja.2010.12.001
NULL

#' @param z.wind Height where measure the wind speed `[m]`. Default 10m.
#' 
#' @rdname ET0_helper
#' @export
cal_U2 <- function(Uz, z.wind = 10) {
    # log, by default natural logarithms
    if (z.wind == 2) return(Uz)
    Uz * 4.87 / log(67.8 * z.wind - 5.42)
}

#' @rdname ET0_helper
#' @export
cal_es <- function(Tair) {
    0.6108 * exp((17.27 * Tair) / (Tair + 237.3))
}

#' @param Tdew dew temperature, ([deg Celsius])
#' 
#' @rdname ET0_helper
#' @export
cal_VPD <- function(Tair, Tdew) {
    ea <- cal_es(Tdew)
    es <- cal_es(Tair)
    es - ea
}

#' @rdname ET0_helper
#' @export
cal_lambda <- function(Tair) {
    (2500 - Tair * 2.2) / 1000
}

#' @rdname ET0_helper
#' @export
cal_slope <- function(Tair) {
    4098 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) /
        (Tair + 237.3)**2
}

# #' @rdname ET0_helper
# #' @export
# delta_es = cal_slope

#' @param Pa land surface Air pressure `[kPa]`.
#' @rdname ET0_helper
#' @export
cal_gamma <- function(Pa, lambda = 2.5) {
    Cp * Pa / (epsilon * lambda)
    # 0.000665 * pres
}

#' @rdname ET0_helper
#' @export
cal_bowen <- function(Tair, Pa = atm) {
    gamma = cal_gamma(Pa, cal_lambda(Tair))
    slope = cal_slope(Tair)
    gamma / slope
}

#' @param z Elevation above sea level `[m]`. Must provided if pres
#'          are not provided.
#' 
#' @rdname ET0_helper
#' @export
cal_Pa <- function(z = NULL) {
    101.3 * ((293.0 - (0.0065 * z)) / 293.0)**5.26 # Eq. 7
}

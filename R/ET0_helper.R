#' @name ET0_helper
#' @title helper functions for potential evapotranspiration (ET0)
#' 
#' @description
#' - `lambda`: latent heat of vaporization, about `[2.5 MJ kg-1]`.
#' - `slope`: The slope of the saturation vapour pressure curve at certain air
#'   temperature `Tair`, `[kPa degC-1]`.
#' - `gamma`: psychrometric constant (`[kPa degC-1]`), `Cp*Pa/(epsilon*lambda)`.
#' - `U2`: 10m wind speed (m/s). According to wind profile relationship, convert U_z to U_2.
#' - `es`: saturation vapor pressure (kPa)
#' - `ea`: actual vapor pressure (kPa)
#'      1. `RH_mean` | `Tmin`, `Tmax`: `(es(Tmax) + es(Tmin)) * RH_mean/200`
#'      2. `Tmin`                    : `es(Tmin)`
#' - `cal_TvK`: vitual temperature (K)
#'      1. Tair + q
#'      2. Tair + ea/Pa
#'      3. `1.01 * (Tair + 273)` , FAO56 Eq. 3-7
#' @references 
#' 1. Allen, R. G., & Luis S. Pereira. (1998). Crop
#'    evapotranspiration-Guidelines for computing crop water requirements-FAO
#'    Irrigation and drainage paper 56. European Journal of Agronomy, 34(3),
#'    144–152. doi:10.1016/j.eja.2010.12.001
#' @examples
#' cal_VPD(10, 5)
#' 
#' par(mfrow = c(2, 2), mar = c(3, 1.8, 2, 1), mgp = c(2, 0.6, 0))
#' Tair = -10:50
#' plot(cal_es(-10:50), main = "es (kPa)", xlab = "Tair")
#' plot(cal_bowen(-10:50), main = "Bowen ratio", xlab = "Tair")
#' plot(cal_slope(-10:50), main = "slope (kPa/degC)", xlab = "Tair")
#' plot(cal_gamma(-10:50), main = "gamma (kPa/degC)", xlab = "Tair")
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


#' @param Tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param Tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param RH_max,RH_mean,RH_min Daily max, mean and min relative humidity `[%]`.
#' 
#' @rdname ET0_helper
#' @export
cal_ea <- function(Tmin, Tmax = NULL, RH_mean = NULL) {
    # if(!is.null(RH_max) && !is.null(RH_min))
    #   return((cal_es(Tmax) * RH_min + cal_es(Tmin) * RH_max)/200)
    # if(!is.null(RH_max))
    #   return(cal_es(Tmin) * RH_max / 100)
    # if(is.null(RH_max) && is.null(RH_mean) && is.null(RH_min))
    #   return(cal_es(Tmin))
    if (!is.null(RH_mean)) {
        return((cal_es(Tmax) + cal_es(Tmin)) * RH_mean / 200)
    }
    return(cal_es(Tmin))
}
# ' 1. `RH_max`, `RH_min` | `Tmin`, `Tmax`: `(es(Tmax) * RH_min + es(Tmin) * RH_max)/200`
# ' 3. `RH_max`           | `Tmin`        : `es(Tmin) * RH_max / 100`

#' @param Tdew dew temperature, (`[deg Celsius]`)
#' @param Tair 2m air temperature, (`[deg Celsius]`)
#' 
#' @rdname ET0_helper
#' @export
cal_VPD <- function(Tair, Tdew = NULL) {
    ea <- cal_es(Tdew)
    es <- cal_es(Tair)
    es - ea
}

#' @rdname ET0_helper
#' @export
cal_lambda <- function(Tair) {
    # MJ kg-1
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
cal_gamma <- function(Tair, Pa = atm) {
    lambda = cal_lambda(Tair)
    Cp * Pa / (epsilon * lambda)
    # 0.000665 * pres
}

#' @rdname ET0_helper
#' @export
cal_bowen <- function(Tair, Pa = atm) {
    gamma = cal_gamma(Tair, Pa)
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


#' @rdname ET0_helper
#' @export
cal_rH <- function(U2, h = 0.12) {
    k <- 0.41 # Karman's constant
    d <- 2 / 3 * h
    z_om <- 0.123 * h
    z_oh <- 0.1 * z_om
    z_m <- z_h <- 2

    log((z_m - d) / z_om) * log((z_h - d) / z_oh) /
        (k^2 * U2)
}

#' @rdname ET0_helper
#' @export
cal_rH2 <- function(U2, Tair, Pa = atm) {
    # `f(U2) = 2.6 * (1 + 0.54U2)` is equivalent to Shuttleworth1993  
    # rou_a * Cp / rH = f(U2)
    f_U2 = 2.6 * (1 + 0.54*U2)
    rou_a = 3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5, kg m-3
    rH = rou_a * Cp / f_U2 * 86400
    rH
}

# vitual temperature
#' @rdname ET0_helper
#' @export
cal_TvK <- function(Tair, q = NULL, ea = NULL, Pa = atm) {
    if (is.null(q) && !is.null(ea)) {
        # ea
        (Tair + K0) * (1 + (1 - epsilon) * ea / Pa)
    } else if (!is.null(q) && is.null(ea)) {
        # q
        (Tair + K0) * (1 + (1 - epsilon) / epsilon * q) # to degK
    } else {
        1.01 * (Tair + 273) # Eq. 3-7
    }
}

#' @rdname ET0_helper
#' @export
cal_rou_a <- function(Tair, Pa = atm, q = NULL, ea = NULL) {
    # Tv = cal_TvK(Tair)
    # R = 0.287 # kJ kg-1 K-1
    3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5
}

#' Soil heat flux.
#'
#' Estimate monthly soil heat flux (G) from the mean air temperature of the 
#' previous, current or next month assuming as grass crop.
#'
#' @param t.p Mean air temperature of the previous month `[Celsius]`.
#' @param t.n Mean air temperature of the next month `[Celsius]`.
#' @param t.c Mean air temperature of the current month `[Celsius]`.
#'
#' @return Soil heat flux `[MJ m-2 day-1]`.
#' @export
soil_heat_flux <- function(t.p, t.n = NULL, t.c = NULL) {
  if (!is.null(t.n)) {
    return(0.07 * (t.n - t.p))
  }
  if (!is.null(t.c)) {
    return(0.14 * (t.c - t.p))
  }
  stop("Temperature of the next time step or the current time step must provide one.")
}

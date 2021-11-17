#' Estimate ET0 by FAO-56 Penman-Monteithe quation.
#' 
#' Estimate reference evapotranspiration (ET0) from a hypothetical short grass 
#' reference surface using the FAO-56 Penman-Monteithe quation 
#' (equation 6 in Allen et al. (1998))
#' 
#' @param Rs Incoming shortwave radiation at crop surface `[MJ m-2 day-1]`.
#'           Can be calculated by [cal_Rs].
#' @param tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param ws Wind speed `[m s-1]`.
#' @param G Soil heat flux `[MJ m-2 day-1]`. Normally set to 0.0 when the time
#'          steps are less than 10 days. Should be calculated in monthly or
#'          longer time step. Can be calculated using function [soil_heat_flux].
#'          Default 0.0.
#' @param h.ws Height where measure the wind speed `[m]`. Default 10m.
#' @param albedo Albedo of the crop as the proportion of gross incoming solar
#'               radiation that is reflected by the surface. Default 0.23
#'               (i. e. the value of the short grass reference crop defined
#'               by the FAO).
#' @param delta Slope of the saturation vapour pressure curve. Can be estimated
#'              by air temperature if is not provided.
#' @param gamma Psychrometric constant `[kPa Celsius -1]`. Can be estimated
#'              by the air pressure or elevation if is not provided.
#' @param z Elevation above sea level `[m]`. Must provided if gamma and pres
#'          are not provided.
#' @param ea Actual vapor pressure `[kPa]`. Can be estimated by maximum or minimum
#'           air temperature and mean relative humidity.
#' @param es Saturated vapor pressure `[kPa]`. Can be estimated by maximum or
#'           minimum air temperature if is not provided.
#' @param RH_mean Daily mean relative humidity `[%]`. If not provided then
#'               the actual vapor pressure would be estimated by minimum air
#'               temperature.
#' @param pres Air pressure `[kPa]`. If not provided, must provide z.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#'            calculate net outgoing longwave radiation than Rso.
#' @param Rso Clear sky incoming shortwave radiation, i. e. extraterrestrial
#'            radiation multiply by clear sky transmissivity (i. e. a + b + 2E-5*elev`[m]`,
#'            a + b are coefficients of Angstrom formula. Normally 0.75)
#'            `[MJ m-2 day-1]`. Ext. rad. can be calculated by [rad_ext].
#'            Should be provided if `cld` is not provided.
#'
#' @return Reference evapotranspiration ET0 from a hypothetical grass
#'         reference surface `[mm day-1]`.
#' 
#' @export
PET_pm <- function(Rs, tmax, tmin, ws, G = 0.0, h.ws = 10.0, albedo = 0.23,
                   delta = NULL, gamma = NULL, z = NULL, ea = NULL, es = NULL,
                   RH_mean = NULL, pres = NULL, Rso = NULL, cld = NULL,
                   tall_crop = FALSE){

  tas <- (tmax + tmin) / 2
  if(is.null(ea)) {
    if (!is.null(RH_mean)) {
      ea <- cal_ea(tmin, tmax, RH_mean = RH_mean)
    } else {
      # if no ea, let ea = es
      ea <- cal_es(tmin)
    }
  }

  if (is.null(es)) es <- (cal_es(tmax) + cal_es(tmin)) / 2
  if(is.null(delta)) delta <- delta_es(tas)

  Rln <- cal_Rln(tmax, tmin, ea, Rs, Rso, cld)
  Rn <- Rs * (1 - albedo) - Rln
  Rn[Rn < 0] <- 0
  
  if(is.null(gamma)) gamma <- cal_gamma(pres, z)

  if(tall_crop) {
    p1 <- 1600
    p2 <- 0.38
  } else {
    p1 <- 900
    p2 <- 0.34
  }

  u2 <- ws * (4.87 / log((67.8 * h.ws) - 5.42))

  (0.408 * delta * (Rn - G) +
      gamma * 900 / (tas + 273.15) * u2 * (es - ea)) /
    (delta + (gamma * (1 + 0.34 * u2)))
}

#' Estimate ET0 by Hargreaves equation.
#' 
#' @description Estimate reference evapotranspiration (ET0) from a hypothetical 
#' short grass reference surface using the the Hargreaves equation.
#' 
#' @param tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param tmean Daily mean air temperature at 2m height `[deg Celsius]`. If not
#' provided it would estimated by averaging the tmax and tmin.
#' 
#' @param Ra  Clear sky incoming shortwave radiation, i. e. extraterrestrial
#' radiation multiply by clear sky transmissivity (i. e. a + b, a and b are
#' coefficients of Angstrom formula. Normally 0.75) `[MJ m-2 day-1]`. If not
#' provided, must provide lat and dates. @param lat Latitude `[degree]`. @param
#' dates A R Date type of a vector of Date type. If not provided, it will Regard
#' the ssd series is begin on the first day of a year.
#' 
#' @return Reference evapotranspiration ET0 from a hypothetical grass reference
#' surface `[mm day-1]`.
#' 
#' @export
PET_hg <- function(tmax, tmin, tmean = NULL, Ra = NULL, lat = NULL, dates=NULL) {
  if(is.null(tmean)) tmean <- (tmax + tmin) / 2
  if(is.null(Ra)) Ra <- rad_ext(lat, dates)
  0.0023 * 0.408 * Ra * sqrt(tmax - tmin) * (tmean + 17.8)
}

#' Psychrometric constant
#'
#' psychrometric constant.
#'
#' @param pres Air pressure `[kPa]`.
#' @param z Elevation above sea level `[m]`. Must provided if pres
#'          are not provided.
#'
#' @return Psychrometric constant `[kPa degC-1]`.
#' @export
cal_gamma <- function(pres = NULL, z = NULL) {
  if (is.null(pres)) {
    if (is.null(z)) {
      stop("If pres is null, must provide z (elevation).")
    } else {
      pres <- 101.3 * ((293.0 - (0.0065 * z)) / 293.0)**5.26
    }
  }
  0.000665 * pres
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

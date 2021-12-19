#' Estimate ET0 by FAO-56 Penman-Monteithe quation.
#' 
#' Estimate reference evapotranspiration (ET0) from a hypothetical short grass 
#' reference surface using the FAO-56 Penman-Monteithe quation 
#' (equation 6 in Allen et al. (1998))
#' 
#' @param Rs Incoming shortwave radiation at crop surface `[MJ m-2 day-1]`.
#'           Can be calculated by [cal_Rs].
#' @param Tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param Tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param Uz Wind speed `[m s-1]`.
#' @param G Soil heat flux `[MJ m-2 day-1]`. Normally set to 0.0 when the time
#'          steps are less than 10 days. Should be calculated in monthly or
#'          longer time step. Can be calculated using function [soil_heat_flux].
#'          Default 0.0.
#' @param z.wind Height where measure the wind speed `[m]`. Default 10m.
#' @param albedo Albedo of the crop as the proportion of gross incoming solar
#'               radiation that is reflected by the surface. Default 0.23
#'               (i. e. the value of the short grass reference crop defined
#'               by the FAO).
#' @param delta Slope of the saturation vapour pressure curve. Can be estimated
#'              by air temperature if is not provided.
#' @param gamma Psychrometric constant `[kPa Celsius -1]`. Can be estimated
#'              by the air pressure or elevation if is not provided.
#' @param z Elevation above sea level `[m]`. Must provided if gamma and Pa
#'          are not provided.
#' @param ea Actual vapor pressure `[kPa]`. Can be estimated by maximum or minimum
#'           air temperature and mean relative humidity.
#' @param es Saturated vapor pressure `[kPa]`. Can be estimated by maximum or
#'           minimum air temperature if is not provided.
#' @param RH_mean Daily mean relative humidity `[%]`. If not provided then
#'               the actual vapor pressure would be estimated by minimum air
#'               temperature.
#' @param Pa Air pressure `[kPa]`. If not provided, must provide z.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#'            calculate net outgoing longwave radiation than Rso.
#' @param Rso Clear sky incoming shortwave radiation, i. e. extraterrestrial
#'            radiation multiply by clear sky transmissivity (i. e. a + b + 2E-5*elev`[m]`,
#'            a + b are coefficients of Angstrom formula. Normally 0.75)
#'            `[MJ m-2 day-1]`. Ext. rad. can be calculated by [cal_Ra].
#'            Should be provided if `cld` is not provided.
#'
#' @return Reference evapotranspiration ET0 from a hypothetical grass
#'         reference surface `[mm day-1]`.
#' 
#' @export
PET_pm <- function(date, Rs, Tmax, Tmin, Uz, z.wind = 10.0, 
  ea = NULL, es = NULL, RH_mean = NULL, Pa = NULL, Rso = NULL, cld = NULL,
  Z = 0, lat = NULL, 
  albedo = 0.23, tall_crop = FALSE)
{
  # 1. calculate VPD
  Tair <- (Tmax + Tmin) / 2
  if(is.null(ea)) ea <- cal_ea(Tmin, Tmax, RH_mean = RH_mean)
  if (is.null(es)) es <- (cal_es(Tmax) + cal_es(Tmin)) / 2
  D = es - ea
  
  # 2. calculate Rn
  ET0_FAO98(Rn, Tair, Pa, D, wind, z.wind = z.wind, tall_crop = tall_crop)
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

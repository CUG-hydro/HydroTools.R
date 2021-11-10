# https://github.com/Sibada/sibadaR/blob/master/R/FAO56.R

#' Actual vapour pressure (kPa)
#'
#' @param q specific humidity in kg/kg
#' @param p surface air pressure
#' 
#' @return e in the same unit as p
#' 
#' @references
#' https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html, Eq-17
#'
#' @export
vapour_press <- function(q, p) {
    q * p / (0.622 + 0.378 * q)
}

#' @export
#' @rdname vapour_press
q2ea <- vapour_press

#' @rdname vapour_press
#' @export
q2RH <- function(q, p, t) {
  ea <- vapour_press(q, p)
  es <- cal_es(t)
  ea / es * 100
}

#' @rdname vapour_press
#' @export
RH2q <- function(RH, t) {
  # ea <- vapour_press(q, p)
  es <- cal_es(t)
  ea <- es * RH/100
  # ea / es * 100
}

#' Estimate saturation vapor pressure.
#'
#' Estimate saturation vapor pressure from air temperature.
#' @param t Air temperature `[Celsius]`.
#' @return saturation vapor pressure `[kPa]`.
#' 
#' @export
#' @rdname vapour_press
cal_es <- function(t) {
    0.6108 * exp((17.27 * t) / (t + 237.3))
}

#' Estimate actual vapor pressure.
#'
#' @description Estimate actual vapor pressure by providing daily maximum and minimum
#'              air temperature, daily maximum, mean and minimum relative humidity.
#' 
#' @param tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param rhmax,rhmean,rhmin Daily max, mean and min relative humidity `[%]`.
#'
#' @details tmin must be provided. 
#' If tmax was not provided, avp will only estimated by tmin. If both rhmax and 
#' rhmin are provided, avp will be estimated by tmax, tmin, rhmax and rhmin.
#' 
#' If rhmin was not provided but provided rhmean, avp will estimated by rhmean, 
#' tmax and tmin.
#' 
#' If rhmin and rhmean are not provided, avp will estimated by rhmax, tmax and 
#' tmin.
#' 
#' If rhmax, rhmean and rhmin are not provided, avp will only estimated by tmin.
#'
#' @return Actual vapor pressure (i.e. avp or ea) `[kPa]`.
#' @export
cal_ea <- function(tmin, tmax = NULL, rhmax = NULL, rhmean = NULL, rhmin = NULL) {
  if(is.null(tmax))
    return(cal_es(tmin))
  
  if(!is.null(rhmax) && !is.null(rhmin))
    return((cal_es(tmax) * rhmin + cal_es(tmin) * rhmax)/200)

  if(!is.null(rhmean))
    return((cal_es(tmax) + cal_es(tmin)) * rhmean/200)

  if(!is.null(rhmax))
    return(cal_es(tmin) * rhmax / 100)

  if(is.null(rhmax) && is.null(rhmean) && is.null(rhmin))
    return(cal_es(tmin))
}

#' Estimate the slope of the saturation vapour pressure curve.
#'
#' @description Estimate the slope of the saturation vapour pressure curve
#'              by providing a air temperature.
#'
#' @param t Daily mean air temperature at 2m height `[deg Celsius]`.
#'
#' @return Slope of the saturation vapour pressure curve.
#' @export
delta_es <- function(t) {
  4098 * (0.6108 * exp((17.27 * t) / (t + 237.3))) /
    (t + 237.3)**2
}

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
#' @examples
#' RH = 90
#' p = atm
#' t = 30
#' q = RH2q(RH, p, t)
#' RH2 = q2RH(q, p, t)
#' e = q2ea(q, p)
#' es = cal_es(t)
#' @export
vapour_press <- function(q, p) {
    q * p / (epsilon + (1 - epsilon) * q)
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

#' @references
#' https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
#' 
#' @rdname vapour_press
#' @export
RH2q <- function(RH, p, t) {
  # ea <- vapour_press(q, p)
  es <- cal_es(t)
  ea <- es * RH/100

  # ws = epsilon * es / (p - es)
  w <- epsilon * ea / (p - ea)
  # w = RH/100 * ws 
  w2q(w) # q: g / g
}


#' `w`: mix ratio, m_v / m_d
#' @rdname vapour_press
#' @export
w2q <- function(w, p) w / (w + 1)

#' @rdname vapour_press
#' @export
q2w <- function(q, p) q / (1 - q)


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
#' @param RH_max,RH_mean,RH_min Daily max, mean and min relative humidity `[%]`.
#'
#' @details tmin must be provided. There are 4 options:
#' 1. `RH_max`, `RH_min` | `tmin`, `tmax`: `(es(tmax) * RH_min + es(tmin) * RH_max)/200`
#' 2. `RH_mean`          | `tmin`, `tmax`: `(es(tmax) + es(tmin)) * RH_mean/200`
#' 3. `RH_max`           | `tmin`        : `es(tmin) * RH_max / 100`
#' 4. `tmin`                             : `es(tmin)`
#' 
#' @return Actual vapor pressure (i.e. avp or ea) `[kPa]`.
#' @export
cal_ea <- function(tmin, tmax = NULL, RH_max = NULL, RH_min = NULL, RH_mean = NULL) {
  if(!is.null(RH_max) && !is.null(RH_min))
    return((cal_es(tmax) * RH_min + cal_es(tmin) * RH_max)/200)

  if(!is.null(RH_max))
    return(cal_es(tmin) * RH_max / 100)
  
  if(is.null(RH_max) && is.null(RH_mean) && is.null(RH_min))
    return(cal_es(tmin))

  if(!is.null(RH_mean))
    return((cal_es(tmax) + cal_es(tmin)) * RH_mean/200)
  # if(is.null(tmax))
  return(cal_es(tmin))  
}

#' The slope of the saturation vapour pressure curve.
#' 
#' the slope of the saturation vapour pressure curve by providing a air temperature.
#'
#' @param t Daily mean air temperature at 2m height `[deg Celsius]`.
#'
#' @return Slope of the saturation vapour pressure curve.
#' @export
delta_es <- function(t) {
  4098 * (0.6108 * exp((17.27 * t) / (t + 237.3))) /
    (t + 237.3)**2
}

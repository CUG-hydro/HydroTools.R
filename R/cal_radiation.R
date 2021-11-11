# #' Extraterrestrial radiation
# #' 
# #' @param lat latitude
# #' @param doy doy is the number of the day in the year between 1 (1 January) 
# #' and 365 or 366 (31 December)
# #' 
# #' @return Ra, MJ m-2 d-1
# #' @export 
# cal_Ra <- function(lat, doy) {
#     # solar declination, rad (1 rad = 57.2957795 deg)
#     delta <- 0.409*sin(0.0172*doy-1.39)
#     # relative distance Earth-Sun, []
#     dr <- 1 + 0.033*cos(0.0172*doy)
#     # sunset hour angle, rad
#     latr <- lat/57.2957795 # (180/pi)
#     sset <- -tan(latr)*tan(delta)
    
#     omegas <- sset*0
#     omegas[abs(sset)<=1] <- acos(sset[abs(sset)<=1])
#     # correction for high latitudes
#     omegas[sset<(-1)] <- max(omegas)
#     # Ra, MJ m-2 d-1
#     Ra <- 37.6*dr*(omegas*sin(latr)*sin(delta)+cos(latr)*cos(delta)*sin(omegas))
#     Ra <- ifelse(Ra<0,0,Ra)
#     Ra
# }

#' Estimate daily daily extraterrestrial radiation.
#' 
#' @description Estimate daily daily extraterrestrial radiation `[MJ m-2 day-1]`.
#' 
#' @param lat Latitude `[degree]`.
#' @param dates A R Date type of a vector of Date type. If not provided, it will
#'              Regard the ssd series is begin on the first day of a year.
#' 
#' @return extraterrestrial radiation `[MJ m-2 day-1]`.
#' @export
cal_Ra <- function(lat, dates) {
  lat <- lat * pi/180
  if (is.numeric(dates)) {
    J <- dates
  } else {
    J <- as.double(format(dates, '%j'))
  }

  dr <- 1 + 0.033 * cos(pi * J/182.5)
  delta <- 0.408 * sin(pi * J/182.5 - 1.39)

  tmpv <- tan(lat) * tan(delta)
  tmpv[tmpv > 1.] <- 1.
  tmpv[tmpv < -1.] <- -1.

  ws <- acos(-tmpv)
  Ra <- 118.08 * dr/pi * (ws * sin(lat) * sin(delta) + cos(lat) *
                            cos(delta) * sin(ws))
  Ra
}

#' @export
#' @rdname cal_Ra
ext_rad <- cal_Ra

#' Daily inward shortwave solar radiation.
#'
#' Estimate daily solar radiation at crop surface `[MJ m-2 day-1]` by providing 
#' sunshine duration (SSD) in hours or cloud cover in fraction.
#' 
#' @param ssd sunshine duration [hours]. 
#' @param lat Latitude `[degree]`. 
#' @param dates A R Date type of a vector of Date type. If not provided, it will 
#' regard the ssd series is begin on the first day of a year. 
#' @param a Coefficient of the Angstrom formula. Determine the relationship 
#' between ssd and radiation. Default 0.25. 
#' @param b Coefficient of the Angstrom formula. Default 0.50.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#' calculate solar radiation rather than SSD and parameter a and b.
#'
#' @return Solar radiation at crop surface `[MJ m-2 day-1]`.
#'
#' @references 
#' Martinezlozano J A, Tena F, Onrubia J E, et al. The historical
#' evolution of the Angstrom formula and its modifications: review and
#' bibliography.[J]. Agricultural & Forest Meteorology, 1985, 33(2):109-128.
#' @export
cal_Rs <- function(ssd, lat, dates = NULL, a = 0.25, b = 0.5, cld = NULL) {
  lat <- lat*pi/180
  if (is.null(dates)) {
    J <- rep_len(c(1:365, 1:365, 1:365, 1:366), length(ssd))
  } else {
    J <- as.double(format(dates, '%j'))
  }

  delta <- 0.408 * sin(pi * J/182.5 - 1.39)
  tmpv <- tan(lat*pi/180) * tan(delta)
  tmpv[tmpv > 1.] <- 1.; tmpv[tmpv < -1.] <- -1.
  ws <- acos(-tmpv)

  Ra <- 24*60*0.082 * dr/pi * (ws * sin(lat) * sin(delta) + cos(lat) *
                          cos(delta) * sin(ws))
  if(!is.null(cld)) {
    Rs <- (1. - cld) * Ra
  } else {
    N <- ws * 24/pi
    Rs <- (a + b * ssd/N) * Ra
  }
  Rs
}

#' Net outgoing longwave radiation.
#' 
#' Net outgoing longwave radiation.
#'
#' @param tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param Rs Incoming shortwave radiation at crop surface `[MJ m-2 day-1]`.
#' @param ea Actual vapor pressure `[kPa]`. Can be estimated by maximum or minimum
#'           air temperature and mean relative humidity.
#' @param Rso Clear sky incoming shortwave radiation, i. e. extraterrestrial
#'            radiation multiply by clear sky transmissivity (i. e. a + b,
#'            a and b are coefficients of Angstrom formula. Normally 0.75)
#'            `[MJ m-2 day-1]`.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#'            calculate net outgoing longwave radiation than Rso.
#'
#' @return Net outgoing longwave radiation `[MJ m-2 day-1]`
#'
#' @export
cal_Rln <- function(tmax, tmin, ea, Rs = NULL, Rso = NULL, cld = NULL) {
  if((is.null(Rs) || is.null(Rso)) && is.null(cld))
    stop("Must provide `Rs` and `Rso`, or provide `cld`.")
  if(is.null(cld)){
    cld <- 1. - Rs / Rso
    cld[is.na(cld)] <- 1.
  }
  (4.093e-9 * (((tmax+273.15)**4 + (tmin+273.15)**4) / 2)) *
    (0.34 - (0.14 * sqrt(ea))) *
    (1.35 * (1. - cld) - 0.35)
}

#' Estimate Incoming longwave radiation.
#'
#' @description Estimate Incoming longwave radiation. Not included in FAO56 
#' paper but added for convinience.
#' 
#' @param temp Near surface air temperature `[deg Celsius]`.
#' @param ea Near surface actual vapour pressure `[kPa]`.
#' @param s Cloud air emissivity, the ratio between acutal incomming shortwave
#' radiation and clear sky incoming shortwave radiation.
#' @param method The method to esteimate the air emissivity. Must be one of 
#' 'MAR', 'SWI', 'IJ', 'BRU', 'T', and 'KON'.
#' 
#' @return Incomming longwave radiation `[W /m2]`.
#' @export
cal_Rlin <- function(temp, ea = NULL, s = 1, method = 'KON') {
  if(!(method %in% c('MAR', 'SWI', 'IJ', 'BRU', 'SAT', 'KON')))
    stop("method must be one of 'MAR', 'SWI', 'IJ', 'BRU', 'SAT', and 'KON'.")
  ea <- ea * 10
  temp <- temp + 273.15
  if(method == 'MAR'){
    ep_ac <- 0.5893 + 5.351e-2 * sqrt(ea*10)
  }
  if(method == 'SWI'){
    ep_ac <- 9.294e-6 * temp*temp
  }
  if(method == 'IJ'){
    ep_ac <- 1 - 0.26 * exp(-7.77e-4 * (273 - temp)**2)
  }
  if(method == 'BRU'){
    ep_ac <- 1.24 * (ea/temp)**(1/7)
  }
  if(method == 'SAT'){
    ep_ac <- 1.08 * exp(-ea ** (temp/2016))
  }
  if(method == 'KON'){
    ep_ac <- 0.23 + 0.848 * (ea/temp)**(1/7)
  }
  ep_a <- 1 - s + s*ep_ac
  ep_a * 5.67e-8 * temp**4
}

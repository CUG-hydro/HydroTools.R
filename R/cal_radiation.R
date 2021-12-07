#' @name cal_ssd
#' @title sunset hour angle
#' @description Calculating sunset hour angle (`ws`) according to Allen Eq. 25.
#' @inheritParams cal_Ra
#' @export
NULL

#' @rdname cal_ssd
#' @export
cal_ws <- function(lat, J) {
  lat %<>% deg2rad()
  sigma <- 0.409 * sin(pi * J / 182.5 - 1.39) # Allen, Eq. 24

  tmp <- clamp(-tan(lat) * tan(sigma), c(-1, 1))
  ws <- acos(tmp) # Eq. 25
  ws
}

#' @rdname cal_ssd
#' @export
ws2ssd <- function(ws) ws / pi * 24 # Ge ChaoXiao, Eq. 2-18

#' @rdname cal_ssd
#' @export
cal_ssd <- function(lat, J) cal_ws(lat, J) %>% ws2ssd()

#' Daily extraterrestrial radiation.
#'
#' @description Estimate daily extraterrestrial radiation `[MJ m-2 day-1]`.
#'
#' @param lat Latitude `[degree]`.
#' @param J `doy` vector.
#'
#' @return extraterrestrial radiation `[MJ m-2 day-1]`.
#' @seealso [cal_Rs()]
#' 
#' @examples
#' cal_Ra(20, 1)
#' # [1] 25.84874
#' @importFrom lubridate yday is.Date
#' @export
cal_Ra <- function(lat, J) {
  J %<>% check_doy()
  dr <- 1 + 0.033 * cos(pi * J / 182.5) # Allen, Eq. 23
  sigma <- 0.409 * sin(pi * J / 182.5 - 1.39) # Allen, Eq. 24

  ws <- cal_ws(lat, J)
  # 24 * 60 * 0.082 = 118.08
  lat %<>% deg2rad()
  Ra <- 118.08 * dr / pi * (ws * sin(lat) * sin(sigma) + cos(lat) * cos(sigma) * sin(ws)) # Allen, Eq. 21
  Ra <- ifelse(Ra < 0, 0, Ra)
  Ra
}

#' Daily inward shortwave solar radiation.
#'
#' Daily inward shortwave solar radiation at crop surface `[MJ m-2 day-1]` by 
#' providing sunshine duration (SSD) in hours or cloud cover in fraction.
#' 
#' @param lat Latitude `[degree]`.
#' @param dates A R Date type of a vector of Date type. If not provided, it will 
#' regard the ssd series is begin on the first day of a year. 
#' @param ssd sunshine duration [hours]. If `ssd = NULL`, `Rs` is the clear-sky 
#' solar radiation.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#' calculate solar radiation rather than SSD and parameter a and b.
#' @param a Coefficient of the Angstrom formula. Determine the relationship 
#' between ssd and radiation. Default 0.25. 
#' @param b Coefficient of the Angstrom formula. Default 0.50.
#' 
#' @return A data.table, solar radiation at crop surface `[MJ m-2 day-1]`.
#' - `Ra`: Clear-sky surface downward shortwave radiation
#' - `Rs`: Surface downward shortwave radiation
#' 
#' @references 
#' 1. Martinezlozano J A, Tena F, Onrubia J E, et al. The historical
#'    evolution of the Angstrom formula and its modifications: review and
#'    bibliography.[J]. Agricultural & Forest Meteorology, 1985, 33(2):109-128.
#' 2. https://github.com/sbegueria/SPEI/blob/master/R/penman.R
#' 
#' @seealso [cal_Ra()]
#' @export
cal_Rs <- function(lat, J, ssd = NULL, cld = NULL, Z = 0, 
  albedo = 0.23, 
  a = 0.25, b = 0.5) 
{
  J %>% check_doy()
  Ra = cal_Ra(lat, J)

  if(!is.null(cld)) {
    nN = (1 - cld) 
  } else {
    # N <- ws * 24/pi # Ge ChaoXiao, Eq. 2-18
    N <- cal_ssd(lat, J)
    if (is.null(ssd)) ssd = N
    nN = ssd / N
  }
  # TODO: check the equation of `cld` influence
  # Rs <- (1. - cld) * Ra
  # Rs <- (1. - cld) * Ra * (a + b) # also named asR_so
  Rs <- (a + b * nN) * Ra
  # Rns <- (1 - albedo) * Rs
  Rso <- (0.75 + 2 * Z / (10^5)) * Ra # 计算晴空辐射, FAO56, Eq. 37
  data.table(Rs, Rso, Ra)
}


#' Net outgoing longwave radiation.
#' 
#' Net outgoing longwave radiation.
#' 
#' @param tmax Daily maximum air temperature at 2m height `[deg Celsius]`.
#' @param tmin Daily minimum air temperature at 2m height `[deg Celsius]`.
#' @param ea Actual vapor pressure `[kPa]`. Can be estimated by maximum or minimum
#'           air temperature and mean relative humidity.
#' @param Rs Incoming shortwave radiation at crop surface `[MJ m-2 day-1]`.
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


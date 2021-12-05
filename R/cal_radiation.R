#' @name cal_ssd
#' @title sunset hour angle
#' @description Calculating sunset hour angle (`ws`) according to Allen Eq. 25.
#' @inheritParams rad_ext
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
rad_ext <- function(lat, J) {
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

#' @export
#' @rdname rad_ext
cal_Ra <- rad_ext

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
#' @param a Coefficient of the Angstrom formula. Determine the relationship 
#' between ssd and radiation. Default 0.25. 
#' @param b Coefficient of the Angstrom formula. Default 0.50.
#' @param cld Cloud cover `[fraction]`. If provided it would be directly used to
#' calculate solar radiation rather than SSD and parameter a and b.
#' 
#' @return Solar radiation at crop surface `[MJ m-2 day-1]`.
#' 
#' @references 
#' 1. Martinezlozano J A, Tena F, Onrubia J E, et al. The historical
#'    evolution of the Angstrom formula and its modifications: review and
#'    bibliography.[J]. Agricultural & Forest Meteorology, 1985, 33(2):109-128.
#' 2. https://github.com/sbegueria/SPEI/blob/master/R/penman.R
#' 
#' @seealso [cal_Ra()]
#' @export
rad_surf <- function(lat, J, ssd = NULL, a = 0.25, b = 0.5, cld = NULL) {
  J %>% check_doy()
  Ra = cal_Ra(lat, J)

  if(!is.null(cld)) {
    # Rs <- (1. - cld) * Ra
    # Rs <- (1. - cld) * Ra * (a + b) # also named asR_so
    nN = (1 - cld) 
  } else {
    # N <- ws * 24/pi # Ge ChaoXiao, Eq. 2-18
    N <- cal_ssd(lat, J)
    if (is.null(ssd)) ssd = N
    nN = ssd / N
  }
  Rs <- (a + b * nN) * Ra
  Rs
}

#' @export
#' @rdname rad_surf
cal_Rs <- rad_surf


#' Surface net shortwave radiation
#' 
#' @param J integer, day of year
#' @param lat float, latitude 
#' @param ssd numeric vector, sun shine duration (hour)
#' @param RH numeric vector, relative humidity (%)
#' @param Tmin,Tmax numeric vector, min and max 2m-air temperature (degC)
#' @param cld (optional) cloud coverage (0-1). At least one of `cld` and `ssd` 
#' should be provided. If `cld` is not null, `ssd` will be ignored.
#' 
#' @param Z (optional) elevation (m), for the calculation of `Rso`
#' @param albedo (optional), `Rsn = (1 - albedo) Rs`
#' 
#' @return radiation in `[MJ d-1]`
#' - `Rn`  : Surface net radiation
#' - `Rln` : Surface outward net longwave radiation (negative means outgoing)
#' - `Rsn` : Surface downward net shortwave radiation
#' - `Rs`  : Surface downward shortwave radiation
#' - `Rso` : Clear-sky surface downward shortwave radiation
#' - `Ra`  : Extraterrestrial radiation
#' 
#' @author
#' Xie YuXuan and Kong Dongdong
#' 
#' @examples
#' cal_Rn(lat = 30, J = 1, ssd = 10, RH = 70, Tmin = 20, Tmax = 30)
#' 
#' @export
cal_Rn <- function(lat, J, ssd = 10, RH, Tmin, Tmax, 
  cld = NULL, 
  albedo = 0.23, Z = 0, ...) 
{ 
  dlt <- 0.409 * sin(J * 2 * pi / 365 - 1.39) # 太阳磁偏角delta,J为日序(1-365或366)
  ws <- acos(-tan(lat*pi/180)*tan(dlt))       # 时角; 日出时角/日落时角
  
  N <- 24*ws/pi # 最大可能日照时数(h)  
  dr <- 1+0.033*cos(J*2*pi/365) # 日地平均距离
  
  # 星际辐射/日地球外辐射:大气圈外接受的阳光辐射能量（可查表）
  as    <- 0.25  # 星际辐射回归系数
  bs    <- 0.5   # 星际辐射到地球的衰减系数
  Gsc   <- 0.082 # 太阳常数(MJ/m^2*min)

  Ra  <- 24*60/pi*Gsc*dr*(ws*sin(lat*pi/180)*sin(dlt)+cos(lat*pi/180)*cos(dlt)*sin(ws)) 

  if (is.null(cld)) {
    nN = ssd / N
  } else {
    nN = 1 - cld
  }

  Rs  <- (as + bs * nN) * Ra # Rs为太阳辐射(n/N日照系数)
  Rsn <- (1 - albedo) * Rs # 太阳净辐射
  Rso <- (0.75 + 2 * Z / (10^5)) * Ra # 计算晴空辐射, FAO56, Eq. 37

  ## longwave radiation 
  es = (cal_es(Tmax) + cal_es(Tmin))/2
  ea = es * RH/100
  
  T0 = 273.15
  # net outgoing longwave radiation [MJ m-2 day-1]
  sigma <- 4.903 / (10^9) #  斯蒂芬-玻尔兹曼常数, [MJ K-4 m-2]
  Rln = -sigma * ( ((Tmax + T0)^4 + (Tmin + T0)^4) / 2) * 
    (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rs / Rso - 0.35) 
  Rn = Rsn + Rln 
  data.table(Rn, Rsn, Rln, Rs, Rso, Ra)
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


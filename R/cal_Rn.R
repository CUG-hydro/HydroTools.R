#' Surface net radiation
#' 
#' @param J integer, day of year
#' @param lat float, latitude
#' @param Tmin,Tmax numeric vector, min and max 2m-air temperature (degC)
#'
#' @param ea Actual vapor pressure (kPa)
#' @param RH Relative humidity (%). If `ea` provided, `RH` will be ignored.
#' 
#' @param ssd numeric vector, sun shine duration (hour)
#' @param cld (optional) cloud coverage (0-1). At least one of `cld` and `ssd`
#' should be provided. If `cld` is not null, `ssd` will be ignored.
#' @param Rs (optional) Surface downward shortwave radiation (MJ m-2 d-1). 
#' - If not provided, `Rs` will be calculated by `(as + bs * nN) * Ra`.
#' - If provided, `ssd` and `cld` will be ignored.
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
#' @note
#' `Rn` might <= 0. Users need to constrain the min value by `pmax(Rn, 0)`.
#' 
#' @author
#' Xie YuXuan and Kong Dongdong
#' @seealso [cal_Ra()], [cal_Rs()]
#' 
#' @examples
#' cal_Rn(lat = 30, J = 1, RH = 70, Tmin = 20, Tmax = 30, ssd = 10)
#' @export
cal_Rn <- function(lat, J, Tmin, Tmax,
  ea = NULL, RH,
  ssd = NULL, cld = NULL,
  Rs = NULL, 
  albedo = 0.23, Z = 0, ...)
{
  J %<>% check_doy()
  # J = yday(date)
  dlt <- 0.409 * sin(J * 2 * pi / 365 - 1.39) # 太阳磁偏角delta,J为日序(1-365或366)
  ws <- acos(-tan(lat*pi/180)*tan(dlt))       # 时角; 日出时角/日落时角

  N <- 24*ws/pi # 最大可能日照时数(h)
  dr <- 1+0.033*cos(J*2*pi/365) # 日地平均距离

  # 星际辐射/日地球外辐射:大气圈外接受的阳光辐射能量（可查表）
  as  <- 0.25  # 星际辐射回归系数
  bs  <- 0.5   # 星际辐射到地球的衰减系数
  Gsc <- 0.082 # 太阳常数(MJ/m^2*min)

  Ra  <- 24*60/pi*Gsc*dr*(ws*sin(lat*pi/180)*sin(dlt)+cos(lat*pi/180)*cos(dlt)*sin(ws))
  
  if (is.null(Rs)) {
    if (!is.null(cld)) {
      nN = (1 - cld)
    } else {
      N <- ws/pi * 24 # Ge ChaoXiao, Eq. 2-18
      # N <- cal_ssd(lat, J)
      if (is.null(ssd)) ssd = N
      nN = ssd / N
    }
    
    # Rs <- (1. - cld) * Ra
    # Rs <- (1. - cld) * Ra * (a + b) # also named as R_so
    ## https://github.com/sbegueria/SPEI/blob/master/R/penman.R
    Rs  <- (as + bs * nN) * Ra # Rs为太阳辐射(n/N日照系数)
  }

  Rsn <- (1 - albedo) * Rs # 太阳净辐射
  Rso <- (0.75 + 2 * Z / (10^5)) * Ra # 计算晴空辐射, FAO56, Eq. 37

  ## longwave radiation
  if (is.null(ea)) {
    es = (cal_es(Tmax) + cal_es(Tmin)) / 2
    ea = es * RH / 100
  }

  T0 = 273.15
  # net outgoing longwave radiation [MJ m-2 day-1]
  sigma <- 4.903 / (10^9) #  斯蒂芬-玻尔兹曼常数, [MJ K-4 m-2]
  Rln = -sigma * ( ((Tmax + T0)^4 + (Tmin + T0)^4) / 2) *
    (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rs / Rso - 0.35)
  Rn = Rsn + Rln 
  # Rn = pmax(Rsn + Rln, 0) # make sure Rn not less zero
  # whether to use this constraint, it is determined by the user.
  data.table(Rn, Rsn, Rln, Rs, Rso, Ra)
}

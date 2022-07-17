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
#' @param Rsi (optional) Surface downward shortwave radiation (MJ m-2 d-1). 
#' - If not provided, `Rsi` will be calculated by `(as + bs * nN) * Rsi_toa`.
#' - If provided, `ssd` and `cld` will be ignored.
#' 
#' @param Z (optional) elevation (m), for the calculation of `Rsi_o`
#' @param albedo (optional), `Rsn = (1 - albedo) Rsi`
#' @param ... ignored
#' 
#' @return radiation in `[MJ d-1]`
#' - `Rn`  : Surface net radiation
#' - `Rln` : Surface outward net longwave radiation (negative means outgoing)
#' - `Rsn` : Surface downward net shortwave radiation
#' - `Rsi`  : Surface downward shortwave radiation (Rsi)
#' - `Rsi_o` : Clear-sky surface downward shortwave radiation
#' - `Rsi_toa`  : Extraterrestrial radiation (top of atmosphere, Rsi_toa)
#' 
#' @note
#' `Rn` might <= 0. Users need to constrain the min value by `pmax(Rn, 0)`.
#' 
#' @author
#' Xie YuXuan and Kong Dongdong
#' @seealso [cal_Rsi_toa()], [cal_Rsi()]
#' 
#' @examples
#' cal_Rn(lat = 30, J = 1, RH = 70, Tmin = 20, Tmax = 30, ssd = 10)
#' @export
cal_Rn <- function(lat, J, Tmin, Tmax,
  ea = NULL, RH,
  ssd = NULL, cld = NULL,
  Rsi = NULL, 
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

  Rsi_toa  <- 24*60/pi*Gsc*dr*(ws*sin(lat*pi/180)*sin(dlt)+cos(lat*pi/180)*cos(dlt)*sin(ws))
  
  if (is.null(Rsi)) {
    if (!is.null(cld)) {
      nN = (1 - cld)
    } else {
      N <- ws/pi * 24 # Ge ChaoXiao, Eq. 2-18
      # N <- cal_ssd(lat, J)
      if (is.null(ssd)) ssd = N
      nN = ssd / N
    }
    
    # Rsi <- (1. - cld) * Rsi_toa
    # Rsi <- (1. - cld) * Rsi_toa * (a + b) # also named as R_so
    ## https://github.com/sbegueria/SPEI/blob/master/R/penman.R
    Rsi  <- (as + bs * nN) * Rsi_toa # Rs为太阳辐射(n/N日照系数), Rsi
  }

  Rsn <- (1 - albedo) * Rsi # 太阳净辐射
  Rsi_o <- (0.75 + 2 * Z / (10^5)) * Rsi_toa # 计算晴空辐射, FAO56, Eq. 37

  ## longwave radiation
  if (is.null(ea)) {
    es = (cal_es(Tmax) + cal_es(Tmin)) / 2
    ea = es * RH / 100
  }

  T0 = 273.15
  # net outgoing longwave radiation [MJ m-2 day-1]
  sigma <- 4.903 / (10^9) #  斯蒂芬-玻尔兹曼常数, [MJ K-4 m-2]
  Rln = -sigma * ( ((Tmax + T0)^4 + (Tmin + T0)^4) / 2) *
    (0.34 - 0.14 * sqrt(ea)) * (1.35 * Rsi / Rsi_o - 0.35)
  Rn = Rsn + Rln 
  # Rn = pmax(Rsn + Rln, 0) # make sure Rn not less zero
  # whether to use this constraint, it is determined by the user.
  data.table(Rn, Rsn, Rln, Rsi, Rsi_o, Rsi_toa) # Rsi, 
}

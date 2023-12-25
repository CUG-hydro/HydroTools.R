
#' @name cal_ssd
#' @title sunset hour angle
#' @description Calculating sunset hour angle (`ws`) according to Allen Eq. 25.
#' @inheritParams cal_Rsi_toa
#' @export
NULL

#' @rdname cal_ssd
#' @export
cal_sunset_angle <- function(lat, J) {
  lat %<>% deg2rad()
  sigma <- 0.409 * sin(pi * J / 182.5 - 1.39) # Allen, Eq. 24

  tmp <- clamp(-tan(lat) * tan(sigma), c(-1, 1))
  ws <- acos(tmp) # Eq. 25
  ws
}

#' @param ws sunset hour angle
#'
#' @rdname cal_ssd
#' @export
ws2ssd <- function(ws) ws / pi * 24 # Ge ChaoXiao, Eq. 2-18

#' @rdname cal_ssd
#' @export
cal_ssd <- function(lat, J) cal_sunset_angle(lat, J) %>% ws2ssd()

# ------------------------------------------------------------------------------

# 赤纬角
#' Sun Declination angle
#' 
#' @importFrom lubridate yday dhours dminutes
#' @export 
get_sigma <- function(J) {
  if (inherits(J, "Date")) J %<>% yday()
  sigma <- 0.409 * sin(2 * pi / 365 * J - 1.39)
  rad2deg(sigma)
}

local2UTC <- function() {
}

deltaT_UTC <- function(lon, timeZone = 8) {
  lon_center <- timeZone * 15
  delta_lon <- (lon - 120)
  delta_minute <- round(delta_lon / 15 * 60, 2)
  dminutes(delta_minute)
}

#' Get local time
#'
#' @param time POSIXct object
#' @param lon longitude in deg
#'
#' @export
get_localtime <- function(time, lon, timeZone = 8) {
  time + deltaT_UTC(lon, timeZone)
}

# 太阳高度角
#' @export 
SunAngle_Elevation <- function(lat, sigma, omega) {
  lat %<>% deg2rad()
  sigma %<>% deg2rad()
  omega %<>% deg2rad()

  ans <- cos(lat) * cos(sigma) * cos(omega) + sin(lat) * sin(sigma)
  asin(ans) %>%
    rad2deg() %>%
    round(2)
}

#' @export 
SunAngle_Azimuth <- function(lat, sigma, omega) {
  lat %<>% deg2rad()
  sigma %<>% deg2rad()
  omega %<>% deg2rad()

  dN <- sin(sigma) * cos(lat) - sin(lat) * cos(sigma) * cos(omega)
  dE <- -cos(sigma) * sin(omega)
  ans <- dE / dN
  atan(dE / dN) %>% rad2deg()
  angle <- acos(dN / sqrt(dN^2 + dE^2)) * sign(dN)
  if (angle < 0) {
    angle <- angle + 2 * pi
  }
  angle <- round(rad2deg(angle), 2) # radian to deg
  angle
}


fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

#' suncalc
#' 
#' The default location is ZuoLing.
#' 
#' @export
suncalc <- function(time, lon = 114.6053, lat = 30.49694, ..., 
  year = year(Sys.time()), 
  verbose = TRUE) 
{
  if (is.integer(time)) {
    J = time
    time_str = sprintf("%d%03d", year, J)
    time <- as.Date(time_str, "%Y%j")
  }
  J <- yday(time)
  angle_sigma <- get_sigma(J)

  delta_minute = deltaT_UTC(lon, timeZone = 8)

  # 我们这里比北京时间晚了21.6min
  time_local <- get_localtime(time, lon)
  time_noon <- format(time_local, "%Y-%m-%d 12:00:00") %>% as.POSIXct()

  omega_mins <- as.numeric(difftime(time_local, time_noon, units = "mins"))
  omega <- omega_mins / 60 * 15

  angle_elev <- SunAngle_Elevation(lat, angle_sigma, omega)
  angle_azimuth <- SunAngle_Azimuth(lat, angle_sigma, omega)

  ws <- cal_sunset_angle(lat, yday(time))
  ssd <- ws2ssd(ws)
  time_begin <- time_noon - dhours(ssd / 2)
  time_end <- time_noon + dhours(ssd / 2)

  if (verbose) {
    fprintf("%-13s: %-10.2f\n", "赤纬角", angle_sigma)
    fprintf("%-13s: %-10.2f\n", "太阳高度角", angle_elev)
    fprintf("%-13s: %-10.2f\n", "太阳方位角", angle_azimuth)

    fprintf("---\n")
    fprintf("%-19s: %-10s\n", "北京时间", time)
    fprintf("%-19s: %-10s\n", "当地时间", time_local)

    fprintf("---\n")
    fprintf("%-14s: %-10s\n", "日出时间(local)", time_begin)
    fprintf("%-14s: %-10s\n", "日落时间(local)", time_end)

    fprintf("---\n")
    fprintf("%-14s: %-10s\n", "日出时间(UTC8) ", time_begin - delta_minute)
    fprintf("%-14s: %-10s\n", "日落时间(UTC8) ", time_end - delta_minute)
  }
  
  data.table(
    angle_elev, angle_azimuth, angle_sigma, delta_minute,
    time_local, time_begin, time_end)
}

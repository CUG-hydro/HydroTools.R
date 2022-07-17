#' @name unit_convert
#' @title unit convertion
#' 
#' @description 
#' - `R2Q`: Q(m^3/s) = R(mm) * area (km^2) * 1000/(3600*24)
#' - `W2mm`: 1 MJ /m2/day  =0.408 mm /day, 1 Watt /m2 = 0.0864 MJ /m2/day
#' 
#' @param Q runoff flow (m^3/s)
#' @param R runoff (mm per dt)
#' @param area basin area (km^2), if area = ,
#' @param dt time duration (hour)
NULL

#' @rdname unit_convert
#' @export 
R2Q <- function(R, area = dt*3.6, dt = 24) {
    # area = dt * 3.6
    R*area/(dt*3.6)
}

#' @rdname unit_convert
#' @export 
Q2R <- function(Q, area = dt*3.6, dt = 24) {
    # area = dt * 3.6
    Q * dt * 3.6 / area
}

# W = J/s
#' @rdname unit_convert
#' @export
MJ_2W <- function(x) {
    x / 86400 * 1e6
}

# ' MJ m-2 day-1 equivalent evaporation (mm/day)
# ' 
# ' 1 mm day-1 = 2.45 MJ m-2 day-1
#' @references http://www.fao.org/3/X0490E/x0490e0i.htm
#' @rdname unit_convert
#' @export
MJ_2mm <- function(x) {
    x/2.45
}

#' @rdname unit_convert
#' @export
W2_MJ <- function(x) {
    x / 1e6 * 86400
}

#' @param tmean daily mean temperature
#' 
#' @rdname unit_convert
#' @export 
W2mm <- function(x, tmean = 0) {
    Cp <- 4.2 * 0.242 # specific heat at constant pressure, 1.013 [kJ kg-1 0C-1]
    lamada <- 2500 - 2.2 * tmean
    x / lamada * 86400 * 10^-3 # W M-2 to mm
}

#' @param x scalar or numeric vector
#' 
#' @rdname unit_convert
#' @export
K2T <- function(x) x - 273.15

#' @rdname unit_convert
#' @export
T2K <- function(x) x + 273.15

#' @rdname unit_convert
#' @export
deg2rad <- function(deg) deg / 180 * pi

#' @param rad angle in radian
#' @param deg angle in degree
#' 
#' @rdname unit_convert
#' @export
rad2deg <- function(rad) rad / pi * 180

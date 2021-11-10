#' @name unit_convert
#' @title unit convertion
#' 
#' @description 
#' - `R2Q`: Q(m^3/s) = R(mm) * area (km^2) * 1000/(3600*24)
#' - `W2mm`: 1 MJ /m2/day  =0.408 mm /day, 1 Watt /m2 = 0.0864 MJ /m2/day
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

#' @rdname unit_convert
#' @export
W2_MJ <- function(x) {
    x / 1e6 * 86400
}

#' @rdname unit_convert
#' 
#' @export 
W2mm <- function(x, tmean = 0) {
    Cp <- 4.2 * 0.242 # specific heat at constant pressure, 1.013 [kJ kg-1 0C-1]
    lamada <- 2500 - 2.2 * tmean
    x / lamada * 86400 * 10^-3 # W M-2 to mm
}

#' @rdname unit_convert
#' @export
K2T <- function(x) x - 273.15

#' @rdname unit_convert
#' @export
T2K <- function(x) x + 273.15

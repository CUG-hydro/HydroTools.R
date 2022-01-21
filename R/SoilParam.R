#' Soil properties in SPAW
#' @name SPAW
#' 
#' @description
#' Wilting point, field capacity and saturated moisture
#' 
#' * `wilting point`: 1500 kPa moisture, %v, (Eq. 1)
#' * `field capacity`: 33 kPa moisture, %v, (Eq. 2)
#' * `saturated moisture`: Saturation (0 kPa moisture), %v, (Eq. 5)
#' 
#' @param S the weight ratio of Sand, (weight ratio, 0-1)
#' @param C the weight ratio of Clay, (weight ratio, 0-1)
#' @param OM the percentage of Organic Matter, (%w, 0-100)
#'
#' @note
#' The unit of `OM` is `%w`, which is different from `S` and `C`.
#' @return %v
#' 
#' @references 
#' 1. Saxton, K. E., & Rawls, W. J. (2006). Soil water characteristic estimates
#'    by texture and organic matter for hydrologic solutions. Soil Science
#'    Society of America Journal, 70(5), 1569-1578. <doi:10.2136/sssaj2005.0117>
#'
#' @examples 
#' S = C = 0.20; OM = 2.5 
#' wilting_point(S, C, OM)  # 13.8%
#' field_capacity(S, C, OM) # 32.1% 
#' saturated_mois(S, C, OM) # 48.2%
NULL

#' @rdname SPAW
#' @export
wilting_point <- function(S, C, OM){
    theta_1500t = -0.024*S + 0.487*C + 0.006*OM +
        0.005*(S*OM) - 0.013*(C*OM) +
        0.068*(S*C) + 0.031
    # Eq. 1, theta_1500 (1500 kPa moisture), %v
    theta_1500t + (0.14 * theta_1500t - 0.02) 
}

#' @rdname SPAW
#' @export
field_capacity <- function(S, C, OM){
    theta_33t = -0.251*S + 0.195*C + 0.011*OM +
        0.006*(S*OM) - 0.027*(C*OM) +
        0.452*(S*C) + 0.299
    # Eq. 2, theta_33 (33 kPa moisture), %v
    theta_33t + (1.283 * theta_33t^2 - 0.374 * theta_33t - 0.015) 
}


theta_S_33 <- function(S, C, OM){
    theta_S_33_t = 0.278*S + 0.034*C + 0.022*OM -
        0.018*(S*OM) - 0.027*(C*OM) - 0.584*(S*C) + 0.078
    # Eq. 3, 
    theta_S_33_t + 0.636*theta_S_33_t - 0.107
}

#' @rdname SPAW
#' @export
saturated_mois <- function(S, C, OM){
    theta33 = theta_33(S, C, OM)
    thetaS_33 = theta_S_33(S, C, OM)
    # Eq. 5 (theta_S)
    theta33 + thetaS_33 - 0.097*S + 0.043
}

theta_1500 = wilting_point
theta_33   = field_capacity
theta_S    = saturated_mois

# Eq. 6
rho_norm <- function(S, C, OM){
    (1 - theta_S(S, C, OM))*2.65
}

# Eq. 19
Rv <- function(S, C, OM){
    alpha = (1 - theta_S(S, C, OM))
}

A <- function(S, C, OM){
    theta33   = theta_33(S, C, OM)
    b = B(S, C, OM)
    exp( log(33) + b*log(theta33) )
}

B <- function(S, C, OM){
    theta1500 = theta_1500(S, C, OM)
    theta33   = theta_33(S, C, OM)
    (log(1500) - log(33)) / (log(theta33) - log(theta1500))
}

lambda <- function(S, C, OM) {
    b = B(S, C, OM)
    1/b
}

# eq. 23-24
psi_O <- function(theta, thetaS, EC){
    thetaS/theta*36*EC
}

psi_1500_33 <- function(theta, S, C, OM){
    a = A(S, C, OM)
    b = B(S, C, OM)
    a*theta^(-b)
}

# salinity
#' @importFrom stats uniroot
#' @param EC Electrical conductance of a saturated soil extract (dS/m = mili-mho/cm)
#' @rdname SPAW
wilting_point_salinity <- function(S = 0.2, C = 0.2, OM = 2.5, EC = 3){
    # theta = 14.48/100
    thetaS = saturated_mois(S, C, OM)
    goal <- function(theta) {
        1500 - psi_1500_33(theta, S, C, OM) - psi_O(theta, thetaS, EC)
    }
    uniroot(goal, c(1e-4, 1))$root
}

#' VIC model soil parameters
#' @name VIC_soilParam
#'
#' @description
#' Soil parameters for VIC model (calculated from HWSD database)
#'
#' - `expt` : 3+2b (Eq. 17)
#' - `Kstat`: Saturated conductivity (mm/hr) (Eq. 16)
#' - `bubble`: Bubbling pressure (cm) (Eq. 4)
#' - `Wcr_FT  `: Field Capacity, Fractional soil moisture content at the
#'    critical point (~70% of field capacity) (fraction of maximum moisture)
#' - `Wpwp_FT `: Wilting Point, Fractional soil moisture content at the
#'    wilting point (fraction of maximum moisture)
#'
#' @inheritParams wilting_point
#'
#' @references
#' Saxton, K. E., & Rawls, W. J. (2006). Soil water characteristic estimates by
#'  texture and organic matter for hydrologic solutions. Soil Science Society of
#'  America Journal, 70(5), 1569-1578. https://doi.org/10.2136/sssaj2005.0117
#'
#' @examples
#' S = C = 0.20; OM = 2.5
#' expt(S, C, OM)
#' Ksat(S, C, OM)           # 12.19
#' bubble(S, C, OM)
#' field_capacity(S, C, OM) # 32.1%
#' wilting_point(S, C, OM)  # 13.8%
#' # Fractional soil moisture
#' Wcr_FT(S, C, OM)
#' Wpwp_FT(S, C, OM)
NULL

# Eq. 17,
#' @rdname VIC_soilParam
#' @export
expt <- function(S, C, OM){
    b = B(S, C, OM)
    3 + 2*b
}

#' @details
#' - Kstat: unit has changed from `mm/h` to `mm/day`
#' 
#' @rdname VIC_soilParam
#' @export
Ksat <- function(S, C, OM) {
    theta33 = theta_33(S, C, OM)
    thetaS  = saturated_mois(S, C, OM)
    lambda  = lambda(S, C, OM)

    1930*(thetaS - theta33)^(3 - lambda) / 24
}

#' @rdname VIC_soilParam
#' @export
bubble <- function(S, C, OM){
    thetaS_33 = theta_S_33(S, C, OM)

    psi = -21.67*S - 27.93*C - 81.97*thetaS_33 + 71.12*(S*thetaS_33) +
        8.29*(C*thetaS_33) + 14.05*(S*C) + 27.16

    psi + 0.02*psi^2 - 0.113*psi - 0.70
}

#' @rdname VIC_soilParam
#' @export
Wcr_FT <- function(S, C, OM) {
    full = saturated_mois(S, C, OM)
    field_capacity(S, C, OM)/full
}

#' @rdname VIC_soilParam
#' @export
Wpwp_FT <- function(S, C, OM) {
    full = saturated_mois(S, C, OM)
    wilting_point(S, C, OM)/full
}

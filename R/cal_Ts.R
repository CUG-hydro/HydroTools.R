#' Evaporative surface temperature
#'
#' Evaporative surface temperature (land surface, Tland; or leaf surface Tleaf).
#' Land surface temperature infered by Monteith 1965 Equation.
#'
#' @details
#'- `rH`: aerodynamic resistance of heat
#'- `rs`: stamotal resistance of water
#'
#' @inheritParams ET0_Penman48
#' @inheritParams cal_rH
#' @param rs If `rs = 0`, Monteith 1965 leaf evaporation Equation becomes Penman
#' 1948 water evaporation.
#' ignore the influence of Ts on net cal_radiation
#'
#' @example R/example/ex-cal_Ts.R
#' @export
cal_Ts <- function(Rn, Tair, D, U2, Pa = atm, rH = NULL, rs = 0, ...) {
    # U2 = cal_U2(wind, z.wind)
    rH = cal_rH(U2, h = 0.12)
    rH = cal_rH2(U2, Tair, Pa) # ET_cr推导结果可能会更好

    # gamma_star = gamma * gH / gw
    gamma = cal_gamma(Tair, Pa)
    slope = cal_slope(Tair)

    # gamma_star = gamma * (ra + rs) / rH
    gamma_star = gamma * (1 + rs / rH) # ra ≈ 0.93 rH
    rou_a = 3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5, kg m-3
    # rou_a * Cp * delta_T * gH (in MJ m-2 s-1)
    # = kg m-3 * MJ kg-1 degC-1 * degC * m s-1
    # = MJ m-2 s-1

    # MJ m-2 s-1 * 1e6 = W m-2, then having the same unit as Rn
    # dat_ET = ET0_Monteith65(Rn, Tair, Pa = Pa, D, wind, z.wind, rs = rs)
    # Ts = Tair + dt
    dT = gamma_star / (slope + gamma_star) * ( Rn / (rou_a * Cp / rH * 1e6) - D / gamma_star)
    dT
    # dat_ET$Ts = Ts
    # dat_ET
}

#' wetbulb temperature
#'
#' @inheritParams ET0_Monteith65
#' @inheritParams cal_Rn
#' @export
cal_Tw <- function(ea, Tair, Pa = atm) {
    n <- length(ea)
    ans <- rep(NA_real_, n)
    if (length(Pa) != length(ea) && length(Pa) == 1) Pa <- rep(Pa, n)

    for (i in 1:n) {
        temp <- ea[i] + Tair[i] + Pa[i]
        if (!is.na(temp)) {
            ans[i] <- cal_Tw_default(ea[i], Tair[i], Pa[i])
        }
    }
    ans
}

#' @rdname cal_Tw
#' @export
cal_Tw_default <- function(ea, Tair, Pa = atm) {
    ea = pmin(ea, cal_ea(Tair)) # make sure ea in a reasonable range
    gamma = cal_gamma(Tair, Pa) # lambda changes slightly as Tair changes
    # lambda = cal_lambda(Tair)

    goal <- function(Tw) {
        # rou_a = 1 # ignored
        # f1 = - Cp * rou_a * (Tw - Ta)
        # f2 = lambda * (q_w - q_a) * rou_a
        f1 = cal_es(Tw) - ea
        f2 = -gamma * (Tw - Tair)
        f1 - f2
    }
    uniroot(goal, c(-150,80))$root
}

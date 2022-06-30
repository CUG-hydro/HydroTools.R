#' Evaporative surface temperature
#'
#' Evaporative surface temperature (land surface, Tland; or leaf surface Tleaf).
#' Land surface temperature infered by Monteith 1965 Equation.
#'
#' @details
#' - `rH`: aerodynamic resistance of heat
#' - `rs`: stamotal resistance of water
#'
#' @inheritParams ET0_Penman48
#' @inheritParams cal_rH
#' @param rs If `rs = 0`, Monteith 1965 leaf evaporation Equation becomes Penman
#' 1948 water evaporation.
#' ignore the influence of Ts on net cal_radiation
#' @param method
#' - `simple`: Monteith 1965 Equation
#' - `full` (not finished): Yang 2019
#' - `ma2021`: Ma 2021
#' 
#' @references 
#' 1. Monteith, J. P. (1965). An introduction to environmental physics.
#' 2. Yang, Y., & Roderick, M. L. (2019). Radiation, surface temperature and
#'    evaporation over wet surfaces. Quarterly Journal of the Royal
#'    Meteorological Society, 145(720), 1118–1129.
#'    https://doi.org/10.1002/qj.3481
#' 3. Ma, N., Szilagyi, J., & Zhang, Y. (2021). Calibration-Free Complementary
#'    Relationship Estimates Terrestrial Evapotranspiration Globally. Water
#'    Resources Research, 57(9), 1–27. https://doi.org/10.1029/2021WR029691
#' 
#' @example R/example/ex-cal_Ts.R
#' @export
cal_Ts <- function(Rn, Tair, D, U2, Pa = atm, rH = NULL, rs = 0, 
    method = c("simple", "full", "ma2021"), ...) 
{
    # U2 = cal_U2(wind, z.wind)
    ## ET_cr推导结果可能会更好
    # rH = cal_rH(U2, h = 0.12)  # FAO98
    rH <- cal_rH2(U2, Tair, Pa) # equivalent to Shuttleworth1993

    # gamma_star = gamma * gH / gw
    gamma <- cal_gamma(Tair, Pa)
    slope <- cal_slope(Tair)

    # gamma_star = gamma * (ra + rs) / rH
    gamma_star <- gamma * (1 + rs / rH) # ra ≈ 0.93 rH
    rou_a <- 3.486 * Pa / cal_TvK(Tair) # FAO56, Eq. 3-5, kg m-3
    # rou_a * Cp * delta_T * gH (in MJ m-2 s-1)
    # = kg m-3 * MJ kg-1 degC-1 * degC * m s-1
    # = MJ m-2 s-1

    # MJ m-2 s-1 * 1e6 = W m-2, then having the same unit as Rn
    # Ts = Tair + dt
    method <- match.arg(method)
    if (method == "simple") {
        ## solution1:
        Rln <- 0
        Tw <- gamma_star / (slope + gamma_star) * (
            (Rn - Rln) / (rou_a * Cp / rH * 1e6) - D / gamma_star) + Tair
    } else if (method == "full") {
        ## solution2: considering Rln emitted by the Tw
        goal <- function(Tw) {
            # Campbell and Norman, 1998; 
            # Yang 2018, doi: 10.1002/qj.3481, Eq. A1
            # 这里有写错误，Rln_out重复考了
            emiss = 0.96 
            sigma = 5.67*1e-8 # 5.67*1e-8 Wm−2 K−4
            Rln_out_w <- emiss * sigma * (Tw + K0)^4
            Rln_out_a <- emiss * sigma * (Tair + K0)^4 # Tair
            # Rn = Rn' - Rln
            # delta_Rn = Rln
            ## TODO: bug here, not finished, should use `cal_Rln_yang2019` instead
            Tw_new <- gamma_star / (slope + gamma_star) * (
                (Rn - Rln) / (rou_a * Cp / rH * 1e6) - D / gamma_star) + Tair
            Tw_new - Tw # abs error as the goal function
        }
        Tw <- uniroot(goal, c(-50, 50) + Tair)$root # 不考虑凝结, 则Tw < Tair
    } else if (method == "ma2021") {
        ## solution3: Maning 2021, Eq. 6
        Ep <- ET0_Monteith65(Rn, Tair, Pa = Pa, D, U2, z.wind = 2, rs = rs)$ET0
        # Ep = ET0_Penman48(Rn, Tair, Pa, D, wind = U2, z.wind = 2)$ET0 # Shuttleworth 1993
        # if rs = 0, `ET0_Monteith65` will be degraded to `Shuttleworth 1993`
        beta <- (Rn - Ep) / Ep
        ea <- cal_es(Tair) - D

        goal <- function(Tw) {
            # gamma = cal_gamma(Tw, Pa)
            # gamma_star = gamma * (1 + rs / rH) # ra ≈ 0.93 rH
            f1 <- gamma_star * (Tw - Tair)
            f2 <- beta * (cal_es(Tw) - ea)
            f1 - f2 # abs error as the goal function
        }
        Tw <- uniroot(goal, c(-50, 0) + Tair)$root # 不考虑凝结, 则Tw < Tair
    }
    Tw
    # dat_ET$Ts = Ts
    # dat_ET
}

#' @export
solve_goal <- function(f, goal, interval, ..., tol = 1e-7) {
    func <- function(x) {
        f(x) - goal
    }
    uniroot(func, interval, ..., tol = tol)$root
}


#' wetbulb temperature
#'
#' @inheritParams ET0_Monteith65
#' @inheritParams cal_Rn
#' @export
cal_Tw <- function(ea, Tair, Pa = atm) {
    n <- length(Tair)
    if (length(ea) < n && length(ea) == 1) ea <- rep(ea, n)
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
    ea <- pmin(ea, cal_ea(Tair)) # make sure ea in a reasonable range
    gamma <- cal_gamma(Tair, Pa) # lambda changes slightly as Tair changes
    # lambda = cal_lambda(Tair)

    goal <- function(Tw) {
        # rou_a = 1 # ignored
        # f1 = - Cp * rou_a * (Tw - Ta)
        # f2 = lambda * (q_w - q_a) * rou_a
        f1 <- cal_es(Tw) - ea
        f2 <- -gamma * (Tw - Tair)
        f1 - f2
    }
    uniroot(goal, c(-150, 80))$root
}

wet_bulb <- function(w, Tair, Pa) {
    ea <- w2ea(w, Pa)
    cal_Tw(ea, Tair, Pa)
}

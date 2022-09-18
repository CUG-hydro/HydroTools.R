#' ET complementary
#'
#' @inheritParams ET0_Penman48
#' @param U2 2m wind speed (m/s)
#' @param method method for `Tw` calculation (not used in `ET_CR_Xiao2020`), see 
#' reference for details
#' @param rs stamotal conductance for water on capony scale
#' 
#' @references
#' 1. Xiao, M., Yu, Z., Kong, D., Gu, X., Mammarella, I., Montagnani, L., …
#'    Gioli, B. (2020). Stomatal response to decreased relative humidity
#'    constrains the acceleration of terrestrial evapotranspiration.
#'    Environmental Research Letters, 15(9). doi:10.1088/1748-9326/ab9967
#'
#' 2. Ma, N., Szilagyi, J., & Zhang, Y. (2021). Calibration-Free Complementary
#'    Relationship Estimates Terrestrial Evapotranspiration Globally. Water
#'    Resources Research, 57(9), 1–27. https://doi.org/10.1029/2021WR029691
#'
#' @examples
#' ET_CR_Xiao2020(250, 25, D = 1, U2 = 2)
#' ET_CR_Ma2021(250, 25, D = 1, U2 = 2)
#' @export
ET_CR_Xiao2020 <- function(Rn, Tair, D, U2, Pa = atm, 
    method = c("simple", "full", "ma2021"), rs = 0) 
{
    ea <- cal_es(Tair) - D

    # 湿球温度和蒸发面的温度不同
    # 考虑Rn交互时为Tws，不考虑时为Twb
    # Tw <- cal_Ts(Rn, Tair, D, U2, Pa, method = method, rs = rs) # 蒸发面, 采用的是Penman1948
    Tw <- cal_Tw(ea, Tair, Pa) # wet surface temperature, T_wb

    # b: a coefficient that depicts the proportion of sensible heat flux that increases ETp
    gamma <- cal_gamma(Tair, Pa)
    slope <- (cal_es(Tair) - cal_es(Tw)) / (Tair - Tw)
    b <- slope / gamma # Xiao2020, Eq.3
    beta_wet <- 1 / b # Xiao2020, Eq.9

    ET_p <- ET0_Penman48(Rn, Tair, Pa, D, wind = U2, z.wind = 2)$ET0 # Shuttleworth 1993
    ET_w <- ET0_eq(Rn, Tw, Pa)$Eeq # in the wet surface
    ET_a <- ((1 + b) * ET_w - ET_p) / b
    data.table(Tair, T_wb = Tw, ET_p, ET_w, ET_a)
}

#' @rdname ET_CR_Xiao2020
#' @export
ET_CR_Ma2021 <- function(Rn, Tair, D, U2, Pa = atm, 
    method = c("simple", "full", "ma2021"), rs = 0) 
{
    Ta = Tair
    gamma <- cal_gamma(Tair, Pa)
    ea <- cal_es(Tair) - D

    T_ws <- cal_Ts(Rn, Tair, D, U2, Pa, method = method, rs = rs)
    T_wb <- cal_Tw(ea, Tair, Pa) # wet bulb temperature, Eq. 8
    T_dry <- T_wb + cal_es(T_wb) / gamma # dry air temperature, ea = 0, Eq. 9
    D_dry <- cal_es(T_dry) #

    Ep <- ET0_Penman48(Rn, Tair, Pa, D, wind = U2)$ET0 # Shuttleworth 1993
    Ep_max <- ET0_Penman48(Rn, T_dry, Pa, D_dry, wind = U2)$ET0

    Ew <- ET0_PT72(Rn, T_ws, Pa, alpha = 1.26)$ET0

    X <- (Ep_max - Ep) / (Ep_max - Ew) * (Ew / Ep)
    y <- (2 - X) * X^2
    ET_a <- Ep * y

    ## 检查Ts，是否计算正确，能量是否拟合
    ## TODO: 与马宁老师代码核对
    # # 能量可以完整的闭合, amazing! # kdd, 20220630
    # # 数值实验证明，公式推导完全正确
    # rH = cal_rH(U2, h = 0.12)
    # rH <- cal_rH2(U2, Tair, Pa) # 两种方法计算的rH差别还挺大
    # gH = 1/rH
    # gw = 1/(rH + rs)

    # lambda = cal_lambda(Tair)
    # slope = cal_slope(Tair)
    # # Rn2 = Rn * 86400 / 1e6 # MJ/d/m2 -> W/m2
    # rho_a <- 3.486 * Pa / cal_TvK(Tair)

    # H = (T_ws - Tair) * rho_a * Cp * gH * 1e6
    # d_es1 = cal_es(T_ws) - ea
    # d_es2 = slope * (T_ws - Ta) + D
    # E1 = lambda * epsilon * rho_a * gw * (cal_es(T_ws) - ea) / Pa * 1e6
    # E2 = lambda * epsilon * rho_a * gw / Pa * 1e6 * 
    #     (slope * (T_ws - Ta) + D)
    # tmp = data.table(Rn, H = H, E1, E2)
    # print(tmp)
    data.table(Tair, T_wb, T_ws, T_dry,
        ET_p = Ep, Ep_max, ET_w = Ew, X, y, ET_a
    )
}

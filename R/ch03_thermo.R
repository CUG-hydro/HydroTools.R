#! 这一部分气压用的hPa，为符合大气的习惯

#' Lifted Condensation Level (LCL)
#' 
#' @description 
#' - `LCL`   : Bolton 1980, Eq. 15
#' - `LCL_RH`: Bolton 1980, Eq. 22
#' 
#' @param P0 pressure at surface (hPa)
#' @param T0 temperature at surface (C)
#' @param Td dew point temperature (C)
#' 
#' @export
LCL <- function(P0, T0, Td) {
  ea <- cal_es(Td) * 10
  w <- ea2w(ea, P0) # g / g
  goal <- function(P) {
    T <- adiabat_T_dry(P0, T0, P) - K0 # in Kdeg
    es <- cal_es(T) * 10
    ws <- ea2w(es, P)
    w - ws
  }

  P_lcl <- uniroot(goal, c(20, P0))$root # 1000 hPa to 20 hPa
  T_lcl <- adiabat_T_dry(P0, T0, P_lcl)
  c(P_lcl = P_lcl, T_lcl = T_lcl - K0)
}

#' @rdname LCL
#' @export
LCL_RH <- function(T, RH) {
  TK = T + K0
  term2 = 1 / (TK - 55) - log(RH/100) / 2840 
  1 / term2 + 55 - K0
}

#' @export
adiabat_T_dry <- function(P0, T0, P, w=NULL) {
  if (is.null(w)) {
    m <- Rd / (Cp * 1e6) # 0.283
  } else {
    m <- 0.2854 * (1 - 0.28 * w) # Bolton 1980
  }
  (T0 + K0) * (P / P0)^m
}

#' @rdname theta_wet
#' @export
theta <- function(P0, T0) {
  adiabat_T_dry(P0, T0, 1000) # hPa
}

#' equivalent potential temperature
#' 
#' @param P0 pressure at surface (hPa)
#' @param T0 temperature at surface (Cdeg)
#' @param Td dew point temperature (Cdeg)
#' 
#' @references 
#' 1. https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.equivalent_potential_temperature.html
#' 2. https://github.com/wrf-model/WRF/blob/master/phys/module_diag_functions.F#L122
#' 3. https://github.com/Unidata/MetPy/blob/e0e24d51702787943fc3c0481fa9a6632abe9d20/src/metpy/calc/thermo.py#L1512
#' 
#' @examples 
#' theta_wet(850, 20, 18)
#' theta_wet_bolton(850, 20, 18)
#' 
#' @export
theta_wet <- function(P0, T0, Td) {
  T_lcl = LCL(P0, T0, Td)["T_lcl"]
  K_lcl = T_lcl + K0
  
  Lv = cal_lambda(T_lcl) #* 1e6
  # Lv = 2.5 * 1e6 MJ kg-1
  # Cp =  # MJ kg-1 degC-1
  w = Tdew2w(Td, P0/10) # w守恒
  theta_d = adiabat_T_dry(P0, T0, 1000, w) # hPa
  
  theta_se = theta_d * exp(Lv * w / (Cp * K_lcl))
  c(T_lcl, theta_d = theta_d - K0, theta_se = theta_se[[1]] - K0)
}

#' @rdname theta_wet
#' @export
theta_wet_bolton <- function(P0, T0, Td) {
  ea <- cal_es(Td) * 10

  td <- Td + K0
  tk <- T0 + K0

  r <- ea2w(ea, P0)
  t_l <- 56 + 1. / (1. / (td - 56) + log(tk / td) / 800.)
  th_l <- theta(P0 - ea, T0) * (tk / t_l)^(0.28 * r)

  theta_se <- th_l * exp(r * (1 + 0.448 * r) * (3036. / t_l - 1.78))
  c(T_lcl = t_l - K0, theta_se = theta_se - K0)
}


# Pdry_adiabat <- function(T0, P0, T) {
#   m <- Rd / (Cp * 1e6) # 0.283
#   ((T + K0) / (T0 + K0))^(1 / m) * P0
#   # (T/T0)^1/m * P0
#   # (T0) * (P / P0)^m
#   # (T0 + K0) * (P / P0)^m - K0
# }

# Pwet_adiabat <- function(T0, P0, w0, T) {
#   P = Pdry_adiabat(T0, P0, T)

#   ea = w2ea(w0, P0) # kPa
#   print(ea)

#   cal_Tw(ea, T, P)
#   # wet_bulb()
#   # T_wet = cal_Tw(ea, T) # Pa also known
#   # m <- Rd / (Cp * 1e6) # 0.283
#   # ((T + K0) / (T0 + K0))^(1 / m) * P0
#   # (T/T0)^1/m * P0
#   # (T0) * (P / P0)^m
#   # (T0 + K0) * (P / P0)^m - K0
# }

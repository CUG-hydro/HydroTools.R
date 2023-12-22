#! 这一部分气压用的hPa，为符合大气的习惯

#' Lifted Condensation Level (LCL)
#' 
#' @param T0 temperature at surface (K)
#' @param P0 pressure at surface (hPa)
#' @param Td dew point temperature (K)
#' 
#' @export
LCL <- function(P0, T0, Td) {
  ea <- cal_es(Td) * 10
  
  goal <- function(P) {
    T <- Tdry_adiabat(P0, T0, P) - K0 # in Kdeg
    es <- cal_es(T) * 10
    es - ea
  }

  # 1000 hPa to 20 hPa
  P_lcl <- uniroot(goal, c(20, P0))$root

  T_lcl <- Tdry_adiabat(P0, T0, P_lcl)
  c(P_lcl = P_lcl, T_lcl = T_lcl - K0)
}

#' theta_wet2
#' 
#' @references 
#' 1. https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.equivalent_potential_temperature.html
#' @examples 
#' theta_wet2(850, 20, 18)
#' 
#' @export 
theta_wet2 <- function(P0, T0, Td) {
  ea <- cal_es(Td) * 10
  
  td = Td + K0
  tk = T0 + K0

  r <- ea2w(ea, P0)
  t_l <- 56 + 1. / (1. / (td - 56) + log(tk / td) / 800.)
  th_l <- theta(P0 - ea, T0) * (tk / t_l)^(0.28 * r)
  
  theta_se = th_l * exp(r * (1 + 0.448 * r) * (3036. / t_l - 1.78))

  c(T_lcl = t_l - K0, theta_se = theta_se - K0)
}
    

#' @export
Tdry_adiabat <- function(P0, T0, P) {
  m <- Rd / (Cp * 1e6) # 0.283
  (T0 + K0) * (P / P0)^m
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

#' @export
theta <- function(P0, T0) {
  Tdry_adiabat(P0, T0, 1000) # hPa
}

#' @export
theta_wet <- function(P0, T0, Td) {
  theta_d = Tdry_adiabat(P0, T0, 1000) # hPa
  
  T_lcl = LCL(P0, T0, Td)["T_lcl"]
  Lv = cal_lambda(T0) #* 1e6
  # Lv = 2.5 * 1e6 MJ kg-1
  w = Tdew2w(Td, P0) # 

  theta_d * exp(Lv * w / Cp * T_lcl)
}

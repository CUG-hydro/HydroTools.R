# Lifted Condensation Level (LCL)
#' @export
LCL <- function(T0, P0, Td) {
  ea <- cal_es(Td)
  goal <- function(P) {
    es <- Tdry_adiabat(T0, P0, P) %>% cal_es()
    es - ea
  }
  # 1000 hPa to 20 hPa
  P_lcl <- uniroot(goal, c(20, 1000)/10)$root
  T_lcl <- Tdry_adiabat(T0, P0, P_lcl)
  c(P_lcl = P_lcl, T_lcl = T_lcl)
}

#' @export
Tdry_adiabat <- function(T0, P0, P) {
  m <- Rd / (Cp * 1e6) # 0.283
  # (T0) * (P / P0)^m 
  (T0 + K0) * (P / P0)^m - K0
  # T0 * (P / P0)^m
}

Pdry_adiabat <- function(T0, P0, T) {
  m <- Rd / (Cp * 1e6) # 0.283
  ((T + K0) / (T0 + K0))^(1 / m) * P0
  # (T/T0)^1/m * P0
  # (T0) * (P / P0)^m 
  # (T0 + K0) * (P / P0)^m - K0
}

Pwet_adiabat <- function(T0, P0, w0, T) {
  P = Pdry_adiabat(T0, P0, T)
  
  ea = w2ea(w0, P0) # kPa
  print(ea)

  cal_Tw(ea, T, P)
  # wet_bulb()
  # T_wet = cal_Tw(ea, T) # Pa also known
  # m <- Rd / (Cp * 1e6) # 0.283
  # ((T + K0) / (T0 + K0))^(1 / m) * P0
  # (T/T0)^1/m * P0
  # (T0) * (P / P0)^m 
  # (T0 + K0) * (P / P0)^m - K0
}

#' @export
theta <- function(T0, P0) {
  dry_adiabat(T0, P0, 100) # kPa
}

#' @export
theta_wet <- function(T0, P0, Td) {
  theta_d = dry_adiabat(T0, P0, 100) # hPa
  
  T_lcl = LCL(T0, PO, Td)["T_lcl"]
  Lv = cal_lambda(T0) #* 1e6
  # Lv = 2.5 * 1e6 MJ kg-1
  w = Tdew2w(Td, P0) # 

  theta_d * exp(Lv * w / Cp * T_lcl)
}

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
  if (length(ea) < n && length(ea) == 1) ea = rep(ea, n)
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
  ea = w2ea(w, Pa)
  cal_Tw(ea, Tair, Pa)
}

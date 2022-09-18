#' black body radiation
#' 
#' @param Ta temperature in C deg
#' @param emiss Emissivity coefficient
#' 
#' @return radiation in W m-2
#' @export
cal_Rl_out <- function(Ta, emiss = 0.97) {
  sigma = 5.67 * 1e-8 # W m-2 K-4
  emiss * sigma * (Ta + K0)^ 4 # W m-2
}

#' @rdname cal_radiation_Ta
#' @export 
radiation_Ta <- function(Rs, Rl = 0, albedo = 0.3, emiss = 0.97) {
  sigma = 5.67 * 1e-8 # W m-2 K-4
  Rl_out = (Rs* (1 - albedo) + emiss * Rl ) 
  (Rl_out / (emiss * sigma))^(1/4) - K0
}

#' air temperature under solar radiation
#' 
#' @param Ta scalar, the initial air temperature, in `degC`
#' @param Rs,Rl numeric vector, inward shortwave and longwave radiation, in `W/m2`
#' @param dt delta `t`, in second
#' 
#' @export 
cal_radiation_Ta <- function(Ta, Rs, Rl = 0, dt = 3600, method = c("approx", "exact"), ..., 
  albedo = 0.3, emiss = 0.97) {
  method = match.arg(method)

  FUN = switch(method, 
    "approx" = .cal_radiation_Ta_approx, 
    "exact" = .cal_radiation_Ta_exact)
  n = length(Rs)
  if (length(Rl) == 1) Rl = rep(Rl, n)

  Ta = rep(Ta, n)
  for (i in 2:n) {
    Ta[i] = FUN(Ta[i-1], Rs[i], Rl[i], dt, albedo, emiss)
  }
  Ta
}

# 粗略解
.cal_radiation_Ta_approx <- function(Ta, Rs, Rl = 0, dt = 3600, 
  albedo = 0.3, emiss = 0.97) {
  
  Cp    = 1.013 * 1e3   # J kg-1 degC-1
  # Ta    = Ta + K0       # to kdeg
  rho_a = cal_rho_a(Ta) # kg m-3

  # Rn = (Rs * (1 - albedo) + emiss * (Rl - sigma * Ta^4) )
  Rl_out = cal_Rl_out(Ta) # W m-2
  Rn = Rs * (1 - albedo) + emiss * (Rl - Rl_out) 
  Tnew = Ta + Rn * dt / (rho_a * Cp)
  tibble(Ta, Rn, Rn * dt / (rho_a * Cp), Rl_out, Tnew) %>% print()
  Tnew - K0 # to Cdeg
}

# 精确解
.cal_radiation_Ta_exact <- function(Ta, Rs, Rl = 0, dt = 3600, 
  albedo = 0.3, emiss = 0.97) {
  
  # Ta    = Ta + K0       # to kdeg
  Cp    = 1.013 * 1e3   # J kg-1 degC-1
  rho_a = cal_rho_a(Ta) # km m-3
  sigma = 5.67 * 1e-8 # W m-2 K-4
  goal <- function(Tnew) {
    Rl_out = emiss * sigma * (Tnew^5 - Ta^5) / (5*(Tnew - Ta)) # W m-2
    Rn = Rs * (1 - albedo) + emiss * (Rl - Rl_out) 
    Tnew2 = Ta + Rn * dt / (rho_a * Cp)
    Tnew2 - Tnew
  }
  uniroot(goal, c(-50, 50) + Ta)$root - K0
}

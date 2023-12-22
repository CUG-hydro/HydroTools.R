# https://github.com/Sibada/sibadaR/blob/master/R/FAO56.R

#' @title Helper functions for vapour pressure
#' @name vapour_press 
#' 
#' @param q specific humidity in kg/kg or g/g
#' @param Pa surface air pressure
#' 
#' @return e in the same unit as Pa
#' 
#' @references
#' https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html, Eq-17
#'
#' @examples
#' RH = 90
#' Pa = atm
#' Tair = 30
#' q = RH2q(RH, Pa, Tair)
#' RH2 = q2RH(q, Pa, Tair)
#' e = q2ea(q, Pa)
#' es = cal_es(Tair)
NULL


#' @param w mix ratio, m_w / m_d
#' @rdname vapour_press
#' @export
w2q <- function(w) w / (w + 1)

#' @rdname vapour_press
#' @export
q2w <- function(q, Pa = atm) q / (1 - q)


#' @rdname vapour_press
#' @export
q2ea <- function(q, Pa = atm) {
  q * Pa / (epsilon + (1 - epsilon) * q)
}

#' @rdname vapour_press
#' @export
w2ea <- function(w, Pa = atm) {
  q = w2q(w)
  q2ea(q, Pa)
}

#' @rdname ET0_helper
#' @export
cal_es <- function(Tair) {
    # FAO98 equation
    0.6108 * exp((17.27 * Tair) / (Tair + 237.3))
    # 0.61094 * exp((17.625 * Tair) / (Tair + 243.04))
}

es2T <- function(es) {
  es = es * 10 # to hPa
  (243.5 * log(es) - 440.8) / (19.48 - log(es))
}

# https://github.com/rpkgs/skew-t/blob/master/thermo_scripts.py#L42
#' @export
cal_es_CC <- function(Tair) {
    lv = 2.5 * 1e6
    # lv = cal_lambda(Tair) * 1e6
    # Rv = 461.5
    Tair = Tair + K0
    0.6108 * exp( (lv / Rw) * ((1 / 273.15) - (1 / Tair)) )
    # 0.6108 * exp((17.27 * Tair) / (Tair + 237.3))
}

#' @param ea actual vapor pressure (kPa)
#' 
#' @rdname vapour_press
#' @export
ea2VPD <- function(ea, RH) {
  ea/(RH/100) - ea
}

#' @rdname vapour_press
#' @export
vapour_press <- q2ea

#' @param Tair air temperature, in degC
#' @rdname vapour_press
#' @export
q2RH <- function(q, Tair, Pa = atm) {
  ea <- vapour_press(q, Pa)
  es <- cal_es(Tair)
  ea / es * 100
}

#' @param RH relative humidity, in %
#' @references
#' https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
#' 
#' @rdname vapour_press
#' @export
RH2q <- function(RH, Tair, Pa = atm) {
  # ea <- vapour_press(q, Pa)
  es <- cal_es(Tair)
  ea <- es * RH/100

  # ws = epsilon * es / (Pa - es)
  w <- epsilon * ea / (Pa - ea) # eq. 1.55
  # w = RH/100 * ws 
  w2q(w) # q: g / g
}

#' @rdname vapour_press
#' @export
q_from_RH = RH2q

#' @rdname vapour_press
#' @export
cal_qs <- function(Tair, Pa = atm) {
  es <- cal_es(Tair)
  ea2q(es, Tair, Pa) # q
}

#' @rdname vapour_press
#' @export
ea2w <- function(ea, Pa = atm) {
  # es <- cal_es(Tair)
  w <- epsilon * ea / (Pa - ea)
  w # g / g
}

#' @rdname vapour_press
#' @export
ea2q <- function(ea, Pa = atm) {
  # es <- cal_es(Tair)
  w <- epsilon * ea / (Pa - ea)
  w2q(w) # q: g / g
}

#' @rdname vapour_press
#' @export
 RH2q <- q_from_RH

#' @param Tdew dew temperature (in degC)
#' 
#' @rdname vapour_press
#' @export
Tdew2q <- function(Tdew, Pa = atm) {
  ea <- cal_es(Tdew)
  # es <- cal_es(Tair)
  epsilon * ea / (Pa - (1 - epsilon) * ea)
}

#' @rdname vapour_press
#' @export
Tdew2w <- function(Tdew, Pa = atm) {
  ea <- cal_es(Tdew)
  # es <- cal_es(Tair)
  epsilon * ea / (Pa - ea)
}

#' @rdname vapour_press
#' @export
Tdew2RH <- function(Tdew, Tair) {
  ea <- cal_es(Tdew)
  es <- cal_es(Tair)
  ea / es * 100
}

# Tdew2q(Tdew_from_q(0.1, 100), 100)

#' @rdname vapour_press
#' @export
Tdew_from_q <- function(q, Pa = atm) {
  ea = q2ea(q, Pa)
  solve_goal(cal_es, ea, c(-100, 100))
}

#' @rdname vapour_press
#' @export
Tdew_from_w <- function(w, Pa = atm) {
  ea = w2ea(w, Pa)
  # solve: cal_es(Tair) ≈ ea
  solve_goal(cal_es, ea, c(-100, 100))
}

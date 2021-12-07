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

#' @rdname vapour_press
#' @export
q2ea <- function(q, Pa) {
  q * Pa / (epsilon + (1 - epsilon) * q)
}


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
q2RH <- function(q, Pa, Tair) {
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
RH2q <- function(RH, Pa, Tair) {
  # ea <- vapour_press(q, Pa)
  es <- cal_es(Tair)
  ea <- es * RH/100

  # ws = epsilon * es / (Pa - es)
  w <- epsilon * ea / (Pa - ea)
  # w = RH/100 * ws 
  w2q(w) # q: g / g
}


#' @param w mix ratio, m_w / m_d
#' @rdname vapour_press
#' @export
w2q <- function(w, Pa) w / (w + 1)

#' @rdname vapour_press
#' @export
q2w <- function(q, Pa) q / (1 - q)

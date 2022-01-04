
#' @title constants in meteorology 
#' @name constants
#'
#' @description
#' - `T0`: [Absolute zero](https://en.wikipedia.org/wiki/Absolute_zero) in
#'   Kelvin \eqn{T_0} (K)
#' 
#' - `Es.T0`: \eqn{e_s(T_0) = 6.11hPa} is the saturation vapor pressure at the
#'   absolute zero \eqn{T_0 = 273.15K}.
#'
#' - `L`: [Latent heat](https://en.wikipedia.org/wiki/Latent_heat) of water
#'   vapor \eqn{L = 2.5 \times 10^6J/kg}
#' 
#' - `atm`: standard pressure at sea surface, `101.325 kPa`
#' 
#' - `Cp`: specific heat at constant pressure, `1.013 * 1e-3 MJ kg-1 degC-1`
NULL

#' @rdname constants
#' @export
T0 <- 273.15

#' @rdname constants
#' @export
L <- 2.5e6

#' @rdname constants
#' @export
Es.T0 <- 6.11

#' @export
#' @rdname constants
atm = 101.325

#' Molecular weight
#' 
#' @description
#' `Mw`: Molecular weight of dry air \eqn{M_d = 28.9634g/mol}
#' `Md`: [Molecular weight](https://en.wikipedia.org/wiki/Molar_mass) of water vapor \eqn{M_w = 18.01528g/mol}
#' `epsilon`: the ratio of `Mw` to `Md`
#' @export
Mw <- 18.01528

#' @rdname Mw
#' @export
Md <- 28.9634

#' @rdname Mw
#' @export
epsilon <- Mw / Md

#' Specific gas constants
#' 
#' @description
#' - `R`: [gas constant](https://en.wikipedia.org/wiki/Gas_constant), J/(mol K)
#' - `Rw`: [Specific gas constant](https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant)
#' of water vapor \eqn{R_w = \frac{1000R}{M_w} = 461.52J/(kg K)}.
#' - `Rd`: Specific gas constant of dry air.
#' @export
R <- 8.3144621 # J/(mol K)

#' @rdname R
#' @export
Rw <- R / Mw * 1000

#' @rdname R
#' @export
Rd <- R / Md * 1000

#' @rdname R
#' @export
Cp = 1.013 * 1e-3 # MJ kg-1 degC-1

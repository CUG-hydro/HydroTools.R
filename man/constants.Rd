% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constant.R
\docType{data}
\name{constants}
\alias{constants}
\alias{T0}
\alias{K0}
\alias{L}
\alias{Es.T0}
\alias{atm}
\title{constants in meteorology}
\format{
An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.
}
\usage{
T0

K0

L

Es.T0

atm
}
\description{
\itemize{
\item \code{T0}: \href{https://en.wikipedia.org/wiki/Absolute_zero}{Absolute zero} in
Kelvin \eqn{T_0} (K)
\item \code{Es.T0}: \eqn{e_s(T_0) = 6.11hPa} is the saturation vapor pressure at the
absolute zero \eqn{T_0 = 273.15K}.
\item \code{L}: \href{https://en.wikipedia.org/wiki/Latent_heat}{Latent heat} of water
vapor \eqn{L = 2.5 \times 10^6J/kg}
\item \code{atm}: standard pressure at sea surface, \verb{101.325 kPa}
\item \code{Cp}: specific heat at constant pressure, \verb{1.013 * 1e-3 MJ kg-1 degC-1}
}
}
\keyword{datasets}

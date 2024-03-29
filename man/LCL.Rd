% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ch03_thermo.R
\name{LCL}
\alias{LCL}
\alias{LCL_RH}
\title{Lifted Condensation Level (LCL)}
\usage{
LCL(P0, T0, Td)

LCL_RH(T, RH)
}
\arguments{
\item{P0}{pressure at surface (hPa)}

\item{T0}{temperature at surface (C)}

\item{Td}{dew point temperature (C)}
}
\description{
\itemize{
\item \code{LCL}   : Bolton 1980, Eq. 15
\item \code{LCL_RH}: Bolton 1980, Eq. 22
}
}

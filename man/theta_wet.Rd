% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ch03_thermo.R
\name{theta}
\alias{theta}
\alias{theta_wet}
\alias{theta_wet_bolton}
\title{equivalent potential temperature}
\usage{
theta(P0, T0)

theta_wet(P0, T0, Td)

theta_wet_bolton(P0, T0, Td)
}
\arguments{
\item{P0}{pressure at surface (hPa)}

\item{T0}{temperature at surface (Cdeg)}

\item{Td}{dew point temperature (Cdeg)}
}
\description{
equivalent potential temperature
}
\examples{
theta_wet(850, 20, 18)
theta_wet_bolton(850, 20, 18)

}
\references{
\enumerate{
\item https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.equivalent_potential_temperature.html
\item https://github.com/wrf-model/WRF/blob/master/phys/module_diag_functions.F#L122
\item https://github.com/Unidata/MetPy/blob/e0e24d51702787943fc3c0481fa9a6632abe9d20/src/metpy/calc/thermo.py#L1512
}
}

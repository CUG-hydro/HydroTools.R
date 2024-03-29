% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ET_budyko.R
\name{budyko}
\alias{budyko}
\alias{ET_budyko}
\alias{budyko_goal}
\alias{budyko_fit}
\title{Budyko Curve}
\usage{
budyko(PET.P, par = 2, method = "Fu1981")

ET_budyko(PET, P, par = 2, method = "Fu1981", ...)

budyko_goal(AET.P, PET.P, par = 2, method = "Fu1981", ...)

budyko_fit(AET.P, PET.P, par = 2, method = "Fu1981", ...)
}
\arguments{
\item{PET.P}{PET divide P,
\itemize{
\item \code{PET}: potential evapotranspiration, can be calculated by \code{\link[=ET0_Penman48]{ET0_Penman48()}}.
\item \code{P}: Precipitation
\item \code{AET.P}: \code{AET / P}
\item \code{PET.P}: \code{PET / P}
}}

\item{par}{scalar value, the parameter of Budyko curve}

\item{method}{one of \code{c("budyko", "Fu1981", "Zhang2001", "Yang2008", "Pike1964")},
or \code{all}.}

\item{PET}{potential evapotranspiration (mm)}

\item{P}{precipitation (mm)}

\item{...}{ignored}

\item{AET.P}{\code{AET / P}}
}
\description{
Budyko Curve
}
\references{
\enumerate{
\item Budyko, M. I., Climate and Life, 508 pp., Academic, San Diego, Calif., 1974.
\item Fu, Baw-Puh., 1981. On the Calculation of the Evaporation from Land
Surface. Chinese Journal of Atmospheric Sciences, 5(1), 23-31.
\item Zhang, L., Dawes, W. R., & Walker, G. R. (2001). Response of mean annual
evapotranspiration to vegetation changes at catchment scale. Water
Resources Research, 37(3), 701–708. \url{doi:10.1029/2000WR900325}
\item Yang, H., Yang, D., Lei, Z., & Sun, F. (2008). New analytical derivation
of the mean annual water-energy balance equation. Water Resources
Research. \url{doi:10.1029/2007WR006135}
}
}

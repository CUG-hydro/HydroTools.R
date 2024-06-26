% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calib_XAJ.R
\name{XAJ_calib}
\alias{XAJ_calib}
\title{XAJ model Parameter calibration}
\usage{
XAJ_calib(
  Qobs,
  prcp,
  ET0,
  area,
  dt = 24,
  date = NULL,
  maxn = 1000,
  index = "KGE",
  seed = 1,
  ...
)
}
\arguments{
\item{Qobs}{Observed Total runoff, (m^3/s)}

\item{prcp}{Precipitation (mm/d)}

\item{ET0}{Pan evaporation or potential evapotranspiration (mm/d)}

\item{area}{basin area (km^2).}

\item{dt}{time step (hour)}

\item{date}{(optional) corresponding date of \code{Qobs}}

\item{index}{KGE or NSE}

\item{seed}{(can be ignored) starting number of random number generator,
see \code{\link[base:Random]{base::set.seed()}} for details.
This parameter is to make sure optimization result is same in different tries.}

\item{...}{ignored}
}
\description{
XAJ model Parameter calibration
}

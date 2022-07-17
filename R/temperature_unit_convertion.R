#' Convert from one temperature metric to another
#'
#' This function allows you to convert a vector of temperature values between
#' Fahrenheit, Celsius, and degrees Kelvin.
#'
#' @param temperature A numeric vector of temperatures to be converted.
#' @param old_metric The metric from which you want to convert. Possible options are:
#' - `fahrenheit`, `f`
#' - `kelvin`, `k`
#' - `celsius`, `c`
#' @param new_metric The metric to which you want to convert. The same options
#'    are possible as for `old_metric`.
#'
#' @author
#' Joshua Ferreri \email{joshua.m.ferreri@@gmail.com},
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#' @return A numeric vector with temperature converted to the metric specified
#'    by the argument `new_metric`.
#'
#' @export
convert_temperature <- function(temperature, old_metric, new_metric) {
  possible_metrics <- c("fahrenheit", "celsius", "kelvin", "f", "c", "k")
  if (!(old_metric %in% possible_metrics) ||
    !(new_metric %in% possible_metrics)) {
    stop(paste0(
      "The arguments `old_metric` and `new_metric` can only ",
      "have one of the following values: `",
      paste(possible_metrics, collapse = "`, `"), "`"
    ))
  } else if (old_metric == new_metric) {
    stop("`old_metric` and `new_metric` must have different values.")
  }

  if (old_metric %in% c("fahrenheit", "f")) {
    if (new_metric %in% c("celsius", "c")) {
      func <- fahrenheit.to.celsius
    } else if (new_metric %in% c("kelvin", "k")) {
      func <- fahrenheit.to.kelvin
    }
  } else if (old_metric %in% c("celsius", "c")) {
    if (new_metric %in% c("fahrenheit", "f")) {
      func <- celsius.to.fahrenheit
    } else if (new_metric %in% c("kelvin", "k")) {
      func <- celsius.to.kelvin
    }
  } else { # Kelvin for old_metric
    if (new_metric %in% c("fahrenheit", "f")) {
      func <- kelvin.to.fahrenheit
    } else if (new_metric %in% c("celsius", "c")) {
      func <- kelvin.to.celsius
    }
  }

  out <- func(temperature)
  return(out)
}

#' @param T.kelvin,T.fahrenheit,T.celsius Numeric vector of temperatures in Kelvin, Fahrenheit and Celsius.
#' @return A numeric vector of temperature values in the new unit.
#' @references
#' 1. [online heat index calculator](http://www.wpc.ncep.noaa.gov/html/heatindex.shtml)
#' 2. [online temperature converter](http://www.srh.noaa.gov/epz/?n=wxcalc_tempconvert)
#' @rdname convert_temperature
#' @export
celsius.to.fahrenheit <- function(T.celsius) {
  T.fahrenheit <- (9 / 5) * T.celsius + 32
  # T.fahrenheit <- round(T.fahrenheit, digits = round)
  return(T.fahrenheit)
}

#' @rdname convert_temperature
#' @export
celsius.to.kelvin <- function(T.celsius) {
  T.kelvin <- T.celsius + 273.15
  return(T.kelvin)
}

#' @rdname convert_temperature
#' @export
fahrenheit.to.celsius <- function(T.fahrenheit) {
  T.celsius <- (5 / 9) * (T.fahrenheit - 32)
  return(T.celsius)
}

#' @rdname convert_temperature
#' @export
fahrenheit.to.kelvin <- function(T.fahrenheit) {
  T.kelvin <- (.5556 * (T.fahrenheit - 32)) + 273.15
  return(T.kelvin)
}


#' @rdname convert_temperature
#' @export
kelvin.to.celsius <- function(T.kelvin) {
  T.celsius <- T.kelvin - 273.15
  return(T.celsius)
}

#' @rdname convert_temperature
#' @export
kelvin.to.fahrenheit <- function(T.kelvin) {
  T.fahrenheit <- (1.8 * (T.kelvin - 273.15)) + 32
  return(T.fahrenheit)
}

#' Algorithm for heat.index function.
#'
#' `heat_index` converts a numeric scalar of temperature
#'    (in Fahrenheit) and a numeric scalar of relative humidity (in %)
#'    to heat index (in Fahrenheit). This function is not meant to be used
#'    outside of the [heat.index()] function.
#'
#' @param t Numeric scalar of air temperature, in celsius.
#' @param rh Numeric scalar of relative humidity, in %.
#'
#' @details If an impossible value of relative humidity is given
#'    (below 0% or above 100%), heat index is returned as `NA`.
#'
#' @return A numeric scalar of heat index, in celsius
#' @note Equations are from the source code for the US National Weather
#'     Service's
#'     [online heat index calculator](http://www.wpc.ncep.noaa.gov/html/heatindex.shtml).
#' 
#' @author
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#' 
#' @references
#' 1. Anderson GB, Bell ML, Peng RD. 2013. Methods to calculate the heat index
#'    as an exposure Metric in environmental health research.
#'    Environmental Health Perspectives 121(10):1111-1119.
#' 2. National Weather Service Hydrometeorological Prediction
#'    Center Web Team. Heat Index Calculator. 30 Jan 2015.
#'    <http://www.wpc.ncep.noaa.gov/html/heatindex.shtml>.
#'    Accessed 18 Dec 2015.
#' 3. Rothfusz L. 1990. The heat index (or, more than you ever wanted to know
#'    about heat index) (Technical Attachment SR 90-23). Fort Worth:
#'    Scientific Services Division, National Weather Service.
#' 4. R. Steadman, 1979. The assessment of sultriness. Part I: A
#'    temperature-humidity index based on human physiology and clothing
#'    science. Journal of Applied Meteorology, 18(7):861--873.
#' @example R/example/ex-HI.R
#' @export
heat_index <- function (t, rh) {
    t <- celsius.to.fahrenheit(t)
    hi <- unlist(.mapply(.heat.index, dots=list(t, rh), MoreArgs=NULL)) %>%
        fahrenheit.to.celsius()
    hi
}

#' @rdname heat_index
#' @export
heat_index_julia <- function(t, rh) {
    dim = dim(t)
    julia_init()
    julia_call("heat_index",
        t %>% set_dim(NULL),
        rh %>% set_dim(NULL)
    ) %>% set_dim(dim)
}

#' @rdname heat_index
#' @export
.heat.index <- function (t = NA, rh = NA) {
    if (NA %in% c(t, rh)) return(NA)
    if (t <= 40) return(t)

    hi <- -10.3 + 0.047 * rh + 1.1 * t
    if (hi > 79) {
        hi <- -42.379 + 2.04901523 * t + 10.14333127 * rh -
            0.22475541 * t * rh - 6.83783 * 10^-3 * t^2 -
            5.481717 * 10^-2 * rh^2 + 1.22874 * 10^-3 * t^2 *
            rh + 8.5282 * 10^-4 * t * rh^2 - 1.99 * 10^-6 *
            t^2 * rh^2
        if (rh <= 13 & t >= 80 & t <= 112) {
            adj1 <- (13 - rh)/4
            adj2 <- sqrt((17 - abs(t - 95))/17)
            hi <- hi - adj1 * adj2
        } else if (rh > 85 & t >= 80 & t <= 87) {
            adj1 <- (rh - 85)/10
            adj2 <- (87 - t)/5
            hi <- hi + adj1 * adj2
        }
    }
    hi
}

#' @rdname heat_index
#' @export
heat_index_vec <- function(t = NA, rh = NA) {
    dim = dim(t)
    t = celsius.to.fahrenheit(t)
    dim(t) <- NULL
    dim(rh) <- NULL

    hi <- -10.3 + 0.047 * rh + 1.1 * t
    con = which(t <= 40)
    hi[con] = t[con]

    # hi[hi > 79] = hi2[hi > 79]
    hi2 <- -42.379 + 2.04901523 * t + 10.14333127 * rh -
        0.22475541 * t * rh - 6.83783 * 10^-3 * t^2 -
        5.481717 * 10^-2 * rh^2 + 1.22874 * 10^-3 * t^2 *
            rh + 8.5282 * 10^-4 * t * rh^2 - 1.99 * 10^-6 *
            t^2 * rh^2

    # rh <= 13 & t >= 80 & t <= 112
    adj1 <- (13 - rh) / 4
    suppressWarnings(adj2 <- sqrt((17 - abs(t - 95)) / 17))
    delta_1 = -adj1 * adj2
    # rh > 85 & t >= 80 & t <= 87
    adj1 <- (rh - 85) / 10
    adj2 <- (87 - t) / 5
    delta_2 = adj1 * adj2

    con1 = which(rh <= 13 & t >= 80 & t <= 112)
    con2 = which(rh > 85 & t >= 80 & t <= 87)

    hi2[con1] %<>% add(delta_1[con1])
    hi2[con2] %<>% add(delta_2[con2])

    # 37.00379
    con = which(hi > 79)
    hi[con] = hi2[con]
    hi %>% fahrenheit.to.celsius() %>% set_dim(dim)
}

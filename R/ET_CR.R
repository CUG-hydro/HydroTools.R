#' ET complementary
#'
#' @inheritParams ET0_Penman48
#' @param U2 2m wind speed (m/s)
#'
#' @references
#' 1. Xiao, M., Yu, Z., Kong, D., Gu, X., Mammarella, I., Montagnani, L., …
#'    Gioli, B. (2020). Stomatal response to decreased relative humidity
#'    constrains the acceleration of terrestrial evapotranspiration.
#'    Environmental Research Letters, 15(9). doi:10.1088/1748-9326/ab9967
#'
#' @examples
#' ET_CR_Xiao2020(250, 25, D = 1, U2 = 2)
#' @export
ET_CR_Xiao2020 <- function(Rn, Tair, D, U2, Pa = atm) {
    ea = cal_es(Tair) - D
    Tw = cal_Tw(ea, Tair, Pa)

    # b: a coefficient that depicts the proportion of sensible heat flux that increases ETp
    gamma = cal_gamma(Tair, Pa)
    slope = (cal_es(Tair) - cal_es(Tw)) / (Tair - Tw)
    b = slope/gamma # Xiao2020, Eq.3
    beta_wet = 1/b  # Xiao2020, Eq.9

    ET_p = ET0_Penman48(Rn, Tair, Pa, D, wind = U2, z.wind = 2)$ET0 # Shuttleworth 1993
    ET_w = ET0_eq(Rn, Tw, Pa)$Eeq # in the wet surface
    ET_a = ( (1 + b) * ET_w - ET_p ) / b
    data.table(Tair, Tw, ET_p, ET_w, ET_a)
}

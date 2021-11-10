#' Run Xinanjiang (XAJ) model (three sources, lumped).
#'
#' @description An R implementation of three-source Xinanjiang model
#'              by Renjun Zhao, used for daily streamflow simulation.
#' @param PREC Time series of precipitation (daily)
#' @param EVAP Time series of pan evaporation or potential evapotranspiration 
#' (daily), length must equal to `PREC`
#' @param params Parameters (see below)
#' @param area Basin area.
#' 
#' @param full.UH Use the unit hydrograph defined by user, rather than the
#'                instantaneous unit hydrograph (IUH) of Nash, for routing
#'                of surface runoff. Default FALSE.
#' 
#' @return 
#' This function returns a data frame of some common variables of the XAJ 
#' model at each time step, such as evaporation, soil moisture, surface and 
#' underground runoff.
#'
#' The variables of the table including:
#'
#' E,   Total evaporation (mm)
#' EU,  Evaporation (mm) of upper soil layer
#' EL,  Evaporation (mm) of lower soil layer
#' ED,  Evaporation (mm) of deep soil layer
#' W,   Total soil moisture (mm)
#' WU,  Soil moisture (mm) of upper soil layer
#' WL,  Soil moisture (mm) of lower soil layer
#' WD,  Soil moisture (mm) of deep soil layer
#' R,   Total runoff (mm) of each time step
#' RS,  Surface runoff (mm) of each time step
#' RI,  Interflow (mm) of each time step
#' RG,  Underground runoff (mm) of each time step
#' Q,   Total runoff (m^3/s) at the outlet of the basin
#' QS,  Surface runoff (m^3/s) at the outlet of the basin
#' QI,  Interflow runoff (m^3/s) at the outlet of the basin
#' QG,  Underground runoff (m^3/s) at the outlet of the basin
#'
#' @references 
#' Zhao and Liu, 1995. The Xinanjiang model, Computer Models of Watershed Hydrology, 
#' Water Resources Publication, Highlands Ranch, CO (1995), pp. 215-232
#' @export
XAJ <- function(PREC, EVAP, params, area = dt*3.6, dt = 24, full.UH = FALSE) {
  if(full.UH){
    UH <- params[14:length(params)]
  } else {
    # Create instantaneous unit hydrograph (IUH)
    UH <- XAJ:::IUH(params[14], params[15], 24)
  }

  # Run XAJ model.
  out <- data.frame(XAJ:::XAJrun(PREC, EVAP, params[1:13], UH, area, dt))
  names(out) <- c("E", "EU", "EL", "ED", "W", "WU", "WL", "WD",
                  "R", "RS", "RI", "RG", "Q", "QS", "QI", "QG")
  out
}


XAJ.param.names <- c("KC", "IM", "WUM", "WLM", "WDM", "C", "B", "SM", "EX", "KI", "KG", "CI", "CG", "N", "NK")

#' The lumped XAJ model has 13 parameters, including:
#'
#' The parameter `params` must be a numeric vector looks like:
#' `c(KC, IM, WUM, WLM, WDM, C, B, SM, EX, KI, KG, CI, CG, N, NK)`
#'
#' when use the instantaneous unit hydrograph of Nash, or looks like:
#' `c(KC, IM, WUM, WLM, WDM, C, B, SM, EX, KI, KG, CI, CG, UH_1, UH_2, ..., UH_n)`
#' UH_1, UH_2, ..., UH_n means the series of the unit hydrograph.
#'
#' @description
#' 1.  KC,   Ratio of potential evap to pan evap
#' 2.  IM,   Fraction of impermeable area
#' 3.  WUM,  Soil moisture capacity of upper layer
#' 4.  WLM,  Soil moisture capacity of lower layer
#' 5.  WDM,  Soil moisture capacity of deep layer
#' 6.  C,    Coefficient of deep evap
#' 7.  B,    Exponent of the soil moisture storage capacity curve
#' 8.  SM,   Areal mean free water capacity of the surface soil layer
#' 9.  EX,   Exponent of the free water capacity curve
#' 10. KI,   outflow coefficients of the free water storage to interflow
#' 11. KG,   outflow coefficients of the free water storage to groundwater
#' 12. CI,   recession constant of the lower interflow storage
#' 13. CG,   recession constant of groundwater storage.
#'      If use the instantaneous unit hydrograph (IUH) of Nash for routing
#'      of surface runoff, should provided two other parameters:
#'   - N,    number of reservoirs in the instantaneous unit hydrograph
#'   - NK,   common storage coefficient in the instantaneous unit hydrograph
#'      Else should provided the whole unit hydrograph defined by user,
#'      and set `full.UH` become `TRUE`.
XAJ.param.range <- data.frame(
    lower = c(0.20, 0.00, 5.0, 10.0, 10.0, 0.05, 0.1, 10.0, 0.50, 0.01, 0.01, 0.50, 0.95, 0.1, 1.0),
    upper = c(2.0, 0.2, 20., 90.0, 60.0, 0.20, 0.6, 100.0, 2.00, 0.70, 0.70, 0.90, 0.998, 5.0, 6.0),
    row.names = XAJ.param.names
)
# XAJ.param.range <- data.frame(
#     lower = c(0.20, 0.00, 5.0, 10.0, 10.0, 0.05, 0.1, 10.0, 0.50, 0.01, 0.01, 0.50, 0.95, 0.1, 1.0),
#     upper = c(1.50, 0.05, 20., 90.0, 60.0, 0.20, 0.6, 60.0, 2.00, 0.70, 0.70, 0.90, 0.998, 5.0, 6.0),
#     row.names = XAJ.param.names
# )
XAJ.param.range$par <- with(XAJ.param.range, lower + upper) / 2

XAJ.opts <- list(
    XAJ.param.range = XAJ.param.range,
    XAJ.param.names = XAJ.param.names
)

library(hydroTools)
library(sf)
library(data.table)
library(magrittr)

indir = "F:/CUG-Hydro/ChinaRunoff/basins/China/"
shp <- read_sf(glue("{indir}/ChinaBasins_(塔河+雅江+长江+黄河+珠江)-sp375_v0.1.0.shp"))
lst <- readRDS(glue("{indir}/ChinaBasins_气象驱动-sp375.RDS"))                  # meteorological forcing
load(glue("{indir}/ChinaBasins_(塔河+雅江+长江+黄河+珠江)-sp375_v0.1.0.rda"))   # runoff data

get_input <- function(site) {
  shp_i <- shp %>% subset(name == site)
  lat <- st_coordinates(shp_i) %>% {mean(.[, 2]) }
  area <- shp_i$area_w_km2 * 1e4 # km^2

  forcing <- lst[[site]][year(date) >= 2012] %>%
    set_names(c("date", "Rl", "prcp", "Pa", "q", "Rsi", "Tavg", "Tmax", "Tmin", "U10"))
  forcing = forcing %>%
    mutate(
      prcp = round(prcp * 24, 2), # mm hr-1, mm/d
      Pa = Pa / 1e3, # Pa, kPa
      Tavg = K2T(Tavg),
      Tmax = K2T(Tmax),
      Tmin = K2T(Tmin),
      # RH = q2RH(q, Tavg, Pa), # kg/kg, kPa
      ea = q2ea(q, Pa),
      # D2 = cal_es(Tavg) - ea,
      D = (cal_es(Tmin) + cal_es(Tmax))/2 - ea,
      Rn = cal_Rn(lat, date, Tmin, Tmax, ea = ea, Rsi = Rsi)$Rn,
      ET0 = ET0_FAO98(Rn, Tavg, Pa, D, U10, z.wind = 10)$ET0
    )
  data = df_river[name == site][year(date) >= 2002, -(1:2)] %>%
    dplyr::rename(Qobs = Q) %>%
    merge(forcing[, .(date, prcp, ET0)])
  list(data = data, area = area, lat = lat, site = site)
}
# srad: `W m-2`
# lrad: `W m-2`
# MJ m-2 d-1

# site <- "安康"
sites = c("府谷", "吴堡", "龙门", "河津", "潼关", "三门峡", "花园口")

# tmp = foreach(site = sites, i = icount()) %do% {
site = "河津"
    input = get_input(site)
    ## XAJ, daily scale
    df = input$data
    res_daily <- XAJ_calib(df$Qobs, df$prcp, df$ET0, date = df$date,
                          input$area, dt = 24, maxn = 1000)
    plot_calib(res_daily$data, main = site)

    outfile = glue("XAJ_daily_{site}.pdf")
    write_fig({
      par(family = "rTimes")
      plot_calib(res_daily$data, main = site)
    }, outfile)
# }

library(hydroTools)
library(sf)
library(data.table)
library(magrittr)

indir = "F:/CUG-Hydro/ChinaRunoff/basins/China/"
shp <- read_sf(glue("{indir}/ChinaBasins_(塔河+雅江+长江+黄河+珠江)-sp375_v0.1.0.shp"))
lst <- readRDS(glue("{indir}/ChinaBasins_气象驱动-sp375.RDS"))                  # meteorological forcing
load(glue("{indir}/ChinaBasins_(塔河+雅江+长江+黄河+珠江)-sp375_v0.1.0.rda"))   # runoff data

# site <- "安康"
site <- "河津"
shp_i <- shp %>% subset(name == site)
lat <- st_coordinates(shp_i) %>% {mean(.[, 2]) }
area <- shp_i$area_w_km2 * 1e4 # km^2


# srad: `W m-2`
# lrad: `W m-2`
# MJ m-2 d-1
{
  forcing <- lst[[site]][year(date) >= 2002] %>%
        set_names(c("date", "Rl", "prcp", "Pa", "q", "Rs", "Tavg", "Tmax", "Tmin", "U10"))
  dat = forcing %>%
    mutate(
      prcp = round(prcp * 24, 2), # mm hr-1, mm/d
      Pa = Pa / 1e3, # Pa, kPa
      Tavg = K2T(Tavg),
      Tmax = K2T(Tmax),
      Tmin = K2T(Tmin),
      # RH = q2RH(q, Pa, Tavg), # kg/kg, kPa
      ea = q2ea(q, Pa),
      # D2 = cal_es(Tavg) - ea,
      D = (cal_es(Tmin) + cal_es(Tmax))/2 - ea,
      Rn = cal_Rn(lat, date, Rs, Tmin, Tmax, ea = ea)$Rn,
      ET0 = ET0_FAO98(Rn, Tavg, Pa, D, U10, z.wind = 10)$ET0
    )
  df <- df_river[name == site][year(date) >= 2002, -(1:2)] %>%
    merge(dat[, .(date, prcp, ET0)])
}

## version 1, daily scale
r_day <- df %$% XAJ_calib(Q, prcp, ET0, area, dt = 24)
plot_calib(r_day, "XAJ_daily.pdf")


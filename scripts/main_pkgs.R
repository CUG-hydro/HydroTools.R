library(dplyr)
library(hydroTools)
library(lubridate)

#! /usr/bin/Rscript --no-init-file
# Dongdong Kong ----------------------------------------------------------------
# Copyright (c) 2021 Dongdong Kong. All rights reserved.
# source('scripts/main_pkgs.R')

# - `Rn`      : 地表净辐射，(MJ m-2 d-1)
# - `EVP`     : 蒸发皿蒸发，(0.1mm)
# - `TG`      : 地表温度，(0.1℃)
# - `Tair`    : 地表大气温度，(0.1℃)
# - `Prcp`    : 降水，(0.1mm)
# - `Pa`      : 大气压，(10 Pa)
# - `RH`      : 相对湿度，(%)
# - `SSD`     : 光照时数，(0.1 hour)
# - `WIN_avg` : 地表平均风速，(0.1 m s-1)

tidy_input <- function(d_raw, lat = 30) {
  d = d_raw[date <= "2019-12-31"] %>%
    # .[year(date) == 2019] %>%
    dplyr::select(site:date, EVP_sm:TG_min, prcp = `Prcp_20-20`,
           SSD = SSD, Tair_avg:Tair_min, U2 = WIN_Avg,
           Pa = Pa_avg, RH = RH_avg) %>%
    mutate(across(EVP_sm:U2, ~ .x %>% divide_by(10)),
           Pa = Pa/1e2)
    d %>% mutate(
        ssd0 = cal_ssd(lat, yday(date)),
        Rn = cal_Rn(lat, yday(date), RH = RH, ssd = SSD, Tmin = Tair_min, Tmax = Tair_max)$Rn %>% MJ_2W())
}

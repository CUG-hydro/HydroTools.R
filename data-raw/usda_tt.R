## code to prepare `usda_sf` dataset goes here

usda_tt = soiltexture::TT.env$TT.par$USDA.TT #%>% str()
usda_sf = usda_list2sf(usda_tt)

usethis::use_data(usda_sf, overwrite = TRUE)

#' delineate hydrological basin according to AecGIS tools
#'
#' @description
#' - add "c:/Python27/ArcGIS10.2/python2.exe" into system path
#' @export
fast_basin <- function(dem = "dem.tif") {
  tool <- system.file("python/fastbasin.bat", package = "hydroTools")
  cmd <- glue::glue("{tool} {dem}")
  shell(cmd)
}

# fill dem in R
# dem2 = topmodel::sinkfill(as.matrix(r_dem), 10^3, degree = 0.00001)
# values(r_dem) = as.numeric(dem2)
# raster::plot(r_dem)

# direc_grid = set_direc(direc_grid, 19, 11, 1)
# direc_grid = set_direc(direc_grid, 14, 3, 3)
# direc_grid = set_direc(direc_grid, 17, 8, 9)
# plot_flow_direc(direc_grid, river_shp=river_shp)

#write_arc_grid(direc_grid,'~/IMERG/direc.txt')

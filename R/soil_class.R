#' USDA soil class classification
#'
#' @docType data
#'
#' @references
#' 1. https://en.wikipedia.org/wiki/Soil_texture
"usda_sf"

# ' @importFrom sf st_polygon st_intersects st_sfc st_sf
#' @importFrom purrr map_int
usda_list2sf <- function(usda_tt) {
    data <- usda_tt$tt.points %>%
        as.matrix() %>%
        multiply_by(100) %>%
        set_colnames(c("clay", "silt", "sand"))
    polys <- map(usda_tt$tt.polygons, function(l) {
        ind <- l$points
        x <- data[ind, c(1, 3)] %>% rbind(.[1, ])
        sf::st_polygon(list(x))
    })
    g <- sf::st_sfc(polys)
    sf::st_sf(soilclass = names(polys), g) # st_usda
}

#' soil texture class, based on USDA triangle
#' 
#' @param clay Vector, percent of clay
#' @param sand Vector, percent of sand
#' 
#' @references 
#' 1. Julien Moeys (2016). soiltexture: Functions for Soil Texture Plot,
#'    Classification and Transformation. R package version 1.5.1.
#'    https://CRAN.R-project.org/package=soiltexture
#' 
#' @return Character vector of soil texture class, one of 
#' `c("Cl", "SiCl", "SaCl", "ClLo", "SiClLo", "SaClLo", "Lo", "SiLo", "SaLo", 
#' "Si", "LoSa", "Sa")`.
#' 
#' @examples
#' clay = c(05, 60, 15, 05, 25, 05, 25, 45, 65, 75, 13, 47)
#' sand = c(90, 32, 70, 70, 20, 10, 10, 10, 20, 10, 70, 10)
#' soil_class(clay, sand)
#' # c("Sa", "Cl", "SaLo", "SaLo", "SiLo", "Si", "SiLo", "SiCl", "Cl", "Cl", "SaLo", "SiCl")
#' @export
soil_class <- function(clay, sand) {
    points <- data.frame(clay, sand) %>% sf::st_as_sf(coords = c("clay", "sand"))
    # r2 = st_intersection(points, usda_sf)
    ans <- sf::st_intersects(points, usda_sf)
    ind <- map_int(ans, dplyr::first)
    usda_sf$soilclass[ind]
}

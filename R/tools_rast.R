rast_info <- function (file) {
    # subfix = stringr::str_extract(basename(file), "(?<=\\.).{2,4}$")
    # if (subfix == "nc") {
    #     ncdim_get(file)
    # }
    # else {
        suppressWarnings(x <- rgdal::GDALinfo(file))
        lon = x %>% {
            seq(.["ll.x"] + .["res.x"]/2, by = .["res.x"],
                length.out = .["columns"])
        }
        lat = x %>% {
            seq(.["ll.y"] + .["res.y"]/2, by = .["res.y"],
                length.out = .["rows"])
        }
        listk(lon, lat)
    # }
}

#' @importFrom raster raster area as.array
rast_coord <- function (r, .area = TRUE){
    if (is.character(r))
        r %<>% raster()
    r %<>% raster()
    pos <- raster::coordinates(r) %>% 
        cbind(values(area(r)),  values(r)) %>% 
        set_colnames(c("lon", "lat", "area", "value")) %>% data.table()
    r_temp <- r
    values(r_temp) <- 1:prod(dim(r)[1:2])
    I_grid <- rast_array(r_temp) %>% as.numeric()
    nrow <- ncol(r)
    pos[I_grid, ] %>% mutate(I = 1:nrow(.), col = ceiling(I/nrow),   
        row = I - (col - 1) * nrow) %>% reorder_name("I",
        tailvars = c("row", "col"))
}

rast_array <- function(r) {
    if (is.character(r)) {
        r %<>% raster()
    }
    aperm(as.array(r), c(2, 1, 3)) %>% flipud()
}

flipud <- function(x, ...) {
    I <- ncol(x):1
    ndim <- length(dim(x))
    if (ndim == 2) {
        x[, I]
    } else if (ndim == 3) {
        x[, I, ]
    }
}

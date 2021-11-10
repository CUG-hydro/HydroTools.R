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



#' Read the rows and columns of grids of a basin from a rout data file.
#'
#' @description Read the rows and columns of gridcells of a basin from a rout data file
#'              created by VIC Hime.
#'
#' @param data_path File path of the rou data file.
#' @return A data frame contains the rows and columns of the gridcells in the basin.
#'
#' @export
read_basin <- function(data_path) {
  basin <- fromJSON(readLines(data_path))$basin
  basin <- data.frame(t(data.frame(basin)))
  rownames(basin) <- 1:nrow(basin)
  colnames(basin) <- c('col', 'row')
  return(basin)
}

#' Find out the gridcells in the basin controled by the station.
#'
#' @description Offer the grid of flow direction and the site (row and column)
#'              of the hydrological station and returns a table of the coordinates
#'              of the gridcells in the basin controled by the station.
#'
#' @param dir Direction grid created by create_flow_direction() or a matrix.
#' @param stn_col Column of the hydrological station
#' @param stn_row Row of the hydrological station
#' @param arc_code If use the ArcInfo direction code. Default True.
#' @param xcor X coordinate of the corner of the grid. Should offered if dir is a matrix.
#' @param ycor Y coordinate of the corner of the grid. Should offered if dir is a matrix.
#' @param csize Size of the gridcells. Should offered if dir is a matrix.
#' @param get_coords If return the coordinates of the gridcells of the basin, or return the rows and columns.
#' @return A data frame contains the coordinates of the gridcells in the basin.
#'
#' @export
detect_basin <- function(dir, stn_col, stn_row, arc_code=TRUE, xcor=NULL, ycor=NULL, csize=NULL, get_coords=FALSE) {
  if(is.list(dir)) {
    xcor <- dir$xllcorner
    ycor <- dir$yllcorner
    csize <- dir$cellsize
    dir <- dir$grid
  }

  grid_satu <- dir
  grid_satu[!is.na(grid_satu)] <- 0
  grid_satu[stn_row, stn_col] <- 1

  dx <- c()
  dy <- c()
  if (arc_code) {
    dx[c(1,2,4,8,16,32,64,128)] <- c(1,1,0,-1,-1,-1,0,1)
    dy[c(1,2,4,8,16,32,64,128)] <- c(0,-1,-1,-1,0,1,1,1)
  } else {
    dx[c(3,4,5,6,7,8,1,2)] <- c(1,1,0,-1,-1,-1,0,1)
    dy[c(3,4,5,6,7,8,1,2)] <- c(0,-1,-1,-1,0,1,1,1)
  }

  bx <- c()
  by <- c()

  nr <- nrow(dir)
  nc <- ncol(dir)
  for(rw in 1:nr) {
    for(cl in 1:nc) {
      if(is.na(dir[rw, cl])) next
      if(grid_satu[rw, cl] < 0) next

      x <- cl
      y <- rw

      # Anti roop
      len_det <- 10
      prex <- rep(0, len_det)
      prey <- rep(0, len_det)

      while(TRUE) {
        in_loop <- FALSE
        for(p in 1:len_det){
          if(prex[p] == x & prey[p] == y) {
            in_loop <- TRUE
          }
        }
        if(in_loop) {print(paste("Warn: loop detect at", x, ",", y))
          grid_satu[y, x] <- -2
          break
        }

        prex[1:(len_det-1)] <- prex[2:len_det]
        prey[1:(len_det-1)] <- prey[2:len_det]

        prex[5] <- x
        prey[5] <- y

        if(x < 1 | y < 1 | x > nc | y > nr) {
          grid_satu[rw, cl] <- -1
          break
        }
        if(is.na(dir[y, x]) | grid_satu[y, x] < 0) {
          grid_satu[rw, cl] <- -1
          break
        }
        if(grid_satu[y, x] == 1) {
          grid_satu[rw, cl] <- 1
          bx <- append(bx, cl)
          by <- append(by, rw)
          break
        }
        d <- dir[y, x]
        nx <- x + dx[d]
        ny <- y + dy[d]
        x=nx
        y=ny
      }
    }
  }
  if(get_coords) {
    bx <- bx*csize+xcor-csize/2
    by <- by*csize+ycor-csize/2
  }
  basin <- data.frame(x=bx, y=by)
  return(basin)
}

# direc_grid = set_direc(direc_grid, 19, 11, 1)
# direc_grid = set_direc(direc_grid, 14, 3, 3)
# direc_grid = set_direc(direc_grid, 17, 8, 9)
# plot_flow_direc(direc_grid, river_shp=river_shp)

#write_arc_grid(direc_grid,'~/IMERG/direc.txt')

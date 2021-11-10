.prj84 = sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

#' Convert the grid data (save as matrix) to points.
#'
#' @description Conver the grid data to points, showed as a table, with the columns of
#'              coordinates for each point and a column for their value.
#'
#' @param grid A matrix of the gridded data.
#' @param x X coordinates for each column of the grid.
#' @param y Y coordinates for each column of the grid.
#' @param csize Size of each gridcell.
#' @param xcor X coordinate of the southwest corner of the grid.
#' @param xcor Y coordinate of the southwest corner of the grid.
#'
#' @return A table of the points, including x and y coordinate and the values, converted
#'         from the grid.
#' @export
grid2points <- function(grid, x=NULL, y=NULL, csize=NULL, xcor=NULL, ycor=NULL, NA_value=NULL, y_rev=FALSE) {
  if((is.null(x) | is.null(y)) & (is.null(csize) | is.null(xcor) | is.null(ycor))) {
    stop("Must provide x, y or csize, xcor, ycor")
  }

  if(is.null(x) | is.null(y)) {
    x <- (1:ncol(grid)) * csize - csize/2 + xcor
    y <- (1:nrow(grid)) * csize - csize/2 + ycor
  }
  if (y_rev) y <- rev(y)

  xs <- c()
  ys <- c()
  vs <- c()
  if(!is.null(NA_value)) grid[grid == NA_value] = NA
  for(row in 1:nrow(grid)) {
    for(col in 1:ncol(grid)) {
      if(is.na(grid[row, col])) next

      xs <- append(xs, x[col])
      ys <- append(ys, y[row])
      vs <- append(vs, grid[row, col])
      }
  }
  points <- data.frame(xs, ys, vs)
  names(points) <- c('x', 'y', 'value')
  return(points)
}


#' Convert the point data (a table including columns of the coordinates and value) to
#' grid.
#'
#' @description Conver the point data, usually a table with columns including coordinates
#'              and values to grids.
#'
#' @param grid A matrix of the gridded data.
#' @param x Vector of x coordinates or hich column store the x cordinate.
#' @param y Which column store the y cordinate.
#' @param val Which column store the point value.
#'
#' @return A matrix of the gridcell value.
#' @export
points2grid <- function(points, x=NULL, y=NULL, val=NULL) {

  if (is.null(points) & (is.null(x) | is.null(y)))
    stop("Must provide points or x and y.")
  if (is.null(x)) {
    xs <- points[, 1]
  }else {
    xs <- x
  }
  if (is.null(y)) {
    ys <- points[, 2]
  } else {
    ys <- y
  }
  if(is.null(val)) {
    if(nrow(points >= 3)) {
      v = points[ ,3]
    }else{
      v=rep(0,length(xs))
    }
  } else {
    v = val
  }
  ux <- sort(unique(xs))
  uy <- sort(unique(ys))
  lux <- length(ux)
  luy <- length(uy)
  itvx <- unique(ux[2:lux] - ux[1:(lux - 1)])[1]
  itvy <- unique(uy[2:luy] - uy[1:(luy - 1)])[1]

  cellsize <- mean(c(itvx[1], itvy[1]))
  xcor <- min(ux) - cellsize/2
  ycor <- min(uy) - cellsize/2
  cols <- round((xs - xcor + cellsize/2)/cellsize)
  rows <- round((ys - ycor + cellsize/2)/cellsize)
  nrow <- max(rows)
  ncol <- max(cols)
  grid <- matrix(nrow = ncol, ncol = nrow)
  for (p in 1:nrow(points)) {
    grid[cols[p], rows[p]] <- v[p]
  }
  grid
}

# Assistant function.
get_border <- function(nc, x, y, xsize, ysize) {
  flons <- nc$dim$x$vals
  flats <- nc$dim$y$vals
  if(is.null(flons) | is.null(flats)) {
    flons <- nc$dim$lon$vals
    flats <- nc$dim$lat$vals
  }
  l <- x - xsize/2
  r <- x + xsize/2
  t <- y + ysize/2
  b <- y - ysize/2
  cols <- which(flons > l & flons < r)
  rows <- which(flats > b & flats < t)
  ncl <- length(cols)
  nrw <- length(rows)
  lf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[1]), c(1, nrw))
  rf <- ncvar_get(nc, nc$var[[1]], c(cols[ncl], rows[1]), c(1, nrw))
  tf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[1]), c(ncl, 1))
  bf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[nrw]), c(ncl, 1))
  lf[is.na(lf)] <- 0
  rf[is.na(rf)] <- 0
  tf[is.na(tf)] <- 0
  bf[is.na(bf)] <- 0

  lf <- rev(lf)
  bf <- rev(bf)
  maxps <- c(which.max(tf), which.max(rf), which.max(bf), which.max(lf))
  maxs <- c(max(tf), max(rf), max(bf), max(lf))
  maxside <- which.max(maxs)

  if(maxside == 1 | maxside == 3) {
    maxp <- maxps[maxside]/ncl
  } else {
    maxp <- maxps[maxside]/nrw
  }
  return(c(maxside, maxp))
}

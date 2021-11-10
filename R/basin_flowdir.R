
#' Create the flow direction of each gridcell for VIC routing model.
#'
#' @description Create the flow direction data of each VIC gridcell for the routing model
#'              of VIC by providing a high resolution flow accumulation data (about 90m
#'              usually caculated from SRTM DEM by ArcGIS) and the domain file of VIC model.
#'              WARNING: This is a simple heuristic method and more fit for small gridcell.
#'              It's usually needs further manual calibration by set_direc and plot_flow_direc.
#'
#' @param flow_file An netCDF4 file of the flow accumulation data from high resolution DEM
#'                  (about 90m). Must containing all the area.
#' @param domain_file Domain file (netCDF4 format) of VIC model.
#' @param arc_code If the direction code is ArcInfo type. If FALSE, it would set to 1 to 8.
#'
#' @param diag_tsh Diag threshould. Determine if flow to the diag gridcell.
#'
#' @return A list including the meta parameters (num of rows and columns, size of the cells,
#'         x and y corner) and the matrix of the flow direction data.
#' @import ncdf4
#'
#' @export
create_flow_direction <- function(flow_file, domain_file, arc_code=TRUE, diag_tsh=0.382) {
  if(arc_code) {
    dir_code <- c(64, 128, 1, 2, 4, 8, 16, 32)
  } else {
    dir_code <- 1:8
  }

  grids_nc <- nc_open(domain_file)
  grids_var <- ncvar_get(grids_nc, grids_nc$var[[1]])
  glons <- grids_nc$dim$lon$vals
  glats <- grids_nc$dim$lat$vals
  nc_close(grids_nc)

  nrows <- length(glons)
  ncols <- length(glats)
  csizex <- abs(mean(glons[2:nrows]-glons[1:(nrows-1)]))
  csizey <- abs(mean(glats[2:ncols]-glats[1:(ncols-1)]))

  flow_nc <- nc_open(flow_file)

  direc <- grids_var
  direc[!is.na(direc)] <- 0
  for(col in 1:ncols) {
    for(row in 1:nrows) {
      if(is.na(direc[row, col]) | direc[row, col] > 0) next

      glon <- glons[row]
      glat <- glats[col]

      mf <- get_border(flow_nc, glon, glat, csizex, csizey)
      mfd <- mf[1]
      mfs <- mf[2]

      nextrow <- row
      nextcol <- col
      if(mfd == 1) {
        nextcol <- col - 1
      } else if(mfd == 2) {
        nextrow <- row + 1
      } else if(mfd == 3) {
        nextcol <- col + 1
      } else {
        nextrow <- row - 1
      }

      if(mfd == 1) {
        old <- dir_code[8]
        osd <- dir_code[1]
        ord <- dir_code[2]
      } else if(mfd == 2) {
        old <- dir_code[2]
        osd <- dir_code[3]
        ord <- dir_code[4]
      } else if(mfd == 3) {
        old <- dir_code[4]
        osd <- dir_code[5]
        ord <- dir_code[6]
      } else if(mfd == 4) {
        old <- dir_code[6]
        osd <- dir_code[7]
        ord <- dir_code[8]
      }

      if(mfd == 1 & nextcol < 1 | mfd == 2 & nextrow > nrows |
         mfd == 3 & nextcol > ncols | mfd == 4 & nextrow < 1) {
        direc[row, col] <- osd
        next
      }

      nglon <- glons[nextrow]
      nglat <- glats[nextcol]

      nmf <- get_border(flow_nc, nglon, nglat, csizex, csizey)
      nmfd <- nmf[1]
      nmfs <- nmf[2]

      nmfd <- (nmfd + 4 - mfd) %% 4 + 1

      if(nmfd == 3) {
        direc[row, col] <- osd
        direc[nextrow, nextcol] <- osd
      } else if(nmfd == 4 & mfs < diag_tsh & nmfs < diag_tsh) {
        direc[row, col] <- old
      } else if(nmfd == 2 & mfs > 1-diag_tsh & nmfs > 1-diag_tsh) {
        direc[row, col] <- ord
      } else {
        direc[row, col] <- osd
      }
    }
  }
  direc <- t(direc)
  nc_close(flow_nc)
  arcgrid <- list('ncols'=nrows,
                  'nrows'=ncols,  # The row and col is inversed of nc file.
                  'xllcorner'=min(glons)-csizex/2,
                  'yllcorner'=min(glats)-csizey/2,
                  'cellsize'=(csizex+csizey)/2,
                  'grid'=direc)
  return(arcgrid)
}


# river_shp='~/IMERG/predata/dense_river.shp'

# flow_file <- '~/IMERG/predata/flow_bj.nc'
# domain_file <- '~/IMERG/predata/grids.nc'

# direc_grid=create_flow_direction(flow_file, domain_file, diag_tsh=0.3)

#' Plot the flow direction of the flow direction data.
#'
#' @description Plot the flow direction of a flow direction data (created by
#'              create_flow_direction()) for visullization and calibration. Can
#'              provide a shp file of the river network to plot helping the
#'              mannual calibration.
#'
#' @param direc Flow direction data created by create_flow_direction().
#' @param arc_code If use the ArcInfo direction code.
#' @param river_shp Path of the river shp file.
#' @import shape maptools
#'
#' @export
plot_flow_direc <- function(direc, xlim = NULL, ylim = NULL, arc_code = TRUE, river_shp = NULL, row_rev = FALSE) {
    if (arc_code) {
        dir_code <- c(64, 128, 1, 2, 4, 8, 16, 32)
    } else {
        dir_code <- 1:8
    }

    if (!is.data.frame(direc)) {
        csize <- direc$cellsize
        xcor <- direc$xllcorner
        ycor <- direc$yllcorner
        direc <- grid2points(direc$grid, xcor = xcor, ycor = ycor, csize = csize, y_rev = row_rev)
    }

    xticks <- sort(unique(direc[, 1]))
    yticks <- sort(unique(direc[, 2]))
    cx <- xticks[2] - xticks[1]
    cy <- yticks[2] - yticks[1]
    xlabs <- round((xticks - min(xticks)) / (cx) + 1)
    ylabs <- round((yticks - min(yticks)) / (cy) + 1)

    if (is.null(xlim)) {
        xlim <- c(min(xticks), max(xticks))
    } else {
        xlim <- xticks[xlim]
    }
    if (is.null(ylim)) {
        ylim <- c(min(yticks), max(yticks))
    } else {
        ylim <- yticks[ylim]
    }

    plot(direc[, 1:2],
        cex = 0, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim,
        xlab = NA, ylab = NA, asp = 1
    )
    axis(1, xticks, xlabs)
    axis(3, xticks, xlabs)
    axis(2, yticks, ylabs)
    axis(4, yticks, ylabs)
    abline(h = yticks, col = "grey")
    abline(v = xticks, col = "grey")

    if (!is.null(river_shp)) {
        rs <- readShapeLines(river_shp)
        lines(rs, col = "blue")
    }
    points(direc[, 1:2], pch = 20, cex = 0.8)

    xcz <- sort(unique(direc[, 1]))
    xcz <- mean(xcz[2:length(xcz)] - xcz[1:(length(xcz)) - 1])
    ycz <- sort(unique(direc[, 2]))
    ycz <- mean(ycz[2:length(ycz)] - ycz[1:(length(ycz)) - 1])

    for (p in 1:nrow(direc)) {
        dx <- 0
        dy <- 0
        x <- direc[p, 1]
        y <- direc[p, 2]
        d <- direc[p, 3]
        if (d == dir_code[1]) {
            dy <- ycz
        } else if (d == dir_code[2]) {
            dy <- ycz
            dx <- xcz
        } else if (d == dir_code[3]) {
            dx <- xcz
        } else if (d == dir_code[4]) {
            dy <- -ycz
            dx <- xcz
        } else if (d == dir_code[5]) {
            dy <- -ycz
        } else if (d == dir_code[6]) {
            dy <- -ycz
            dx <- -xcz
        } else if (d == dir_code[7]) {
            dx <- -xcz
        } else if (d == dir_code[8]) {
            dx <- -xcz
            dy <- ycz
        }
        # arrows(c(x, x+dx), c(y, y+dy))
        Arrows(x, y, x + dx, y + dy, arr.length = 0.15, arr.adj = 0.15)
    }
}

#' Revise the flow direction of the flow direction data.
#'
#' @description Revise the flow direction of the flow direction data created by
#'              create_flow_direction().
#'
#' @param direc_grid Flow direction data created by create_flow_direction().
#' @param row Row of the gridcell to be revised.
#' @param col Column of the gridcell to be revised.
#' @param direc New flow direction of the gridcell. For more convinient calibration,
#'              the direction code is sat as the keypad. It means that 1 is southwest,
#'              2 is south, 3 is southeast, and so on. set keypad_dir to FALSE can
#'              revise the direction code directly.
#' @param arc_code If use the ArcInfo flow direction code.
#' @param keypad_dir If use keypad flow direction code.
#'
#' @return The flow direction data after revise.
#'
#' @export
set_direc <- function(direc_grid, row, col, direc, arc_code = TRUE, keypad_dir = TRUE) {
    if (!is.list(direc_grid) | is.null(direc_grid$nrow | is.null(direc_grid$grid))) {
          stop("direc_grid is not a ArcInfo grid.")
      }
    if (arc_code) {
        dir_map <- c(8, 4, 2, 16, NA, 1, 32, 64, 128)
    } else {
        dir_map <- c(6, 5, 4, 7, NA, 3, 8, 1, 2)
    }
    nrow <- direc_grid$nrows
    if (keypad_dir) di <- dir_map[direc] else di <- direc
    direc_grid$grid[row, col] <- di
    return(direc_grid)
}

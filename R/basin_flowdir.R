flowdir <- function(file_dir, arcinfo = TRUE) {
    r <- raster::raster(file_dir)
    grid <- rast_array(r)
    # grid <- rast_array(file_dir)

    to_x <- c(1, 1, 0, -1, -1, -1, 0, 1)
    to_y <- c(0, -1, -1, -1, 0, 1, 1, 1)

    if (arcinfo) {
        from <- c(1, 2, 4, 8, 16, 32, 64, 128)
    } else {
        from <- c(3:8, 1:2)
    }

    grid_dx <- mapvalues(grid, from, to_x)
    grid_dy <- mapvalues(grid, from, to_y)

    pos <- rast_coord(file_dir)
    info <- data.table(dir = as.numeric(grid), dx = as.numeric(grid_dx), dy = as.numeric(grid_dy)) %>%
        cbind(pos, .)
    info[, angle := asin(dy / sqrt(dx^2 + dy^2))]
    info
}

flowdir_find_next <- function(info){
    csize = info$lon %>% diff() %>% median()
    # ng <- length(ic)
    ng = nrow(info)
    nextg <- rep(0, ng)
    dist <- rep(0., ng)
    # find next grid
    for (i in 1:ng) {
        # idir <- paste(grid[ir[i], ic[i]]) # note idir is string
        # if (idir <= 0) {
        #     nextg[i] <- 0
        #     next
        # }
        if (is.na(info$dir[i]) || info$dir[i] <= 0) next()

        lat <- info$lat[i]
        lat2 <- lat + info$dy[i] * csize
        # https://en.wikipedia.org/wiki/Great-circle_distance
        dist[i] <- 6371229 * acos(sin(lat * pi / 180) * sin(lat2 * pi / 180) +
            cos(lat * pi / 180) * cos(lat2 * pi / 180) * cos(info$dx[i] * csize * pi / 180))

        nx <- info$row[i] + info$dx[i]
        ny <- info$col[i] + info$dy[i]
        inextg <- info[row == nx & col == ny, I]
        # inextg <- which(ir == nx & ic == ny)
        nextg[i] <- ifelse(length(inextg) == 0, 0, inextg[1])
    }
    ans = data.table(nextg, dist, angle_next = NA_real_) %>% cbind(info, .)
    ans[nextg > 0, angle_next := angle[nextg]]
    ans
}

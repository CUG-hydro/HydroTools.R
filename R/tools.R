.prj84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

# remove NA_real_ and Inf values in Y_sim, Y_obs and w
valid <- function(x) !is.na(x) & is.finite(x)

valindex <- function(obs, sim, ...) {
    if (length(obs) != length(sim)) {
        stop("Invalid argument: 'length(sim) != length(obs)' !! (", length(sim), "!=", length(obs), ") !!")
    } else {
        index <- which(valid(obs) & valid(sim))
        if (length(index) == 0) {
            warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
        }
        return(index)
    }
}

is_empty <- function(x) {
    is.null(x) || (is.data.frame(x) && nrow(x) == 0) || length(x) == 0
}

clamp <- function(x, lims = c(0, 1), fill.na = FALSE) {
    if (fill.na) {
        x[x < lims[1]] <- NA_real_
        x[x > lims[2]] <- NA_real_
    } else {
        x[x < lims[1]] <- lims[1]
        x[x > lims[2]] <- lims[2]
    }
    x
}

check_doy <- function(x) {
    if (is.Date(x)) x %<>% yday()
    x
}

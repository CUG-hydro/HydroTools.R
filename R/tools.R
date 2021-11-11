.prj84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

# remove NA_real_ and Inf values in Y_sim, Y_obs and w
valid <- function(x) !is.na(x) & is.finite(x)

valindex <- function(obs, sim, ...) {
    if (length(obs) != length(sim)) {
        stop("Invalid argument: 'length(sim) != length(obs)' !! (", length(sim), "!=", length(obs), ") !!")
    } else {
        index <- which(valid(obs) & valid(sim))
        if (length(index) == 0) 
            warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
        return(index)
    } # ELSE end
}

is_empty <- function(x) {
    is.null(x) || (is.data.frame(x) && nrow(x) == 0) || length(x) == 0
}

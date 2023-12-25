# .prj84 <- "epsg:4326"

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

listk <- function(...) {
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_noname <- if (is.null(vars)) {
    seq_along(cols)
  } else {
    which(vars == "")
  }
  if (length(Id_noname) > 0) {
    vars[Id_noname] <- sapply(cols[Id_noname], deparse)
  }
  x <- setNames(list(...), vars)
  return(x)
}

set_dim <- function(x, dim) {
  dim(x) <- dim
  x
}


#' @importFrom lubridate make_date year month day
make_monthdate <- function(x) {
  x %>% {lubridate::make_date(year(.), month(.), 1)}
}

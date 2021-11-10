runningId <- function(i, step = 1, N, prefix = "") {
    perc <- ifelse(missing(N), "", sprintf(", %.1f%% ", i/N*100))
    if (mod(i, step) == 0) cat(sprintf("[%s] running%s %d ...\n", prefix, perc, i))    
}

apply_row <- function(mat, by, FUN = rowMeans2, ...) {
    if (length(by) != ncol(mat)) {
        stop('Length of by is not equal to ncol of mat')
    }
    grps <- unique(by) %>% sort()

    foreach(grp = grps, .combine = cbind) %do% {
        I <- which(by == grp)
        FUN(mat[, I], na.rm = TRUE, ...)
    } %>% set_colnames(grps) %>% 
    set_rownames(rownames(mat))
}

# ' @importFrom base intersect setdiff union
#' @importFrom data.table data.table is.data.table as.data.table
reorder_name <- function(
    d,
    headvars = c("site", "date", "year", "doy", "d8", "d16"),
    tailvars = "")
{
    names <- names(d)
    headvars %<>% intersect(names)
    tailvars %<>% intersect(names)
    varnames <- c(headvars, setdiff(names, union(headvars, tailvars)), tailvars)

    if (is.data.table(d)) {
        d[, varnames, with = F]
    } else if (is.data.frame(d)) {
        d[, varnames]
    } else if (is.list(d)){
        d[varnames]
    } else{
        stop("Unknown data type!")
    }
}

# pacman::p_load(data.table, dplyr, lubridate)
#' @importFrom lubridate ddays
get_event_info <- function(date, Q, extend = ddays(5)) {
  date_beg <- min(date) - extend
  date_end <- max(date) + extend
  days <- date_end - date_beg + 1
  data.table(date_beg, date_end, days,
    Q_min = min(Q, na.rm = TRUE),
    Q_max = max(Q, na.rm = TRUE),
    Q_mean = mean(Q, na.rm = TRUE)
  )
}

#' @param gap_max distance smaller than gap_max is considered as a same group
#'
#' @rdname flood_divide
#' @export
detect_groups <- function(df, inds, gap_max = 5, extend = ddays(5)) {
  grps <- cumsum(c(0, diff(inds) > gap_max)) # 分组
  info <- cbind(group = grps, df[inds, ]) # 分组
  info[, get_event_info(date, Q, extend), group]
}

#' detect_flood_events
#'
#' @param Q_min minimum discharge to detect flood events
#' @param Q_peak peak discharge to detect flood events
#'
#' @param gap_max Default `5`. if `index gap > gap_max`, events will be regarded as two.
#' If Q is hourly, gap_max can be set `5*24`.
#' @param extend Default `ddays(5)`. Extend `nday` in the left and right of a event
#'
#' @rdname flood_divide
#' @export
detect_flood_events <- function(date, Q, Q_min = 2, Q_peak = 10, gap_max = 5, extend = ddays(5)) {
  df <- data.table(date, Q)
  inds <- df[, which(Q > Q_min)] #

  info_group <- detect_groups(df, inds, gap_max, extend) %>%
    .[Q_max > Q_peak, ]

  ## 数据压缩
  lgl <- rep(FALSE, nrow(df))
  for (i in 1:nrow(info_group)) {
    info <- info_group[i, ]
    lgl[date >= info$date_beg & date <= info$date_end] <- TRUE
  }

  # 由于已经扩展了5天，这里对洪水事件重新进行编号，不需要二次扩展
  inds <- which(lgl)
  info_group <- detect_groups(df, inds, gap_max, extend = 0)
  info_group
  # listk(group = info_group, index = inds) #
}

#' flood_divide
#' @param df A data.table with date and Q columns
#' @param ... parameters passed to [detect_flood_events()]
#' @export
flood_divide <- function(df, ...) {
  date <- df$date
  if ("time" %in% names(df)) date <- df$time
  df <- df[order(date)]
  detect_flood_events(date, df$Q, ...) # info_group
}

#' @rdname flood_divide
#' @export
merge_flood <- function(df, info_flood, format = "%Y.%m") {
  if ("time" %in% names(df)) time <- df$time
  r <- map(1:nrow(info_flood), function(i) {
    info <- info_flood[i, ]
    date_beg <- info$date_beg
    date_end <- info$date_end
    df[time >= date_beg & time <= date_end]
  }) %>% melt_list("group")
  r[, group_name := format(time[1], "%Y.%m.%d"), .(group)]
  relocate(r, group, group_name)
}

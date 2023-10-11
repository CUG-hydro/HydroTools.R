# set_lims <- function(date_begin = 2012, date_end = date_begin) {
#   if (is.numeric(date_begin)) date_begin %<>% make_datetime(., 1, 1)
#   if (is.numeric(date_end)) date_end %<>% make_datetime(., 12, 31, 23)
# 
#   lims <- c(date_begin, date_end)
#   list(scale_x_datetime(limits = lims))
# }

# year vertical grids
#' grid_date
#' 
#' @param ... others to [abline()]
#' 
#' @keywords internal
#' @importFrom lubridate make_date make_datetime
#' @export
grid_date <- function(x, col = "grey60", lty = 3, lwd = 0.4, ...) {
  year_min = year(min(x)) - 1
  year_max = year(max(x)) + 1

  if ("POSIXct" %in% class(x)) {
    date_begin = make_datetime(year_min)
    date_end = make_datetime(year_max, 12, 31, 23)
  } else if ("Date" %in% class(x)) {
    date_begin = make_date(year_min)
    date_end = make_date(year_max, 12, 31)
  }
  is_short = (year_max - year_min + 1 - 2) <= 4
  by = ifelse(is_short, "3 month", "year")
  t_grids <- seq(date_begin, date_end, by = by)
  abline(v = t_grids, col = col, lty = lty, lwd = lwd, ...)
  if (is_short) {
    abline(v = seq(date_begin, date_end, by = "year"), col = col, lty = 2, lwd = 0.8, ...)
  }
}

#' plot_runoff
#'
#' @param df_prcp A data.frame with the columns of `date`, `prcp`
#' @param df_q A data.frame with the columns of `date`, `Qsim`, `Qobs`. 
#' At least one of `Qsim` and `Qobs` should be provided.
#' @param xlim limit of x axis
#' @param ylim2 limit of second y axis (for precipitation)
#' @param ... other parameters passed to [plot()] for precipitation plot
#' 
#' @export
plot_runoff <- function(df_q, df_prcp = NULL, xlim = NULL, ylim2 = c(50, 0), 
  legend.position = c(1, 1), legend.justification = c(1, 1),
  ...) {
  if (is.null(xlim)) {
    # 确保prcp和Q的xlim相同
    lim1 = range(df_prcp$date, na.rm = TRUE)
    lim2 = range(df_q$date, na.rm = TRUE)
    # 以df_q的日期为主
    # xlim = c(min(lim1[1], lim2[1]), max(lim1[2], lim2[2]))
    xlim = lim2
  }

  str_gof = ""
  if (!is.null(df_q$Qobs) && !is.null(df_q$Qsim)) {
    gof = GOF(df_q$Qobs, df_q$Qsim) %>% as.list()
    str_gof = do.call(sprintf, c(
      "NSE=%.3f, KGE=%.3f, RMSE=%.1f, \nBias=%.1f, R2=%.3f, n=%d",
      gof[c("NSE", "KGE", "RMSE", "Bias", "R2", "n_sim")]
    ))
  }

  lab_Q = expression(bold(Q * "(" * m^3 * "/s)"))
  par(mar = c(3.5, 3.5, 2, 3.5), mgp = c(1.8, 0.6, 0), xpd = TRUE)
  if (!is.null(df_q$Qobs)) {
    plot(Qobs~date, df_q, type = "l", col = "black", xlab = "Date",
      ylab = lab_Q,
      # font.lab = 2,
      # xaxt = "n", 
      xlim = xlim)
    lines(df_q$date, df_q$Qsim, col = "red")
  } else {
    plot(df_q$date, df_q$Qsim, type = "l", col = "red", xlab = "Date",
      ylab = lab_Q,
      xaxt = "n", 
      xlim = xlim)
  }

  # mtext("Precipitation (mm/h)", side = 4, col = "blue")
  corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  legend(corners[1], corners[4], legend = c("Qobs", "Qsim"),
         yjust = 0, col = c("black", "red"), bty = "n",
         lty = 1, ncol = 3)

  # TODO: 调整gof的位置
  norm_trans <- function(xy) {
    usr = par("usr") # (x1, x2, y1, y2)
    xrange = diff(usr[1:2])
    yrange = diff(usr[3:4])
    c(xy[1] * xrange + usr[1], xy[2] * yrange + usr[3])
  }
  legend.position %<>% norm_trans()
  xjust = legend.justification[1]
  yjust = legend.justification[2]
  legend(legend.position[1], legend.position[2], str_gof, 
    xjust = xjust, yjust = yjust,
    text.col = "red", bty = "n")
  # ADD grid
  par(xpd = FALSE)
  grid_date(df_q$date)

  if (!is.null(df_prcp)) {
    par(new = TRUE, xpd = FALSE)
    plot(prcp ~ date, df_prcp,
      ylim = ylim2, xlim = xlim, type = "h", col = "blue",
      yaxs = "i",
      xaxt = "n", yaxt = "n",
      ylab = "", xlab = "", ...
    )
    axis(side = 4, at = pretty(range(df_prcp$prcp, na.rm = TRUE)),
         col = "blue", col.ticks = "blue", col.axis = "blue")
    corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
    par(xpd = TRUE) #Draw outside plot area
    text(x = corners[2], y = mean(corners[3:4]), "Precipitation (mm/h)",
         col = "blue", font = 2, srt = 270, adj = c(0.95, -3.5))
    # mtext("Precipitation (mm/h)", side = 4, col = "blue")
    par(xpd = FALSE)
  }
}

#' @export
plot_gof <- function(ysim, yobs, date = NULL) {
  par(mar = c(3, 2, 2, 1), mgp = c(3, 0.6, 0), mfrow = c(1, 2))
  if (is.null(date)) date <- seq_along(ysim)
  
  plot(date, yobs, type = "l", xlab = "") # , ylab = "Discharges [mm/day]"
  lines(date, ysim, col = 2)

  info <- GOF(yobs, ysim) %>% as.list()
  gof <- sprintf("NSE = %.3f, R2 = %.3f, Bias(%%) = %.1f%%", info$NSE, info$R2, info$Bias_perc * 100)
  legend("topleft", legend = c("Observations", "Simulations"), col = c(1, 2), lty = 1, bty = "n")

  range1 <- range(yobs, na.rm = TRUE)
  range2 <- range(ysim, na.rm = TRUE)
  lims <- c(min(range1[1], range2[1]), max(range1[2], range2[2]))

  plot(yobs, ysim,
    xlim = lims, ylim = lims, pch = 20,
    xlab = "Observations", ylab = "Simulations",
    col = scales::alpha("black", 0.4)
  )
  legend("topright", legend = gof, col = "red", bty = "n")
  abline(a = 0, b = 1, col = "red")
}

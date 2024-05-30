### XAJ default parameter ------------------------------------------------------
param.names <- c("KC", "IM", "WUM", "WLM", "WDM", "C", "B", "SM", "EX", "KI", "KG", "CI", "CG", "N", "NK")
XAJ.param.range <- data.frame(
    lower = c(0.20, 0.00, 5.0, 10.0, 10.0, 0.05, 0.1, 10.0, 0.50, 0.01, 0.01, 0.50, 0.95, 0.1, 1.0),
    upper = c(1.50, 0.05, 20., 90.0, 60.0, 0.20, 0.6, 60.0, 2.00, 0.70, 0.70, 0.90, 0.998, 5.0, 6.0),
    row.names = param.names
)
XAJ.param.range$par <- with(XAJ.param.range, lower + upper) / 2

d_par <- XAJ.param.range
lb <- d_par$lower
ub <- d_par$upper
par0 <- d_par$par
# ------------------------------------------------------------------------------


XAJ_goal <- function(par, Qobs, prcp, ET0, area, date = NULL, dt = 24, index = "KGE") {
  r = VIC5::XAJ(prcp, ET0, par, area = area, dt = dt)

  Q = r$Q
  if (length(Qobs) < length(prcp)) {
    Q = aggregate(Q, by = list(make_monthdate(date)), mean)[, 2] # convert to monthly
  }

  if (index == "KGE") {
    -KGE(Qobs, Q)
  } else if (index == "NSE") {
    -NSE(Qobs, Q)
  }
}

#' XAJ model Parameter calibration
#'
#' @inheritParams VIC5::XAJ
#'
#' @param Qobs Observed Total runoff, (m^3/s)
#' @param prcp Precipitation (mm/d)
#' @param ET0 Pan evaporation or potential evapotranspiration (mm/d)
#' @param date (optional) corresponding date of `Qobs`
#'
#' @param area basin area (km^2).
#' @param dt time step (hour)
#'
#' @param index KGE or NSE
#' @param seed (can be ignored) starting number of random number generator,
#' see [base::set.seed()] for details.
#' This parameter is to make sure optimization result is same in different tries.
#' @param ... ignored
#'
#' @export
XAJ_calib <- function(Qobs, prcp, ET0, area, dt = 24, date = NULL,
  maxn = 1000, index = "KGE", seed = 1, ...)
{
  set.seed(1)
  l_opt = rtop::sceua(XAJ_goal, par0, lb, ub,
    Qobs = Qobs, prcp = prcp, ET0 = ET0, area = area, dt = dt, date = date,
    maxn = maxn, index=index)
  # The optimized best parameter
  opt = l_opt$par %>% set_names(rownames(d_par))

  r = VIC5::XAJ(prcp, ET0, opt, area = area, dt = dt)

  # output data
  data = data.table(prcp, ET0, Qobs, Qsim = r$Q)
  if (!is.null(date)) {
    data %<>% cbind(date, .)
    r %<>% cbind(date, .)
  }
  # Qsim = aggregate(Qsim, by = list(make_monthdate(dates)), mean)[, 2]
  gof = GOF(Qobs, data$Qsim)
  list(data = data, gof = gof, model = listk(par = opt, area, dt), intermediate = r)
}

#' XAJ_predict
#'
#' @param model XAJ model returned by [XAJ_calib()]
#' @param newdata A data.frame, with the column of `prcp`, `ET0`, `date` (optional).
#'
#' @export
XAJ_predict <- function(model, newdata) {
  res = VIC5::XAJ(newdata$prcp, newdata$ET0, model$par, area = model$area, dt = model$dt)
  # Qsim = res$Q
  cbind(newdata, Qsim = res$Q)
  # date = newdata$date
  # if (!is.null(date)) {
  #   res %<>% cbind(date, .)
  # }
  # res
}

#' plot_calib
#'
#' @param data with the columns of `date`, `prcp`, `Qsim`,  `Qobs`
#' @param ... others to plot precipitation
#' 
#' @export
plot_calib <- function(data, ...) {
  date = data$date
  if (is.null(date)) date = 1:nrow(data)

  gof = GOF(data$Qobs, data$Qsim) %>% as.list()
  # TODO: update plot_runoff
  str_gof = do.call(sprintf,
    c("NSE=%.3f, KGE=%.3f, RMSE=%.1f, \nBias=%.1f, R2=%.3f, n=%d",
    gof[c("NSE", "KGE", "RMSE", "Bias", "R2", "n_sim")] %>% as.list()))
  
  # write_fig({
    par(mar = c(3.5, 3.5, 2, 1), mgp = c(1.8, 0.6, 0))
    plot(date, data$Qobs, type = "l", col = "black",
         xlab = "Date", ylab = expression(Q * "(" * m^3 * ")"))
    lines(date, data$Qsim, col = "red")

    par(new = TRUE)
    plot(date, data$prcp,
        ylim = c(140, -25), type = "h", col = "blue", xaxt = "n", yaxt = "n",
        ylab = "", xlab = "", ...
    )

    legend("topleft", legend = c("obs", "sim", "prcp"),
           col = c("black", "red", "blue"), lty = 1, ncol = 3)
    legend("topright", str_gof, text.col = "red", bty = "n")
  # }, file)
}

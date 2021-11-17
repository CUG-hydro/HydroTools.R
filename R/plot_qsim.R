
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

    # smoothScatter
    plot(yobs, ysim,
        xlim = lims, ylim = lims, pch = 20, 
        xlab = "Observations", ylab = "Simulations",
        col = scales::alpha("black", 0.4)
    )
    legend("topright", legend = gof, col = "red", bty = "n")
    abline(a = 0, b = 1, col = "red")
}

library(hydroTools)
library(Ipaper) # for `write_fig` function, `remotes::install_github("rpkgs/Ipaper")`

infile = system.file("extdata/XAJ_input_Hejin_2012-2018.rda", package = "hydroTools")
load(infile)

df = input$data
res <- XAJ_calib(df$Qobs, df$prcp, df$ET0, date = df$date,
    input$area, dt = 24, maxn = 1000)

write_fig({
  par(family = "rTimes")
  plot_calib(res$data, main = input$site)
}, "XAJ_Hejin.pdf")

XAJ_predict(res$model, df)

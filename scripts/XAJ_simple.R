library(hydroTools)
library(XAJ)

load("XAJ_input_河津_2012-2018.rda")

df = input$data
res <- XAJ_calib(df$Qobs, df$prcp, df$ET0, date = df$date,
                        input$area, dt = 24, maxn = 1000)
plot_calib(res$data, main = input$site)

XAJ_predict(res$model, df)

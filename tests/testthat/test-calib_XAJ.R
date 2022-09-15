
test_that("XAJ_calib works", {

  infile = system.file("extdata/XAJ_input_Hejin_2012-2018.rda", package = "hydroTools")
  load(infile)

  df = input$data
  res <- XAJ_calib(df$Qobs, df$prcp, df$ET0, date = df$date,
      input$area, dt = 24, maxn = 200)

  # print(str(res), 2)
  df_q = res$data[, c("date", "Qobs", "Qsim")]
  df_prcp = res$data[, c("date", "prcp")]
  # plot_calib(res$data, main = input$site)
  expect_silent(plot_runoff(df_q, df_prcp, ylim2 = c(60, 0), legend.position = c(1, 0.6)))

  res_valid = XAJ_predict(res$model, newdata = df)
  head(res_valid)
  expect_equal(colnames(res_valid), c("date", "Qobs", "prcp", "ET0", "Qsim"))
})

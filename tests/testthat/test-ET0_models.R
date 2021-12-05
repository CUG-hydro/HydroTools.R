test_that("ET0 models work", {
    Ra = cal_Rn(lat = 30, J = 1, ssd = 10, RH = 70, Tmin = 20, Tmax = 30, cld = 0.2)

    et0_pt = ET0_PT72(250, 25, 1)
    et0_pm = ET0_PM93(250, 25, 1, wind = 2)
    et0_pm98 = ET0_PM98(250, 25, 1, wind = 2)
    
    expect_true(et0_pm$ET0_PM93 > et0_pt$ET0_PM72/1.26)
    expect_true(et0_pm98$Eeq_pm < et0_pm98$Eeq)
})

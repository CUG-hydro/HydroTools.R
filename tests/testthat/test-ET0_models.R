test_that("ET0 models work", {
    Ra = cal_Rn(lat = 30, J = 1, ssd = 10, RH = 70, Tmin = 20, Tmax = 30, cld = 0.2)

    et0_pt72 = ET0_PT72(250, 25, D = 1)
    et0_pm54 = ET0_Penman54(250, 25, D = 1, wind = 2)
    et0_fao98 = ET0_FAO98(250, 25, D = 1, wind = 2)
    
    expect_true(et0_pm54$ET0 > et0_pt72$ET0 / 1.26)
    expect_true(et0_fao98$Eeq < et0_pm54$Eeq)
})

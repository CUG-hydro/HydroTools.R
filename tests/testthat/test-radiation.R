
test_that("radiation function works", {
    # ssd
    ssd = cal_ssd(20, 1)
    expect_true(ssd > 10 && ssd < 24)
    expect_true(cal_ssd(20, 1) < cal_ssd(20, 150))
    
    # Rsi_toa
    Rsi_toa = cal_Rsi_toa(20, 1)
    Rsi = cal_Rsi(20, 1)
    expect_equal(Rsi_toa * 0.75, Rsi$Rsi) # `[MJ m-2 day-1]`
})


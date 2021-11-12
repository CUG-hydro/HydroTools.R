
test_that("radiation function works", {
    # ssd
    ssd = cal_ssd(20, 1)
    expect_true(ssd > 10 && ssd < 24)
    expect_true(cal_ssd(20, 1) < cal_ssd(20, 150))
    
    # Ra
    Ra = cal_Ra(20, 1)
    Rs = cal_Rs(20, 1)
    expect_equal(Ra * 0.75, Rs) # `[MJ m-2 day-1]`
})


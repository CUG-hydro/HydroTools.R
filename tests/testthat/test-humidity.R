
test_that("humidity function works", {
    RH = 90
    p = atm
    t = 30
    q = RH2q(RH, t, p)
    RH2 = q2RH(q, t, p)
    e = q2ea(q, p)

    max_MAE(q2RH(q, t, p), RH, 1e-7)

    # If no RH, it is `ea = es(Tmin)`
    e1 = cal_ea(20) 
    expect_equal(e1, cal_ea(20, 30))

    e2 = cal_ea(20, 25, 60)
    # expect_equal(e2, cal_ea(20, 25, 60))
    # expect_true(e2 < cal_ea(20, 25, RH_mean = 80))

    # e3 = cal_ea(20, RH_max = 80)
    # expect_true(e1 > e3)
    # delta_es = delta_es(20)
})

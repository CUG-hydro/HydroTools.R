maxError <- function(yobs, ysim, maxE = 0.005){
    re = abs(ysim - yobs) #/yobs
    # print(re)
    expect_lte(re, maxE)
}

test_that("humidity function works", {
    RH = 90
    p = atm
    t = 30
    q = RH2q(RH, p, t)
    RH2 = q2RH(q, p, t)
    e = q2ea(q, p)

    maxError(q2RH(q, p, t), RH, 1e-7)

    e1 = cal_ea(20)
    expect_equal(e1, cal_ea(20, 30))

    e2 = cal_ea(20, 25, 60, 80)
    expect_equal(e2, cal_ea(20, 25, 60, 80, 65))
    expect_true(e2 < cal_ea(20, 25, RH_mean = 80))

    e3 = cal_ea(20, RH_max = 80)
    expect_true(e1 > e3)

    delta_es = delta_es(20)
})

yobs <- rnorm(100)
#' ysim = yobs + rnorm(100)/4
#' GOF(yobs, ysim)

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
})

maxError <- function(yobs, ysim, maxE = 0.005){
    re = abs(ysim - yobs)/yobs
    # print(re)
    expect_lte(re, maxE)
}

test_that("SPAW function works", {
    S = 20/100; C = 20/100; OM = 2.5

    maxError(0.137, wilting_point(S, C, OM))
    maxError(0.321, field_capacity(S, C, OM))
    maxError(0.482, saturated_mois(S, C, OM))
    maxError(12.19/24, Ksat(S, C, OM))
    maxError(1.37 , rho_norm(S, C, OM))
    maxError(0.145 , wilting_point_salinity(S, C, OM, 3))
})

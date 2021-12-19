test_that("SPAW function works", {
    S = 20/100; C = 20/100; OM = 2.5

    max_Bias(0.137, wilting_point(S, C, OM))
    max_Bias(0.321, field_capacity(S, C, OM))
    max_Bias(0.482, saturated_mois(S, C, OM))
    max_Bias(12.19/24, Ksat(S, C, OM))
    max_Bias(1.37 , rho_norm(S, C, OM))
    max_Bias(0.145 , wilting_point_salinity(S, C, OM, 3))
})

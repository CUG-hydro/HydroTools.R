
test_that("plot works", {
    Y_obs <- c(1, 3, 4, 2, NA, 6)
    Y_sim <- c(2, 2, 4, 1, NA, 5)
    
    expect_true({
        plot_gof(Y_obs, Y_sim)
        TRUE
    })
    # expect_equal(dim(UH$UH), c(96, 17))
})

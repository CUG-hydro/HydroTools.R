
# absolute value
max_MAE <- function(yobs, ysim, maxE = 0.005) {
    re <- abs(ysim - yobs) # /yobs
    # print(re)
    expect_lte(re, maxE)
}

# percentage
max_Bias <- function(yobs, ysim, maxE = 0.005) {
    re <- abs(ysim - yobs) / yobs
    # print(re)
    expect_lte(re, maxE)
}

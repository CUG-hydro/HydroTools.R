n <- 1e3
t <- seq(-50, 100, length.out = n)
rh <- seq(0, 100, length.out = length(t))

## bug exist
i <- 88
heat_index(t[i], rh[i])
heat_index_vec(t[i], rh[i])
# heat_index_julia(t[i], rh[i])

# julia_setup()
r <- heat_index(t, rh)
r_vec <- heat_index_vec(t, rh)
# r_jl <- heat_index_julia(t, rh)

# add tests for julia NA values
# r1 <- heat_index_julia(t*NA, rh*NA)

all.equal(r, r_vec)
# all.equal(r_jl, r_vec)

# microbenchmark::microbenchmark(
#     r2 <- heat_index(t, rh),
#     r1 <- heat_index_julia(t, rh)
#     # r_vec <- heat_index_vec(t, rh)
# )

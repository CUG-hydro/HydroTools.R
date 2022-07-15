test_that("heat_index_vec works", {
  n <- 1e3
  t <- seq(-10, 50, length.out = n)
  rh <- seq(0.1, 1, length.out = length(t))

  r0 <- heat_index(t, rh)
  r_vec <- heat_index_vec(t, rh)
  expect_equal(r0, r_vec)
})

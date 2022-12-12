test_that("cal_angle and cal_azimuth works", {
  p1 <- c(1, 0)
  p2 <- c(0, 1)

  expect_equal(cal_angle(c(1, 0), c(2, 0)), 0)
  expect_equal(cal_angle(c(0, 1), c(0, 2)), 0)
  expect_equal(cal_angle(p1, p2), 90) # counterclockwise 90
  expect_equal(cal_angle(p2, p1), -90) # clockwise -90

  expect_equal(cal_azimuth(c(1, 0)), 90)
  expect_equal(cal_azimuth(c(0, 1)), 0)
  expect_equal(cal_azimuth(c(-1, 0)), 270)
  expect_equal(cal_azimuth(c(0, -1)), 180)
})


test_that("cal_angle works for matrix", {
  p_east <- c(1, 0)
  p_west <- c(-1, 0)
  p_north <- c(0, 1)
  p_south <- c(0, -1)
  
  mat1 = rbind(p_east, p_west, p_north, p_south)
  mat2 = rbind(p_west, p_north, p_south, p_east)

  expect_equal(
    setNames(cal_azimuth(mat1), NULL), 
    c(90, 270, 0, 180)
  )
  expect_equal(
    setNames(cal_angle(mat1, mat2, clockwise = 1), NULL),
    c(-180, 90, 180, 90)
  )
})

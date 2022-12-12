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

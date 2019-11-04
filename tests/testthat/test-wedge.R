context("Results from the wedge function")
library(treenomial)


test_that("Check wedging two cherries result", {
  cherry <- matrix(c(0, 1, 0, 0, 1, 0), nrow = 2) # x^2 + y
  cherryWedged <- matrix(c(0,1,1,0,0,0,0,2,0,0,0,0,1,0,0), nrow = 3) # x^4 + 2x^2y + y^2 + y
  expect_equal(wedge(cherry, cherry), cherryWedged)
})

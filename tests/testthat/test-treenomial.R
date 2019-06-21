library(ape)
library(apTreeshape)
library(Matrix)

test_that("Check wedging two cherries result", {
  cherry <- sparseMatrix(c(1, 2), c(3, 1), x = c(1, 1))
  cherryWedged <- sparseMatrix(c(3, 2, 2, 1), c(1, 1, 3, 5), x = c(1, 1, 2, 1))
  expect_equal(wedge(cherry, cherry), cherryWedged)
})

test_that("Test consistency of coefficientMatrix runs tips = 500", {
  numTips <- 500
  tree <- rtree(n = numTips)
  expect_equal(coefficientMatrix(tree), coefficientMatrix(tree))
})

# test with labelled data

# test on plotting

# test on error conditions

# test on all trees

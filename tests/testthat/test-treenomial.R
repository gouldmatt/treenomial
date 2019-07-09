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

test_that("test complex all trees up to 16 tips", {
  wadNum <- c(1,1,1,2,3,6, 11, 23, 46, 98, 207, 451,
              983, 2179, 4850, 10905)


  allTrees <- allTrees(16, complex = TRUE)
  # verify against wad number
  allTrees <- unique(allTrees)

  expect_equal(wadNum,lengths(allTrees))
})



# test with labelled data

# test on plotting

# test on error conditions

# R0 tests

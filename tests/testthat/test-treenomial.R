library(apTreeshape)
library(Matrix)

#####
# tests to add
#   ===> error cases
#   ===> alltrees
#   ===> wedge
#   ===>
#


##### tests involving construction of coefficient matrices ######
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

##### tests involving metric on coefficient matrices ######
test_that("test distance matrix of same trees is zero", {
  tree <- list(rtree(10))
  identicalForest <- rep(tree,10)
  distanceMatrix <- coefficientDist(coefficientMatrix(identicalForest))
  expect_true(all(distanceMatrix == 0))
})

test_that("ensure distance matrix is symmetric", {
  numTrees <- 100
  numTips <- 20

  pdaTrees <- rtreeshape(numTrees, tip.number = numTips, model = "pda")
  yuleTrees <- rtreeshape(numTrees, tip.number = numTips, model = "yule")
  aldousTrees <- rtreeshape(numTrees, tip.number = numTips, model = "aldous")
  biasedTrees <- rtreeshape(numTrees, tip.number = numTips, model = "biased")

  modelTrees <- list(pda = pdaTrees, yule = yuleTrees,aldous =aldousTrees, biased = biasedTrees)
  modelTrees <- lapply(unlist(modelTrees, recursive=FALSE), as.phylo)

  coeffMats <- coefficientMatrix(modelTrees)

  distanceMatrix <- coefficientDist(coeffMats)

  expect_true(all(distanceMatrix == t(distanceMatrix)))
})


context("Distances between coefficient matrices")
library(treenomial)
library(ape)


## tests on polyToDistMat ##
test_that("test that different distance method have correct result (small manual check)", {


  coeffs <- allTrees(6)[[6]]

  dLogDiff <- polyToDistMat(coeffs, method = "logDiff")

  dWLogDiff <- polyToDistMat(coeffs, method = "wLogDiff")

  testWlogL1 <- function(coeffA,coeffB){
    logDiffMat <- rowSums(log(1 + abs(coeffA-coeffB)))

    weightVect <- c(1,(1:(nrow(coeffA)-1)))^(-2)

    sum(logDiffMat*weightVect)

  }

  expect_equal(dWLogDiff[2,4], testWlogL1(coeffs[[2]], coeffs[[4]]))

})

## tests on treeToDistMat ##
test_that("test distance matrix of same trees is zero for each method", {
  tree <- list(rtree(10))
  identicalForest <- rep(tree,10)

  distanceMatrix <- treeToDistMat(identicalForest, method = "logDiff", type = "complex")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- treeToDistMat(identicalForest, method = "logDiff")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- treeToDistMat(identicalForest, method = "wLogDiff")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- treeToDistMat(identicalForest, method = "pa")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- treeToDistMat(identicalForest, method = "ap")
  expect_true(all(distanceMatrix == 0))


})

test_that("testing naming carry through", {
  trees <- c(rmtree(2,10),rmtree(2,131),rtree(2),rtree(3))
  names(trees) <- c(rep(c("smaller","larger"), each = 2),"twoTip","threeTip")
  d <- treeToDistMat(trees)
  expect_equal(rownames(d),names(trees))
  expect_equal(colnames(d),names(trees))
})


test_that("ensure distance matrix is symmetric", {

  distanceMatrix <- treeToDistMat(rmtree(100,20))

  expect_true(all(distanceMatrix == t(distanceMatrix)))
})


## tests on treeToDistMatPlot ##
# test_that("ensure correct max/min dist value is being used", {
#   d <- matrix(data = 7, nrow = 100, ncol = 100)
#   diag(d) <- 0
#   d <- d + t(d)
#
#   d[50,49] <- 42
#   d[50,36] <- 12
#
#   res <- treeToDistMatPlot(rmtree(100,2), d,50, "min", returnNearestInfo = T)
#   expect_equal(res$distance, 12)
#
#   res <- treeToDistMatPlot(rmtree(100,2), d,50, "max", returnNearestInfo = T)
#   expect_equal(res$distance, 42)
#
#   forest <- rmtree(100, 15)
#
#   d <- treeToDistMat(forest)
#
#   res <- treeToDistMatPlot(forest, d,50, "min", returnNearestInfo = T)
#
#   expect_equal(res$distance, min(d[50,-50]))
#
#
#   res <- treeToDistMatPlot(forest, d,25, "max", returnNearestInfo = T)
#
#   expect_equal(res$distance, max(d[25,-25]))
#
#   forest <- rmtree(100, 5)
#
#   d <- treeToDistMat(forest, method = "pa")
#
#   res <- treeToDistMatPlot(forest, d,50, "min", returnNearestInfo = T)
#
#   expect_equal(res$distance, min(d[50,-50]))
#
#
#   res <- treeToDistMatPlot(forest, d,25, "max", returnNearestInfo = T)
#
#   expect_equal(res$distance, max(d[25,-25]))
#
#
#
#
# })

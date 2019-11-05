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


## tests on plotExtremeTrees ##
test_that("ensure correct min/max trees are being found", {
  forestTen <- rmtree(100,10)
  forestSixty <- rmtree(42,60)
  threeTip <- rtree(3)

  minTrees <- plotExtremeTrees(threeTip, c(forestSixty,threeTip,forestTen), 2)

  expect_equal(minTrees[[1]]$distance, 0)
  expect_equal(treeDist(minTrees[[1]]$tree,threeTip), 0)

  maxTrees <- plotExtremeTrees(threeTip, c(forestSixty,threeTip,forestTen), 2, comparison = "max")

  expect_equal(maxTrees[[1]]$distance, treeDist(threeTip,maxTrees[[1]]$tree))
})


## other tests ##
test_that("tests with single and lists arguments", {
  forestTen <- rmtree(100,10)
  forestSixty <- rmtree(42,60)
  threeTip <- rtree(3)

  expect_silent(treeDist(threeTip,threeTip))
  expect_silent(treeDist(threeTip,forestSixty))
  expect_error(treeDist(list(threeTip),forestSixty))
  expect_error(treeDist(list(threeTip),list(forestSixty)))

  coeffsSixty <- treeToPoly(forestSixty)
  coeffsThree <- treeToPoly(threeTip)

  expect_silent(polyDist(coeffsThree,coeffsThree))
  expect_silent(polyDist(coeffsThree,coeffsSixty))
  expect_error(polyDist(list(coeffsThree),coeffsSixty))
  expect_error(polyDist(list(coeffsThree),list(coeffsSixty)))

  expect_silent(plotExtremeTrees(threeTip, threeTip, 1))
  expect_silent(plotExtremeTrees(threeTip, forestSixty, 1))
  expect_error(plotExtremeTrees(list(threeTip),forestSixty,1))
  expect_error(plotExtremeTrees(list(threeTip),list(forestSixty),1))
})



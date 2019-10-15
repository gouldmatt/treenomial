library(apTreeshape)
library(Matrix)
library(ape)
library(pracma)

#####
# tests to add
#   ===> error cases
#   ===> alltrees
#   ===> wedge
#   ===>
#

# test to generate warnings

# test incorrect inputs

# test different list inputs

# test repeated calls of single

# test for alignments

# test for wedging phylo objects

# test calculating distance for single and lists of trees

# test that the formats coming out of functions is correct like matrix or list or what is expected


############## tests on wedge
test_that("Check wedging two cherries result", {
  cherry <- sparseMatrix(c(1, 2), c(3, 1), x = c(1, 1))
  cherryWedged <- sparseMatrix(c(3, 2, 2, 1), c(1, 1, 3, 5), x = c(1, 1, 2, 1))
  expect_equal(wedge(cherry, cherry, type = "real"), cherryWedged)
})


############## tests on coeffMatrix
test_that("Test consistency of coeffMatrix runs tips = 500", {
  numTips <- 500
  tree <- rtree(n = numTips)
  expect_equal(coeffMatrix(tree), coeffMatrix(tree))
})

############## tests on coeffDist
test_that("test distance matrix of same trees is zero for each method", {
  tree <- list(rtree(10))
  identicalForest <- rep(tree,10)

  distanceMatrix <- phyloDist(identicalForest, method = "logL1", type = "complex")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- phyloDist(identicalForest, method = "logL1")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- phyloDist(identicalForest, method = "wLogL1")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- phyloDist(identicalForest, method = "b")
  expect_true(all(distanceMatrix == 0))

  distanceMatrix <- phyloDist(identicalForest, method = "c")
  expect_true(all(distanceMatrix == 0))


})

test_that("ensure distance matrix is symmetric", {

  distanceMatrix <- phyloDist(rmtree(100,20))

  expect_true(all(distanceMatrix == t(distanceMatrix)))
})


############## tests on allBinaryTreeShapes
test_that("ensure that number of trees from allBinaryTreeShapes match Wedderburn-Etherington numbers up to 13 tips", {
  wadNum <- c(1,1,1,2,3,6, 11, 23, 46, 98, 207, 451,
              983)

              # 2179, 4850, 10905)

  # check real version
  allTrees <- allBinaryTreeShapes(13, type = "real")
  allTrees <- unique(allTrees)
  expect_equal(wadNum,lengths(allTrees))

  # check complex version
  allTrees <- allBinaryTreeShapes(13, type = "complex")
  allTrees <- unique(allTrees)
  expect_equal(wadNum,lengths(allTrees))

  allTrees <- allBinaryTreeShapes(13, type = "phylo")
  allTrees <- unique(allTrees)
  expect_equal(wadNum,lengths(allTrees))
})


############## tests on distPlot
test_that("ensure correct max/min dist value is being used", {
  d <- matrix(data = 7, nrow = 100, ncol = 100)
  diag(d) <- 0
  d <- d + t(d)

  d[50,49] <- 42
  d[50,36] <- 12

  res <- phyloDistPlot(rmtree(100,2), d,50, "min", returnNearestInfo = T)
  expect_equal(res$distance, 12)

  res <- phyloDistPlot(rmtree(100,2), d,50, "max", returnNearestInfo = T)
  expect_equal(res$distance, 42)

  forest <- rmtree(100, 15)

  d <- phyloDist(forest)

  res <- phyloDistPlot(forest, d,50, "min", returnNearestInfo = T)

  expect_equal(res$distance, min(d[50,-50]))


  res <- phyloDistPlot(forest, d,25, "max", returnNearestInfo = T)

  expect_equal(res$distance, max(d[25,-25]))

  forest <- rmtree(100, 5)

  d <- phyloDist(forest, method = "b")

  res <- phyloDistPlot(forest, d,50, "min", returnNearestInfo = T)

  expect_equal(res$distance, min(d[50,-50]))


  res <- phyloDistPlot(forest, d,25, "max", returnNearestInfo = T)

  expect_equal(res$distance, max(d[25,-25]))




})





context("Building coefficient matrices")
library(treenomial)
library(ape)

test_that("Test consistency of treeToLattice runs tips = 120", {
  numTips <- 120
  tree <- rtree(n = numTips)
  expect_equal(treeToLattice(tree, numThreads = 0), treeToLattice(tree, numThreads = 0))
})

test_that("Test lattice with two tip tree", {
  numTips <- 2
  tree <- rtree(n = numTips)
  expect_silent(treeToLattice(tree, numThreads = 0))
})

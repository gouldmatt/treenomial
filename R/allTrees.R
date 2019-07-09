#' Calculate the coefficient matrices of all trees up to n tips
#'
#' Return coefficient matrices for all possible full binary trees with 2 to n tips.
#' @param numtips the max number of tips to find all coefficient matrices for
#' @return list containing all possible coefficient matrices
#' @import gtools
#' @examples
#' # example goes here
#' @export
allTrees <- function(numTips, complex = FALSE) {
  wadNum <- c(1,1,1,2,3,6, 11, 23, 46, 98, 207, 451,
              983, 2179, 4850, 10905, 24631, 56011, 127912,
              293547, 676157, 1563372, 3626149, 8436379, 19680277,
              46026618, 107890609)

  nTipsTrees <- vector("list", length = numTips)
  if(complex){
    nTipsTrees[[1]] <- list(c(0,1))
  } else {
    nTipsTrees[[1]] <- list(sparseMatrix(1, 2, x = 1))
  }
  currTips <- 2

  while (currTips != (numTips + 1)) {

    # determine the possible combinations of smaller trees that could be used to generate a tree with currtipss
    treesToWedge <- subset_sum(currTips)

    # utilize all possible previous combinations of the smaller trees
    treeCount <- 0
    for (i in 1:nrow(treesToWedge)) {
      lengthA <- length(nTipsTrees[[treesToWedge[i, 1]]])
      lengthB <- length(nTipsTrees[[treesToWedge[i, 2]]])

      treeOptions <- expand.grid(1:lengthA, 1:lengthB)

      # remove order swapped trees
      if (treesToWedge[i, 2] == treesToWedge[i, 1]) {
        treeOptions <- treeOptions[treeOptions[, 2] >= treeOptions[, 1], ]
      }

      # perform the wedge
      for (j in 1:nrow(treeOptions)) {
        treeCount <- treeCount + 1
        tree1 <- nTipsTrees[[treesToWedge[i, 1]]][[treeOptions[j, 1]]]
        tree2 <- nTipsTrees[[treesToWedge[i, 2]]][[treeOptions[j, 2]]]
        if(complex){
          nTipsTrees[[currTips]][[treeCount]] <- wedgeC(nTipsTrees[[treesToWedge[i,1]]][[treeOptions[j,1]]], nTipsTrees[[treesToWedge[i,2]]][[treeOptions[j,2]]])
        } else {
          nTipsTrees[[currTips]] <- append(nTipsTrees[[currTips]], wedge(tree1, tree2))
        }
      }
    }

    currTips <- currTips + 1
  }
  return(nTipsTrees)
}

subset_sum <- function(targetTips) {
  allComb <- combinations(n = targetTips - 1, r = 2, repeats.allowed = TRUE)
  res <- allComb[rowSums(allComb) == targetTips, ]
  if (!is.matrix(res)) {
    res <- t(res)
  }
  return(res)
}

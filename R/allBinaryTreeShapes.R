#' Calculate all full unordered binary trees up to n tips
#'
#' Return complex coefficient matrices, real coefficient matrices or phylo objects for all possible full binary trees up to n tips.
#' @param numtips max number of tips
#' @param type "real", "complex" or "phylo"
#' @importFrom gtools combinations
#' @return list of lists containing the trees
#' @examples
#'
#' # phylo type example
# allBinSix <- allBinaryTreeShapes(6, type = "phylo")[[6]]
# par(mfrow=c(1,6))
# lapply(allBinSix, function(t) plot.phylo(ladderize(t), direction = "downwards", show.tip.label = FALSE))
#'
#' @export
allBinaryTreeShapes <- function(numTips, type = "real", progressBar = FALSE, startingPoint) {
  # wadNum <- c(1,1,1,2,3,6, 11, 23, 46, 98, 207, 451,
  #             983, 2179, 4850, 10905, 24631, 56011, 127912,
  #             293547, 676157, 1563372, 3626149, 8436379, 19680277,
  #             46026618, 107890609)

  if(missing(startingPoint)){

    nTipsTrees <- vector("list", length = numTips)
    nTipsTrees <- lapply(1:numTips, function(i){ nTipsTrees[[i]] <- list()})

    if(type == "complex"){
      nTipsTrees[[1]] <- list(c(0,1))
    } else if(type == "real") {
      nTipsTrees[[1]] <- list(sparseMatrix(1, 2, x = 1))
    } else if(type == "phylo"){
      nTipsTrees[[1]] <- "leaf"
    } else {
      # error
    }
    currTips <- 2

    if(progressBar){
      # progress bar
      pb <- txtProgressBar(min = 0, max = numTips + 1, style = 3)
    }
  } else {
    currTips <- length(startingPoint) + 1
    nTipsTrees <- vector("list", length = numTips)
    nTipsTrees[1:length(startingPoint)] <- startingPoint

    if(progressBar){
      # progress bar
      pb <- txtProgressBar(min = currTips, max = numTips + 1, style = 3)
    }
  }

  while (currTips != (numTips + 1)) {

    # determine the possible combinations of smaller trees that could be used to generate a tree with currtips
    treesToWedge <- tipOptions(currTips)

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

        nTipsTrees[[currTips]][[treeCount]] <- wedge(nTipsTrees[[treesToWedge[i,1]]][[treeOptions[j,1]]], nTipsTrees[[treesToWedge[i,2]]][[treeOptions[j,2]]], type = type)

      }
    }

    currTips <- currTips + 1

    if(progressBar){
      setTxtProgressBar(pb, currTips)
    }
  }

  if(progressBar){
    close(pb)
  }

  return(nTipsTrees)
}

tipOptions <- function(targetTips) {
  allComb <- combinations(n = targetTips - 1, r = 2, repeats.allowed = TRUE)
  res <- allComb[rowSums(allComb) == targetTips, ]
  if (!is.matrix(res)) {
    res <- t(res)
  }
  return(res)
}

#' Converts trees (objects of class "phylo") to real coefficient matrices or complex coefficient vectors
#'
#' @param trees an object or list of class "phylo"
#' @param type "real" to return a coefficient matrix of the real polynomial or "complex" for the complex coefficietn vector of the real polynomial with y = 1 + i
#' @param progressBar add a progress bar
#' @param cl A cluster object created by makeCluster, from the parallel package
#' @return The resulting coefficient matrix or matrices.
#' @import Matrix
#' @import pbapply
#' @importFrom ape as.phylo
#' @examples
#' require(ape)
#' # generate a tree
#' tree <- rtree(n = 30, rooted = TRUE)
#' # a real coefficient matrix
#' coefficientMatrix(tree)
#'
#' # complex coefficient vector for the tree
#' coefficientMatrix(tree, complex = TRUE)
#'
#' # for a list of trees with a progress bar
#' coeffMatrix(rmtree(10, 30), progressBar = TRUE)
#'
#' # use the parallel package to speed up conversion for a large list of phylo objects
#' require(parallel)
#' cl <- makeCluster(detectCores())
#' largeForest <- coeffMatrix(rmtree(400, 500), progressBar = TRUE, cl = cl)
#' stopCluster(cl)
#' @export
coeffMatrix <- function(trees, type = "real", progressBar = FALSE, cl = NULL) {

  # check input format
  if (class(trees) == "phylo") {
    singleTree <- TRUE
  } else if (class(trees) == "list" || class(trees) == "multiPhylo") {
    tryCatch({
      trees <- lapply(trees, as.phylo)
      singleTree <- FALSE
    }, error = function(e) {
      stop("incorrect input format, trees must be phylo or list of phylo objects")
    })
  } else {
    tryCatch({
      trees <- as.phylo(trees)
      singleTree <- TRUE
    }, error = function(e) {
      stop("incorrect input format, trees must be phylo or list of phylo objects")
    })
  }

  # set the max times the progress bar is updated and check if user wants a progress bar
  pboptions(nout = 50, type = ifelse(progressBar, "timer", "none"))


  ### tipLabel option
  if(type == "tipLabel"){
    if (singleTree) {
      if (progressBar) {
        warning("progress bar only available for list inputs")
      }

      coefficientMat <- tipLabelCoeffMat(trees)
    } else {
      coefficientMat <- vector("list", length = length(trees))
      coefficientMat <- pblapply(trees, tipLabelCoeffMat)
      names(coefficientMat) <- names(trees)
    }

    return(coefficientMat)
  }


  # call appropiately for the number and size of the input trees
  if (singleTree) {
    if (progressBar) {
      warning("progress bar only available for list inputs")
    }

    coefficientMat <- singleCoeffMat(trees, type = type)
  } else {
    coefficientMat <- vector("list", length = length(trees))
    coefficientMat <- pblapply(trees, singleCoeffMat, type = type, cl = cl)
    names(coefficientMat) <- names(trees)
  }

  return(coefficientMat)
}

singleCoeffMat <- function(tree, type = "real") {

  # construct the child matrix of the tree
  nnodes <- tree$Nnode # number of internal nodes
  nTips <- length(tree$tip.label)
  nodeids <- (nTips + 1):(nTips + nnodes) # label starting after last tip to the total nodes
  ChildMatrix <- t(sapply(nodeids, function(x) tree$edge[tree$edge[, 1] == x, 2]))
  ChildMatrix <- cbind(nodeids, ChildMatrix) #  first col are the internal nodes

  # determine the order of wedges needed to generate the trees coefficient matrix
  numTips <- ChildMatrix[1, 1] - 1
  wedgeOrder <- ifelse(unique(rev(as.vector(t(ChildMatrix)))) <= numTips, "0", "1")

  done <- F

  # intialize first two entries for leaves and cherries
  subCoeffMats <- list()

  if (type == "complex") {
    subCoeffMats[["0"]] <- c(0, 1) # leaf
    subCoeffMats[["001"]] <- c(1i + 1, 0, 1) # cherry
  } else {
    subCoeffMats[["0"]] <- sparseMatrix(1, 2, x = 1) # leaf
    subCoeffMats[["001"]] <- sparseMatrix(c(1, 2), c(3, 1), x = c(1, 1)) # cherry
  }

  # continue while there is still stuff left to wedge
  while (!(length(wedgeOrder) == 1)) {

    # find the next two operands to wedge
    j <- 3
    while (TRUE) {
      if (wedgeOrder[j] == 1) {
        operand1 <- wedgeOrder[j - 2]
        operand2 <- wedgeOrder[j - 1]
        break
      }
      j <- j + 1
    }

    # calculate the pattern
    subPattern <- paste(operand1, operand2, "1", sep = "")
    old <- paste(wedgeOrder, collapse = " ", sep = ",")
    oldTest <- paste(wedgeOrder, collapse = "")

    # calculate wedge
    subCoeffMats[[as.character(subPattern)]] <- wedge(subCoeffMats[[operand1]], subCoeffMats[[operand2]], type = type)

    # replace subsequent operations of the same wedge
    wedgeOrder <- patternCheck(wedgeOrder, subPattern, operand1, operand2, "1")
  }
  if (type == "complex") {
    resPolyMat <- subCoeffMats[[wedgeOrder]]
  } else {
    resPolyMat <- subCoeffMats[[wedgeOrder]]
  }

  return(resPolyMat)
}

tipLabelCoeffMat <- function(tree) {

  # construct the child matrix of the tree
  nnodes <- tree$Nnode # number of internal nodes
  nTips <- length(tree$tip.label)
  nodeids <- (nTips + 1):(nTips + nnodes) # label starting after last tip to the total nodes
  ChildMatrix <- t(sapply(nodeids, function(x) tree$edge[tree$edge[, 1] == x, 2]))
  ChildMatrix <- cbind(nodeids, ChildMatrix) #  first col are the internal nodes

  # determine the order of wedges needed to generate the trees coefficient matrix
  numTips <- ChildMatrix[1, 1] - 1

  tipLabs <- tree$tip.label

  wedgeOrder <- unique(rev(as.vector(t(ChildMatrix))))
  wedgeOrder <- tipLabs[wedgeOrder]
  wedgeOrder[is.na(wedgeOrder)] <- 1

  done <- F

  subCoeffMats <- list()
  # intialize first two entries for leaves and cherries
  tipLabs <- unique(tipLabs)
  if (length(tipLabs) == 2) {
    subCoeffMats[[as.character(tipLabs[[1]])]] <- matrix(data = c(0, 1 + 0i), nrow = 1, ncol = 2)
    subCoeffMats[[as.character(tipLabs[[2]])]] <- matrix(data = c(0, 1 + 0i), nrow = 2, ncol = 1)
  } else if (length(tipLabs) == 1) {
    subCoeffMats[[as.character(tipLabs[[1]])]] <- matrix(data = c(0, 1 + 0i), nrow = 1, ncol = 2)
    subCoeffMats[[as.character(tipLabs[[1]])]] <- matrix(data = c(0, 1 + 0i), nrow = 2, ncol = 1)
  } else {
    stop("only phylo objects with two tip labels are currently supported")
  }

  # continue while there is still stuff left to wedge
  while (!(length(wedgeOrder) == 1)) {

    # find the next two operands to wedge
    j <- 3
    while (TRUE) {
      if (wedgeOrder[j] == 1) {
        operand1 <- wedgeOrder[j - 2]
        operand2 <- wedgeOrder[j - 1]
        break
      }
      j <- j + 1
    }

    # calculate the pattern
    subPattern <- paste(operand1, operand2, "1", sep = "")
    old <- paste(wedgeOrder, collapse = " ", sep = ",")
    oldTest <- paste(wedgeOrder, collapse = "")

    # calculate wedge
    subCoeffMats[[as.character(subPattern)]] <- tipLabelWedge(subCoeffMats[[operand1]], subCoeffMats[[operand2]])

    # replace subsequent operations of the same wedge
    wedgeOrder <- patternCheck(wedgeOrder, subPattern, operand1, operand2, "1")
  }

  return(subCoeffMats[[wedgeOrder]])
}

# function replaces all 3 object patterns (including patterns with the first 2 objects order swapped) in a vector of objects
patternCheck <- function(completePattern, replacement, firstNum, secondNum, thirdNum) {

  # check the pattern for the sequence with the option of having the first and second number swapped
  testPoints <- seq(1, length(completePattern) - 2, 1)
  matchingLocationsForward <- completePattern[testPoints] == firstNum & completePattern[(testPoints + 1)] == secondNum & completePattern[(testPoints + 2)] == thirdNum
  matchingLocationsBackward <- completePattern[testPoints] == secondNum & completePattern[(testPoints + 1)] == firstNum & completePattern[(testPoints + 2)] == thirdNum
  matchingLocations <- matchingLocationsForward | matchingLocationsBackward
  pointsToSimplify <- testPoints[matchingLocations]

  # assign the replacement to the first position and remove the other two positions
  completePattern[pointsToSimplify] <- replacement
  pointsToSimplify <- c((pointsToSimplify + 1), (pointsToSimplify + 2))
  completePattern <- completePattern[-pointsToSimplify]

  return(completePattern)
}

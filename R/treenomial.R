#' Converts trees (objects of class "phylo") to coefficient matrices
#'
#' @param tree an object or list of class "phylo"
#' @return The resulting coefficient matrix or matrices.
#' @import Matrix
#' @import pbapply
#' @import parallel
#' @import ape
#' @examples
#' # example goes here
#' @export
coefficientMatrix <- function(tree, complexMode = FALSE) {
  # determine what the incoming trees look like
  tryCatch({
    tree <<- as.phylo(tree)
    singleTree <<- TRUE
  }, error = function(e) {
    singleTree <<- FALSE
    tree <<- lapply(tree, as.phylo)
  })

  # call appropiately for the number and size of the input trees
  if (singleTree) {
    # print("singleTree")
    coefficientMat <- singleCoeffMat(tree, loadingOn = TRUE)
  } else {
    # need to test more for best time to go multicore
    if (length(tree) >= 300 && tree[[1]]$Nnode > 400) {
      cl <- makeCluster(detectCores())
      # print("multicore")
      coefficientMat <- vector("list", length = length(tree))
      coefficientMat <- pblapply(tree, singleCoeffMat, loadingOn = FALSE, cl = cl)

      stopCluster(cl)
    } else {
      # print("singleMulti")
      coefficientMat <- vector("list", length = length(tree))
      coefficientMat <- pblapply(tree, singleCoeffMat, loadingOn = FALSE)
    }
  }

  return(coefficientMat)
}

singleCoeffMat <- function(tree, loadingOn) {

  # construct the child matrix of the tree
  nnodes <- tree$Nnode # number of internal nodes
  nTips <- length(tree$tip.label)
  nodeids <- (nTips + 1):(nTips + nnodes) # label starting after last tip to the total nodes
  ChildMatrix <- t(sapply(nodeids, function(x) tree$edge[tree$edge[, 1] == x, 2]))
  ChildMatrix <- cbind(nodeids, ChildMatrix) #  first col are the internal nodes

  # determine the order of wedges needed to generate the trees polynomial matrix
  numTips <- ChildMatrix[1, 1] - 1
  wedgeOrder <- ifelse(unique(rev(as.vector(t(ChildMatrix)))) <= numTips, "0", "1")

  if (loadingOn) {
    ## loading bar stuff
    maxLength <- length(wedgeOrder)
    pb <- txtProgressBar(min = 0, max = maxLength - 1, initial = 0, style = 3)
  }

  done <- F
  hash <- 4

  # intialize first two entries for leaves and cherries
  subPolys <- list()


  subPolys[["0"]] <- sparseMatrix(1, 2, x = 1) # leaf
  subPolys[["001"]] <- sparseMatrix(c(1, 2), c(3, 1), x = c(1, 1)) # cherry


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
    hash <- paste(operand1, operand2, "1", sep = "")
    old <- paste(wedgeOrder, collapse = " ", sep = ",")
    oldTest <- paste(wedgeOrder, collapse = "")

    # calculate wedge
    subPolys[[as.character(hash)]] <- wedge(subPolys[[operand1]], subPolys[[operand2]])


    # replace subsequent operations of the same wedge
    wedgeOrder <- patternCheck(wedgeOrder, hash, operand1, operand2, "1")

    if (loadingOn) {
      ## loading bar stuff
      setTxtProgressBar(pb, maxLength - length(wedgeOrder))
    }
  }
  resPolyMat <- as.matrix(subPolys[[wedgeOrder]])

  return(resPolyMat)
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

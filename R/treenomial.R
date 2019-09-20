#' Converts trees (objects of class "phylo") to coefficient matrices
#'
#' @param tree an object or list of class "phylo"
#' @param complex whether or not to return the polynomial with y = i
#' @param tipLabel
#' @param normalize whether or not to use a row sum normilization of the coefficients
#' @param progressBar add a progress bar to track performance
#' @param cl
#' @return The resulting coefficient matrix or matrices.
#' @import Matrix
#' @import pbapply
#' @importFrom ape as.phylo
#' @examples
#' library(ape)
#' # generate a tree
#' tree <- rtree(n = 500, rooted = TRUE)
#' # get the coefficient matrix describing the polynomial for this tree
#' coefficientMatrix(tree)
#' @export
coefficientMatrix <- function(tree, complex = FALSE, tipLabel = FALSE, normalize = FALSE, progressBar = FALSE, cl = NULL) {

  # determine what the incoming trees look like
  tryCatch({
    tree <<- as.phylo(tree)
    singleTree <<- TRUE
  }, error = function(e) {
    singleTree <<- FALSE
    tree <<- lapply(tree, as.phylo)
  })

  # set the max times the progress bar is updated and check if user wants a progress bar
  pboptions(nout = 50, type = ifelse(progressBar, "timer", "none"))

  if(!tipLabel){
    # call appropiately for the number and size of the input trees
    if (singleTree) {
      coefficientMat <- singleCoeffMat(tree,complex = complex)
    } else {
      coefficientMat <- vector("list", length = length(tree))
      coefficientMat <- pblapply(tree, singleCoeffMat,complex = complex, cl = cl)
    }

    # normalize if selected
    if(normalize && !complex){
      coefficientMat <- lapply(coefficientMat, function(x){
        rowSum <- rowSums(x)
        x[rowSum != 0,] <- x[rowSum != 0,]/rowSum[rowSum != 0]
        return(x)
      })
    }

  } else {
    # call appropiately for the number and size of the input trees
    if (singleTree) {
      coefficientMat <- tipLabelCoeffMat(tree)
    } else {
      coefficientMat <- vector("list", length = length(tree))
      coefficientMat <- pblapply(tree, tipLabelCoeffMat, cl = cl)

      # tip labels coeff mats will be of different size depending on tip labels, this adds zeroes to align
      maxRowCol <- sapply(coefficientMat, function(x) c(nrow(x),ncol(x)))
      maxRowCol <- c(max(maxRowCol[1,]),max(maxRowCol[2,]))
      coefficientMat <- lapply(coefficientMat, function(x){
        amountSmaller  <- maxRowCol - dim(x)
        x <- cbind(x,matrix(data = 0, ncol = amountSmaller[2], nrow = nrow(x)))
        x <- rbind(x,matrix(data = 0, ncol = ncol(x), nrow = amountSmaller[1]))
      })
    }
  }

  return(coefficientMat)
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

  wedgeOrder <- unique(rev(as.vector(t(ChildMatrix))))
  wedgeOrder <- sapply(wedgeOrder, function(i) {
    unname(ifelse(i <= numTips, tree$tip.label[i], "1"))
  })

  done <- F

  # intialize first two entries for leaves and cherries
  subCoeffMats <- list()
  subCoeffMats[["t1"]] <- matrix(data = c(0,1 + 0i), nrow = 1, ncol = 2)
  subCoeffMats[["t2"]] <- matrix(data = c(0,1 + 0i),  nrow = 2, ncol = 1)

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

  return( subCoeffMats[[wedgeOrder]])
}

singleCoeffMat <- function(tree,complex = FALSE) {

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

  if(complex){
    subCoeffMats[["0"]] <- c(0,1) # leaf
    subCoeffMats[["001"]] <- c(1i,0,1) # cherry
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
    if(complex){
      subCoeffMats[[as.character(subPattern)]] <- wedgeComplex(subCoeffMats[[operand1]], subCoeffMats[[operand2]])
    } else {
      subCoeffMats[[as.character(subPattern)]] <- wedge(subCoeffMats[[operand1]], subCoeffMats[[operand2]])
    }

    # replace subsequent operations of the same wedge
    wedgeOrder <- patternCheck(wedgeOrder, subPattern, operand1, operand2, "1")

  }
  if(complex){
    resPolyMat <- subCoeffMats[[wedgeOrder]]
  } else {
    resPolyMat <- as.matrix(subCoeffMats[[wedgeOrder]])
  }

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

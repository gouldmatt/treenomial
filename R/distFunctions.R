#' Calculates the distance matrix from multiple coefficient matrices
#'
#' @param coefficientMatrices list of complex or real coefficient matrices
#' @param method method to use when calculating coefficient distances: "logL1", "wLogL1", "b"  or "c"
#' @param progressBar add a progress bar to track progress
#' @param smallDistanceMatrix whether to return a distance matrix (TRUE) or a distance value (FALSE) if only two coefficient matrices are given
#' @return distance matrix calculated from argument coefficient matrices
#' @note the complex coefficient vector only supports the "logL1" method
#' @note by default a list of two coefficient matrices will return a distance value rather than a 2*2 distance matrix
#' @import Matrix
#' @importFrom utils combn
#' @import pbapply
#' @examples
#'
#'
#'
#'
#' @export
coeffDist <- function(coefficientMatrices, method = "logL1", progressBar = FALSE, smallDistanceMatrix = FALSE) {

   # check input arguments
  if(class(coefficientMatrices) == "list"){
    if(class(coefficientMatrices[[1]]) == "dgCMatrix"){
      coefficientMatrices <- lapply(coefficientMatrices, function(i){as.matrix(i)})
      distMat <- coeffDistRcpp(coefficientMatrices, method = method, progressBar = progressBar)

    } else if(class(coefficientMatrices[[1]]) == "matrix"){

      distMat <- coeffDistRcpp(coefficientMatrices, method = method, progressBar = progressBar)

    } else if(class(coefficientMatrices[[1]]) == "complex"){
      distMat <- coeffDistRcpp(coefficientMatrices, method = "logL1Complex", progressBar = progressBar)

    } else {
      stop("input must be a list of coefficient matrices")
    }

  } else {
    stop("input must be a list of coefficient matrices")
  }

  # align
  coefficientMatrices <- coefficientAlign(coefficientMatrices)


  # distMat <- new("dspMatrix", Dim = as.integer(c(numCoeffs,numCoeffs)), x = u, uplo = "U")
  rownames(distMat) <- names(coefficientMatrices)
  colnames(distMat) <- names(coefficientMatrices)

  if(!smallDistanceMatrix && length(coefficientMatrices) == 2){
    return(distMat[1,2])
  }

  return(distMat)
}

#' Calculates the distance matrix from a list of phylo objects
#' @param type type of the polynomial one of:
#' "real" to use real polynomials \cr
#' "complex" to use the real polynomial with y = 1 + i \cr
#' "tipLabel" to use polynomial that utilize binary trait tip labels on the phylo objects \cr
#' @inheritParams coeffDist
#' @inheritParams coeffMatrix
#' @return a distance matrix
#' @examples
#' #' require(ape)
#  # distance matrix for 10 trees of 30 tips
#' phyloDist(trees,method = "logL1", type = "tipLabel")
#'
#' @export
phyloDist <- function(trees, method = "logL1", type = "real", progressBar = FALSE, smallDistanceMatrix = FALSE, cl = NULL){

  if(type == "tipLabel"){
    coeffDist(coefficientMatrices = coeffMatrix(trees, type = "tipLabel" ,progressBar,cl), method = type, progressBar = progressBar, smallDistanceMatrix = smallDistanceMatrix)
  } else {
    coeffDist(coeffMatrix(trees, type = type ,progressBar,cl), method = method, progressBar = progressBar, smallDistanceMatrix = smallDistanceMatrix)
  }
}

#' Plot the min/max distance tree from a target tree in a distance matrix
#'
#' @param trees list of phylo objects
#' @param distMatrix
#' @param target
#' @param comparison find the "min" or the "max" distance tree in the distance matrix
#' @param plotFacing whether to plot the trees with the tips facing each other
#' @param returnNearestInfo whether to return the smallest/max distance phylo object and its distance to target (TRUE) or have the function return no value (FALSE)
#' @return
#' @examples
#'
#' @export
phyloDistPlot <- function(trees, distMatrix, target, comparison = "min", plotFacing = FALSE, returnNearestInfo = FALSE){

  facing = ifelse(plotFacing,"leftwards","rightwards")

  if(class(target) == "character"){
    targetName = target
    target = which(rownames(distMatrix) == target)
  } else {
    targetName = target
  }

  if(comparison == "min"){
    minTreeIndex <- which.min(distMatrix[target,-target])
    if(minTreeIndex > target){
       minTreeIndex <- minTreeIndex + 1
    }

    par(mfrow=c(1,2), oma=c(2,0,2,0))
    plot.phylo(ladderize(as.phylo(trees[[target]])), show.tip.label = FALSE, main = paste("Target:",targetName), use.edge.length = F)
    plot.phylo(ladderize(as.phylo(trees[[minTreeIndex]])), show.tip.label = FALSE, main = paste("Min. Dist. Tree:",names(minTreeIndex)), use.edge.length = F, direction = facing)

    mtext(side=1, paste("Distance Value:",as.character(distMatrix[target,minTreeIndex])), outer=TRUE)

    if(returnNearestInfo){
      list(minTree = trees[[minTreeIndex]], distance = distMatrix[target,minTreeIndex])
    }

  } else if(comparison == "max"){

    maxTreeIndex <- which.max(distMatrix[target,])

    par(mfrow=c(1,2), oma=c(2,0,2,0))
    plot.phylo(ladderize(as.phylo(trees[[target]])), show.tip.label = FALSE, main = paste("Target:",targetName), use.edge.length = F)
    plot.phylo(ladderize(as.phylo(trees[[maxTreeIndex]])), show.tip.label = FALSE, main = paste("Max. Dist. Tree:",names(maxTreeIndex)), use.edge.length = F, direction = facing)

    mtext(side=1, paste("Distance Value:",as.character(distMatrix[target,maxTreeIndex])), outer=TRUE)

    if(returnNearestInfo){
      list(maxTree = trees[[maxTreeIndex]], distance = distMatrix[target,maxTreeIndex])
    }

  }

}

coefficientAlign <- function(coeffMats){

  # check for complex poly
  if(is.complex(coeffMats[[1]])){

    # check for the matrix of tip labelled complex matrices
    if(class(coeffMats[[1]]) == "matrix"){
      # tip labels coeff mats will be of different size depending on tip labels, this adds zeroes to align
      maxRowCol <- sapply(coeffMats, function(x) c(nrow(x),ncol(x)))
      maxRowCol <- c(max(maxRowCol[1,]),max(maxRowCol[2,]))
      coeffMats <- lapply(coeffMats, function(x){
        amountSmaller  <- maxRowCol - dim(x)
        x <- cbind(x,matrix(data = 0, ncol = amountSmaller[2], nrow = nrow(x)))
        x <- rbind(x,matrix(data = 0, ncol = ncol(x), nrow = amountSmaller[1]))
      })
      return(coeffMats)
    }

    coeffsLengths <- lengths(coeffMats)

    # determine the size of the largest matrix in the list
    maxSize <- max(coeffsLengths)

    # check there exists at least one smaller matrix
    if(!any(coeffsLengths < maxSize)){
      return(coeffMats)
    }

    # go through and align smaller complex vectors
    coeffMats <- lapply(coeffMats, function(x){
      # check if alignment is neccessary
      if(length(x) < maxSize){
        # align
        alignedVector <- vector(mode = "complex", length = maxSize)
        alignedVector[(maxSize - length(x) + 1):maxSize] <- x
        return(alignedVector)
      } else {
        return(x)
      }

    })
    return(coeffMats)
  }

  colSizes <- vapply(coeffMats, ncol, FUN.VALUE =  numeric(1))

  # determine the size of the largest matrix in the list
  maxSize <- max(colSizes)
  maxRows <- ceiling(maxSize/2)

  # check there exists at least one smaller matrix
  if(!any(colSizes < maxSize)){
    return(coeffMats)
  }

  # go through and align smaller matrices
  coeffMats <- lapply(coeffMats, function(x){
    # check if alignment is neccessary
    if(ncol(x) < maxSize){
      # align
      alignedMatrix <- matrix(data = 0, nrow = maxRows, ncol = maxSize)
      alignedMatrix[1:nrow(x),(maxSize - ncol(x) + 1):maxSize] <- x
      return(alignedMatrix)
    } else {
      return(x)
    }

  })
  return(coeffMats)
}



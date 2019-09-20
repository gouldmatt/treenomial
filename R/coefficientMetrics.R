#' Aligns a list of coefficent matrices to the max size in the argument list
#'
#'
#' @param coeffMats list of real or complex coefficient matrices of different sizes to be aligned
#' @return the aligned list of coefficient matrices
#' @examples
#' # example goes here
#' @export
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


#' Calculates the distance matrix of multiple coefficient matrices
#'
#'
#' @param coefficientMatrices list of coefficient matrices
#' @param method method to use when calculating coefficient distances: "sumLogDiff", "sumLogDiffComplex","sumLogDiffLabels" or "fager"
#' @param progressBar add a progress bar to track progress
#' @param cl a cluster object created by makeCluster
#' @param smallDistanceMatrix whether to return a distance matrix (TRUE) or a distance value (FALSE) if only two coefficient matrices are given
#' @return distance matrix calculated from argument coefficient matrices
#' @note by default a list of two coefficient matrices will return a distance value rather than a 2*2 distance matrix
#' @import Matrix
#' @importFrom utils combn
#' @import pbapply
#' @examples
#' # example goes here
#' @export
coefficientDist <- function(coefficientMatrices, method = "sumLogDiff", progressBar = FALSE, smallDistanceMatrix = FALSE, cl = NULL) {

  # determine and define the places to calculate in the distance matrix
  upperDist <- t(combn(length(coefficientMatrices), 2))
  distMat <- matrix(data = 0, nrow = length(coefficientMatrices), ncol = length(coefficientMatrices))

  # set the max times the progress bar is updated and check if user wants a progress bar
  pboptions(nout = ceiling(length(coefficientMatrices)/4), type = ifelse(progressBar, "timer", "none"))

  if (method == "sumLogDiff") {
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      sumLogDiff(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    }, cl = cl)
  } else if (method == "sumLogDiffComplex") {
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      sumLogDiffComplex(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    }, cl = cl)

  } else if (method == "sumLogDiffLabels"){
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      sumLogDiffLabels(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    }, cl = cl)

  } else if (method == "fager") {
    coefficientMatrices <- lapply(coefficientMatrices, as.logical)
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      coeffMatA <- coefficientMatrices[[i[1]]]
      coeffMatB <- coefficientMatrices[[i[2]]]

      c <- as.numeric(length(coeffMatA[coeffMatA == TRUE & coeffMatB == FALSE]))
      b <- as.numeric(length(coeffMatA[coeffMatA == FALSE & coeffMatB == TRUE]))
      d <- as.numeric(length(coeffMatA[coeffMatA == FALSE & coeffMatB == FALSE]))
      a <- as.numeric(length(coeffMatA[coeffMatA == TRUE & coeffMatB == TRUE]))

      return((a/sqrt((a+b)*(a+c))) - sqrt(a+c)/2)
      #sqrt(b+c)/(a+b+c+d)**2
      # return(sqrt(b+c)/(a+b+c+d)^2)
      #return(a+b+c+d)
    }, cl = cl)
  } else {
    # error
  }

  if(smallDistanceMatrix == FALSE && length(coefficientMatrices) == 2){
    return(distMat[1,2])
  }

  distMat <- t(distMat) + distMat

  return(distMat)
}


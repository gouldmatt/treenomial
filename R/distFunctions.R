#' Calculates the distance matrix from multiple coefficient matrices
#'
#' @param coefficientMatrices list of complex or real coefficient matrices
#' @param method method to use when calculating coefficient distances: "logL1", "wLogL1", "b"  or "c"
#' @param progressBar add a progress bar to track progress
#' @param smallDistanceMatrix whether to return a distance matrix (TRUE) or a distance value (FALSE) if only two coefficient matrices are given
#' @return distance matrix calculated from argument coefficient matrices
#' @note the complex coefficient vector only support the "logL1" method
#' @note by default a list of two coefficient matrices will return a distance value rather than a 2*2 distance matrix
#' @import Matrix
#' @importFrom utils combn
#' @import pbapply
#' @examples
#' # example goes here
#' @export
coefficientDist <- function(coefficientMatrices, method = "logL1", progressBar = FALSE, smallDistanceMatrix = FALSE) {
  ############ change to using symmetric matrix from Matrix package for distance matrix

  numCoeffs <- as.numeric(length(coefficientMatrices))

  u <-  rep(0, (numCoeffs*numCoeffs-numCoeffs)/2 + numCoeffs)
  u[(1:numCoeffs)*((1:numCoeffs)+1)/2] <- 0

  distMat <- new("dspMatrix", Dim = as.integer(c(numCoeffs,numCoeffs)), x = u, uplo = "U")
  rownames(distMat) <- names(coefficientMatrices)
  colnames(distMat) <- names(coefficientMatrices)

  # determine and define the places to calculate in the distance matrix
  # upperDist <- t(combn(length(coefficientMatrices), 2))
  # distMat <- symm(data = 0, nrow = length(coefficientMatrices), ncol = length(coefficientMatrices), sparse = FALSE, )

  # vapply((1:4), function(i){ c(rep(i,5-i),(i+1):5)}, FUN.VALUE = list())


  i <- Vectorize(FUN = function(i,numCoeffs){ rep(i,numCoeffs-i)})
  j <- Vectorize(FUN = function(i,numCoeffs){ (i+1):numCoeffs})
  i <- unlist(i(1:(numCoeffs-1),numCoeffs = numCoeffs))
  j <- unlist(j(1:(numCoeffs-1),numCoeffs = numCoeffs))

  upperDist <- cbind(i,j)

  # set the max times the progress bar is updated and check if user wants a progress bar
  pboptions(nout = ceiling(length(coefficientMatrices)/4), type = ifelse(progressBar, "txt", "none"))

  # check input arguments
  if(class(coefficientMatrices) == "list"){
    if(class(coefficientMatrices[[1]]) == "dgCMatrix"){
      complex <- FALSE

    } else if(class(coefficientMatrices[[1]]) == "complex"){

      complex <- TRUE
    } else {
      stop("input must be a list of coefficient matrices")
    }

  } else {
    stop("input must be a list of coefficient matrices")
  }


  if (complex) {
    if(method != "logL1"){
      warning("only the logL1 method is available for the complex coefficients, using this method")
    }

    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      logL1Complex(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    })

  } else if (method == "logL1") {
    coefficientMatrices <- lapply(coefficientMatrices, function(i){as.matrix(i)})
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      logL1(coefficientMatrices[[i[1]]],coefficientMatrices[[i[2]]])
    })


  } else if (method == "c") {
    coefficientMatrices <- lapply(coefficientMatrices, function(i){as.logical(i)})
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {

      coeffMatA <- coefficientMatrices[[i[1]]]
      coeffMatB <- coefficientMatrices[[i[2]]]

      # return the c disimliarity
      as.numeric(length(coeffMatA[coeffMatA == TRUE & coeffMatB == FALSE]))
    })

  } else if (method == "b") {
    coefficientMatrices <- lapply(coefficientMatrices, function(i){ as.vector( as.logical(i))})
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {

      coeffMatA <- coefficientMatrices[[i[1]]]
      coeffMatB <- coefficientMatrices[[i[2]]]

      # return the b disiliarity
      as.numeric(length(coeffMatA[coeffMatA == FALSE & coeffMatB == TRUE]))

    })

  } else if(method == "wLogL1"){
    coefficientMatrices <- lapply(coefficientMatrices, function(i){as.matrix(i)})
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {


      logDiffMat <- rowSums(log(1 + abs(coefficientMatrices[[i[1]]]-coefficientMatrices[[i[2]]])))

      weightVect <- c(1,(1:(nrow(coefficientMatrices[[1]])-1))^(-2))

      sum(logDiffMat*weightVect)

    })

  } else {
    stop("missing or incorrect method, see documentation for available methods")
  }

  if(smallDistanceMatrix == FALSE && length(coefficientMatrices) == 2){
    return(distMat[1,2])
  }


  # symmetric
  distMat <- t(distMat) + distMat

  return(distMat)
}

#' Calculates the distance matrix from a list of trees
#'
#' @param trees list of phylo objects
#' @param method method to use when calculating coefficient distances: "logL1", "wLogL1", "b"  or "c"
#' @param progressBar add a progress bar to track progress
#' @param cl a cluster object created by makeCluster
#' @param smallDistanceMatrix whether to return a distance matrix (TRUE) or a distance value (FALSE) if only two coefficient matrices are given
#' @return distance matrix calculated from argument coefficient matrices
#' @note the complex coefficient vector only support the "logL1" method
#' @note by default a list of two coefficient matrices will return a distance value rather than a 2*2 distance matrix
#' @examples
#' @export
phyloDist <- function(trees, method = "logL1", type = "real", progressBar = FALSE, smallDistanceMatrix = FALSE, cl = NULL){

  coefficientDist(coefficientMatrix(trees, type = type ,progressBar,cl), method = method, progressBar = progressBar, smallDistanceMatrix = smallDistanceMatrix)

}

#' @export
distPlot <- function(trees, distanceMatrix, target, comparison = "min"){



  if(comparison == "min"){
    minTreeIndex <- which.min(distanceMatrix[target,(distanceMatrix[target,] != 0)]) + 1
    par(mfrow=c(1,2))

    plot.phylo(ladderize(as.phylo(trees[[target]])), show.tip.label = FALSE, main = target, use.edge.length = F)
    plot.phylo(ladderize(as.phylo(trees[[minTreeIndex]])), show.tip.label = FALSE, main = names(minTreeIndex), use.edge.length = F)
    mtext(paste("Distance Value:",as.character(distanceMatrix[target,minTreeIndex])), outer = FALSE, cex = NA, side = 1)

  } else if(comparison == "max"){

    maxTreeIndex <- which.max(distanceMatrix[target,(distanceMatrix[target,] != 0)]) + 1
    par(mfrow=c(1,2))

    plot.phylo(ladderize(as.phylo(trees[[target]])), show.tip.label = FALSE, main = target, use.edge.length = F)
    plot.phylo(ladderize(as.phylo(trees[[maxTreeIndex]])), show.tip.label = FALSE, main = names(maxTreeIndex), use.edge.length = F)
    mtext(paste("Distance Value:",as.character(distanceMatrix[target,maxTreeIndex])), outer = FALSE, cex = NA, side = 1)
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


#' Calculates the distance between coefficient matrices
#'
#' Calculates the distance between two coefficient matrices or a coefficient matrix and a list of coefficient matrices.
#'
#' @param x single coefficient matrix to find distances to
#' @param Y list or single coefficient matrix
#' @inheritParams polyToDistMat
#' @return vector of distances
#' @note the complex coefficient vector and the complex tip label coefficient matrix only support the \dQuote{logDiff} method
#' @examples
#'
#'
#' @export
polyDist <- function(x, Y, method = "logDiff"){
  # check input arguments

      if(typeof(x) == "list"){
        x <- x[[1]]
      }

      if(typeof(x[[1]]) == "double"){
        coefficientMatrices <- alignCoeffs(c(list(x),Y), type = "real")
        compareCoeffRcpp( coefficientMatrices, method = method)

      } else if(typeof(x[[1]]) == "complex"){
        if(dim(x[[1]])[[1]] == 1){
          if(method != "logDiff") warning("only the logDiff method is available for the complex polynomial")
          coefficientMatrices <- alignCoeffs(c(x,Y), type = "complex")
          compareCoeffRcpp( coefficientMatrices, method = "logDiffComplex")

        } else {
          if(method != "logDiff") warning("only the logDiff method is available for binary trait label polynomial")
          coefficientMatrices <- alignCoeffs(c(x,Y), type = "binTipLabel")
          compareCoeffRcpp( coefficientMatrices, method = "tipLab")
        }
      } else {
        stop("invalid input")
      }

}


#' Calculates the distance between coefficient matrices
#'
#' Calculates the distance between two coefficient matrices or a coefficient matrix and a list of coefficient matrices.
#'
#' @param x single coefficient matrix to find distances to
#' @param Y list or single coefficient matrix
#' @inheritParams treeToDistMat
#' @return vector of distances
#' @note the complex coefficient vector and the complex tip label coefficient matrix only supports the \dQuote{logDiff} method
#' @examples
#'
#'
#' @export
treeDist <- function(x, Y, type = "real", method = "logDiff"){
  coeffs <- treeToPoly(c(x,Y), type = type)
  coeffs <- alignPoly(coeffs)
  compareCoeffRcpp(coeffs, method = method)
}

#' Calculates the distance matrix from multiple coefficient matrices
#'
#'
#' @param coefficientMatrices list of complex or real coefficient matrices
#' @param method method to use when calculating coefficient distances:
#' \describe{
#'   \item{\dQuote{logDiff}}{for two coefficient matrices A and B returns sum(log(1+abs(A-B))}
#'   \item{\dQuote{wLogDiff}}{performs the \dQuote{logDiff} method with weights on the rows}
#'   \item{\dQuote{pa}}{total pairs where the coefficient is present in one matrix and absent in the other (presence-absence)}
#'   \item{\dQuote{ap}}{opposite comparison of pa (absence-presence)}
#' }
#' @return distance matrix calculated from argument coefficient matrices
#' @note \itemize{
#'   \item the complex coefficient vector and the complex tip label coefficient matrix only support the \dQuote{logDiff} method
#'   \item \dQuote{pa} and \dQuote{ap} force symmetry in the output distance matrix
#' }
#' @examples
#'
#'
#'
#'
#' @export
polyToDistMat <- function(coefficientMatrices, method = "logDiff") {

  # check input arguments
  if(class(coefficientMatrices) == "list"){
    if(class(coefficientMatrices[[1]]) == "matrix"){
      if(typeof(coefficientMatrices[[1]]) == "double"){
        coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "real")
        distMat <- coeffDistRcpp(coefficientMatrices, method = method)
      } else if(typeof(coefficientMatrices[[1]]) == "complex"){
        if(dim(coefficientMatrices[[1]])[[1]] == 1){
          if(method != "logDiff") warning("only the logDiff method is available for the complex polynomial")
          coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "complex")
          distMat <- coeffDistRcpp(coefficientMatrices, method = "logDiffComplex")
        } else {
          if(method != "logDiff") warning("only the logDiff method is available for binary trait label polynomial")
          coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "binTipLabel")
          distMat <- coeffDistRcpp(coefficientMatrices, method = "tipLab")
        }
      } else {
        stop("invalid input")
      }

    } else {
      stop("input must be a list of coefficient matrices")
    }

  } else {
    stop("input must be a list of coefficient matrices")
  }

  # distMat <- new("dspMatrix", Dim = as.integer(c(numCoeffs,numCoeffs)), x = u, uplo = "U")
  rownames(distMat) <- names(coefficientMatrices)
  colnames(distMat) <- names(coefficientMatrices)

  return(distMat)
}

#' Calculates the distance matrix from a list of phylo objects
#' @inheritParams polyToDistMat
#' @inheritParams treeToPoly
#' @return a distance matrix
#' @export
#' @examples
#' library(treenomial)
#' library(ape)
#' # distance matrix for 10 trees of 30 tips
#' phyloDist(rmtree(10,30),method = "wLogDiff")
#'
#' @export
treeToDistMat <- function(trees, method = "logDiff", type = "real"){
  polyToDistMat(treeToPoly(trees, type = type), method = method)
}

#' Plot the min/max distance trees from a target tree
#'
#' @param target the phylo object of the tree to calculate the distances to
#' @param trees a list of phylo objects to compare with the \strong{target}
#' @param n the number of trees to find and plot
#' @param comparison whether to find the \dQuote{min} or the \dQuote{max} distance trees from the \strong{target}
#' @return a list of lists containing the \strong{n} min/max distance trees and their distances to \strong{target}
#' @inheritParams treeToDistMat
#' @examples
#'
# trees <- c(rmtree(1000,50),rmtree(10,9))
# target <- rtree(50)
# minTrees <- plotExtremeTrees(target,trees,2, comparison = "min")
#'
#'
#' @export
plotExtremeTrees <- function(target, trees, n, method = "logDiff", type = "real", comparison = "min"){

  distances <- treeDist(target, trees, type = type, method = method)

  # if greater then 16 trees break up over multiple pages
  if(n < 15){
    par(mfrow=c(ceiling(n/4),ifelse(ceiling(n/4) == 1,1+n,4)), oma=c(2,0,2,0))
  } else {
    par(mfrow=c(4,4), oma=c(2,0,2,0))
  }

  target$edge.length <- rep(1,length(target$edge.length))

  plot.phylo(ladderize(target), show.tip.label = FALSE, main = "Target:", use.edge.length = T)

  if(comparison == "min"){
    orderMin <- order(distances)

    minList <- vector("list",n)

    for(i in 1:n){

      minTree <- trees[[orderMin[i]]]
      minTree$edge.length <- rep(1,length(minTree$edge.length))

      currName <- names(trees)[orderMin[i]]
      minTitle <- ifelse(!is.null(currName),paste(currName, "tree distance: ",signif(distances[orderMin[i]],6)),paste("Distance: ",signif(distances[orderMin[i]],6)))

      plot.phylo(ladderize(minTree), show.tip.label = FALSE, main = minTitle, use.edge.length = T)

      minList[[i]] <- list(tree = minTree, distance = distances[[orderMin[i]]])
    }

    return(minList)

  } else if(comparison == "max"){
    orderMax <- order(distances, decreasing = TRUE)

    maxList <- vector("list",n)

    for(i in 1:n){

      maxTree <- trees[[orderMax[i]]]
      maxTree$edge.length <- rep(1,length(maxTree$edge.length))

      currName <- names(trees)[orderMax[i]]
      minTitle <- ifelse(!is.null(currName),paste(currName, "tree distance: ",signif(distances[orderMax[i]],6)),paste("Distance: ",signif(distances[orderMax[i]],6)))

      plot.phylo(ladderize(maxTree), show.tip.label = FALSE, main = minTitle, use.edge.length = T)

      maxList[[i]] <- list(tree = maxTree, distance = distances[[orderMax[i]]])
    }

    return(maxList)

  } else {
    stop("invalid comparison")
  }

}

#' Calculates the distance matrix of multiple coefficient matrices
#'
#' Uses a specified method to construct a distance matrix from coefficient matrices.
#' For a list of methods see \link[parallelDist]{parDist}.
#' @param coefficientMatrices list of coefficient matrices
#' @param method method to use when calculating coefficient distances
#' @return distance matrix calculated from argument coefficient matrices
#' @import Matrix
#' @import parallelDist
#' @import RcppArmadillo
#' @import RcppXPtrUtils
#' @examples
#' # example goes here
#' @export
coefficientDist <- function(coefficientMatrices, method = "canberra") {
  # select and perform method
  if (method == "sumLogdiff") {
    euclideanFuncPtr <- cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) {
                                arma::mat temp = arma::log1p(arma::abs(A-B));
                                return arma::accu(temp); }",
      depends = c("RcppArmadillo")
    )

    disMat <- as.matrix(parDist(coefficientMatrices, method = "custom", func = euclideanFuncPtr, upper = TRUE))
  } else if (method == "geodesic" || method == "chord") {
    # flatten the coefficient matrices
    coefficientMatrices <- lapply(coefficientMatrices, function(i) {
      t(matrix(as.vector(i)))
    })
    disMat <- as.matrix(parDist(coefficientMatrices, method = method, upper = TRUE))
  } else {
    disMat <- as.matrix(parDist(coefficientMatrices, method = method, upper = TRUE))
  }

  return(disMat)
}

#' Builds a 2d or 3d mds plot from tree or coefficient matrix data
#'
#' Either calculates or recieves coefficient matrix data of phylo objects and uses a distance method
#' on the coefficients to construct an mds plot.
#' @param trees data frame of phylo objects used to calculate coefficient matrices
#' @param coefficientMatrices list of coefficient matrices if not passing trees data
#' @param method method to use when caculating coefficient distances
#' @param dim number of dims of resulting mds plot
#' @return 2d or 3d MDS plot
#' @import Matrix
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom plotly %>%
#' @import ggplot2
#' @examples
#' library(apTreeshape)
#' # generate some forests
#' numTrees <- 10
#' numTips <- 500
#'
#' pdaTrees <- rtreeshape(numTrees,tip.number=numTips,model="pda")
#' yuleTrees <- rtreeshape(numTrees,tip.number=numTips,model="yule")
#' aldousTrees <- rtreeshape(numTrees,tip.number=numTips,model="aldous")
#' biasedTrees <- rtreeshape(numTrees,tip.number=numTips,model="biased")
#'
#' # place into a list
#' allTrees <- list(pdaTrees = pdaTrees, yuleTrees = yuleTrees, aldousTrees = aldousTrees, biasedTrees, biasedTrees)
#'
#' # convert to phylo type
#' allTrees <- lapply(unlist(allTrees, recursive=FALSE), as.phylo)
#'
#' # construct the 2d MDS on the trees using the "canberra" distance method
#' mds2D <- coefficientMds(trees = allTrees, method = "canberra", dim = 2)
#'
#' # first construct the coefficient matrices and use them to create two MDS plots
#' coeffMats <- coefficientMatrix(allTrees)
#' mds2dTwoStep <- coefficientMds(coefficientMatrices = coeffMats, method = "canberra", dim = 2)
#' mds3dTwoStep <- coefficientMds(coefficientMatrices = coeffMats, method = "canberra", dim = 3)
#'
#' @export
coefficientMds <- function(trees, coefficientMatrices, method = "canberra", dim = 2) {
  # check what data was supplied
  if (missing(coefficientMatrices)) {
    treeLab <- names(trees)
    coeffMats <- coefficientMatrix(trees)

    distMatrix <- as.matrix(coefficientDist(coeffMats, method))
  } else if (missing(trees)) {
    treeLab <- names(coefficientMatrices)
    distMatrix <- as.matrix(coefficientDist(coefficientMatrices, method))
  } else {
    # error no supplied data
  }
  if (dim == 2) {
    fit <- cmdscale(distMatrix, k = 2)
    if (!is.null(treeLab)) {
      methodResults <- data.frame("method" = gsub("[[:digit:]]+", "", treeLab), fit)
    } else {
      mdsRes <- data.frame(fit)
    }

    ggplot(methodResults, aes(x = X1, y = X2, shape = method, color = method)) +
      geom_point() +
      ggtitle(method) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    fit <- cmdscale(distMatrix, k = 3) # k is the number of dim

    if (!is.null(treeLab)) {
      mdsRes <- data.frame("method" = gsub("[[:digit:]]+", "", treeLab), fit)
    } else {
      mdsRes <- data.frame(fit)
    }

    plot_ly(mdsRes, x = mdsRes$X1, y = mdsRes$X2, z = mdsRes$X3, color = mdsRes$method, type = "scatter3d", size = 2, mode = "markers", showlegend = TRUE) %>%
      layout(
        title = method
      )
  }
}

#' Calculates the distance matrix of multiple coefficient matrices
#'
#'
#' @param coefficientMatrices list of coefficient matrices
#' @param method method to use when calculating coefficient distances
#' @return distance matrix calculated from argument coefficient matrices
#' @import Matrix
#' @importFrom utils combn
#' @import pbapply
#' @examples
#' # example goes here
#' @export
coefficientDist <- function(coefficientMatrices, method = "sumLogDiff") {

  # determine and define the places to calculate in the distance matrix
  upperDist <- t(combn(length(coefficientMatrices), 2))
  distMat <- matrix(data = 0, nrow = length(coefficientMatrices), ncol = length(coefficientMatrices))

  # set the max times the progress bar is updated
  pboptions(nout = 50)

  isDouble <- is.double(coefficientMatrices[[1]])
  isComplex <- is.complex(coefficientMatrices[[1]])

  # determine if coefficients are real or imaginary
  if (isDouble && method == "sumLogDiff") {
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      sumLogDiff(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    })
  } else if (isComplex && method == "sumLogDiff") {
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      sumLogDiffComplex(coefficientMatrices[[i[1]]], coefficientMatrices[[i[2]]])
    })
  } else if (isComplex || isDouble && method == "fager") {
    coefficientMatrices <- lapply(coefficientMatrices, as.logical)
    distMat[upperDist] <- pbapply(upperDist, MARGIN = 1, FUN = function(i) {
      coeffMatA <- coefficientMatrices[[i[1]]]
      coeffMatB <- coefficientMatrices[[i[2]]]

      c <- as.numeric(length(coeffMatA[coeffMatA == TRUE & coeffMatB == FALSE]))
      b <- as.numeric(length(coeffMatA[coeffMatA == FALSE & coeffMatB == TRUE]))
      d <- as.numeric(length(coeffMatA[coeffMatA == FALSE & coeffMatB == FALSE]))
      a <- as.numeric(length(coeffMatA[coeffMatA == TRUE & coeffMatB == TRUE]))

      return((a/sqrt((a+b)*(a+c))) - sqrt(a+c)/2)

    })
  } else {
    # error
  }

  distMat <- t(distMat) + distMat

  return(distMat)
}

#' Builds a 2d or 3d mds plot from tree or coefficient matrix data
#'
#' Either calculates or recieves coefficient matrix data of phylo objects and uses a distance method
#' on the coefficients to construct an mds plot.
#' @param trees list of phylo objects used to calculate coefficient matrices
#' @param coefficientMatrices list of coefficient matrices if not passing list of trees
#' @param method method to use when calculating coefficient distances
#' @param dim number of dims of resulting mds plot
#' @return 2d or 3d MDS plot
#' @import Matrix
#' @rawNamespace import(plotly, except = last_plot)
#' @import ggplot2
#' @examples
#' library(apTreeshape)
#' # generate some forests
#' numTrees <- 10
#' numTips <- 500
#'
#' pdaTrees <- rtreeshape(numTrees, tip.number = numTips, model = "pda")
#' yuleTrees <- rtreeshape(numTrees, tip.number = numTips, model = "yule")
#' aldousTrees <- rtreeshape(numTrees, tip.number = numTips, model = "aldous")
#' biasedTrees <- rtreeshape(numTrees, tip.number = numTips, model = "biased")
#'
#' # place into a list
#' allTrees <- list(pdaTrees = pdaTrees, yuleTrees = yuleTrees, aldousTrees = aldousTrees, biasedTrees, biasedTrees)
#'
#' # convert to phylo type
#' allTrees <- lapply(unlist(allTrees, recursive = FALSE), as.phylo)
#'
#' # construct the 2d MDS on the trees using the "canberra" distance method
#' mds2D <- coefficientMds(trees = allTrees, method = "canberra", dim = 2)
#'
#' # first construct the coefficient matrices and use them to create two MDS plots
#' coeffMats <- coefficientMatrix(allTrees)
#' mds2dTwoStep <- coefficientMds(coefficientMatrices = coeffMats, method = "canberra", dim = 2)
#' mds3dTwoStep <- coefficientMds(coefficientMatrices = coeffMats, method = "canberra", dim = 3)
#' @export
coefficientMds <- function(trees, coefficientMatrices, method = "sumLogDiff", dim = 2, title = method, legendTitle = " ") {
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
      methodResults <- data.frame("method" = treeLab, fit, "color" = treeLab)
    } else {
      methodResults <- data.frame(fit)
    }

    ggplot(methodResults, aes(x = X1, y = X2, color = color)) +
      geom_point() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(color = legendTitle)
  } else {
    fit <- cmdscale(distMatrix, k = 3) # k is the number of dim

    if (!is.null(treeLab)) {
      mdsRes <- data.frame("method" = treeLab, fit, "color" = treeLab)
    } else {
      mdsRes <- data.frame(fit)
    }

    plot_ly(mdsRes, x = mdsRes$X1, y = mdsRes$X2, z = mdsRes$X3, color = mdsRes$color, type = "scatter3d", size = 2, mode = "markers", showlegend = TRUE) %>%
      add_annotations(
        text = legendTitle, xref = "paper", yref = "paper",
        x = 1.02, xanchor = "left",
        y = 0.8, yanchor = "bottom", # Same y as legend below
        legendtitle = TRUE, showarrow = FALSE
      ) %>%
      layout(
        legend = list(y = 0.8, yanchor = "top"),
        title = title
      )
  }
}

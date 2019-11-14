#' Convert trees to coefficient matrices
#'
#' Converts rooted full binary trees to tree distinguishing polynomials described with coefficient matrices.
#'
#' @param trees a single phylo object or a list of phylo objects
#' @param type one of:
#' \describe{
#'   \item{\dQuote{real}}{real coefficient polynomials}
#'   \item{\dQuote{complex}}{complex coefficient polynomials (real polynomials with y = 1 + i)}
#'   \item{\dQuote{tipLabel}}{complex coefficient polynomial that utilize binary trait tip labels on the phylo objects}
#' }
#' @return the resulting coefficient matrix or matrices of the form:
#' \describe{
#'   \item{\dQuote{real}}{a real matrix where the ith row, jth column represents the x^(j-1)*y^(i-1) coefficient}
#'   \item{\dQuote{complex}}{a complex vector where the kth column represents the x^(k-1) coefficient}
#'   \item{\dQuote{tipLabel}}{given trees with two unique tip labels \dQuote{a}, \dQuote{b} a complex matrix where the ith row, jth column represents the a^(i-1)*b^(j-1) coefficient}
#' }
#' @param numThreads number of threads to be used, the default (-1) will use the number of cores in the machine and numThreads = 0 will only use the main thread
#' @importFrom ape as.phylo
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#' library(treenomial)
#' library(ape)
#'
#' # generate a tree
#' tree <- rtree(n = 30, rooted = TRUE)
#'
#' # a real coefficient matrix
#' treeToPoly(tree, numThreads = 0)
#'
#' # complex coefficient vector for the tree
#' treeToPoly(tree, type = "complex", numThreads = 0)
#'
#' # for a list of trees
#' treeToPoly(rmtree(10, 30), numThreads = 0)
#'
#' @export
treeToPoly <- function(trees, type = "real", numThreads = -1) {

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


  if (type == "real" || type == "complex") {
    if (singleTree) {
      inds <- unique(rev(trees$edge[trees$edge[, 1] >= length(trees$tip.label), ]))
      wedgeOrder <- ifelse(inds <= length(trees$tip.label), "0", "1")

      res <- coeffMatList(list(wedgeOrder), type = type, nThreads = numThreads)
      attributes(res) <- NULL
    } else {
      trees <- lapply(trees, function(x) {
        inds <- unique(rev(x$edge[x$edge[, 1] >= length(x$tip.label), ]))
        ifelse(inds <= length(x$tip.label), "0", "1")
      })

      res <- coeffMatList(trees, type = type, nThreads = numThreads)

      attributes(res) <- NULL
      if(!is.null(names(trees))){
        names(res) <- names(trees)
      }
    }
  } else if (type == "tipLabel") {
    if (singleTree) {
      uniqueTipLab <- sort(unique(trees$tip.label))


      if (length(uniqueTipLab) == 1) {
        uniqueTipLab <- rep(uniqueTipLab, 2)
      } else if (length(uniqueTipLab) > 2) {
        stop("only phylo trees with two unique tip labels are currently supported")
      }

      inds <- unique(rev(trees$edge[trees$edge[, 1] >= length(trees$tip.label), ]))

      wedgeOrder <-  ifelse(inds <= length(trees$tip.label), trees$tip.label[inds], "1")

      res <- coeffMatList(list(wedgeOrder), type = type, tipLabA = uniqueTipLab[[1]], tipLabB = uniqueTipLab[[2]], nThreads = numThreads)
      attributes(res) <- NULL
    } else {
      uniqueTipLabFirst <- sort(unique(trees[[1]]$tip.label))

      if (length(uniqueTipLabFirst) == 1) {
        uniqueTipLabFirst <- rep(uniqueTipLabFirst, 2)
      } else if (length(uniqueTipLabFirst) > 2) {
        stop("only phylo trees with two unique tip labels are currently supported")
      }

      lapply(trees, function(i){
        uniqueTipLab <- unique(i$tip.label)

        if (length(uniqueTipLab) == 1) {
          uniqueTipLab <- rep(uniqueTipLab, 2)
        } else if (length(uniqueTipLab) > 2) {
          stop("only phylo trees with two unique tip labels are currently supported")
        }
      })

      wedgeOrders <- lapply(1:length(trees), function(i) {
        x <- trees[[i]]

        inds <- unique(rev(x$edge[x$edge[, 1] >= length(x$tip.label), ]))
        ifelse(inds <= length(x$tip.label), x$tip.label[inds], "1")
      })

      res <- coeffMatList(wedgeOrders, type = type, tipLabA = uniqueTipLabFirst[[1]], tipLabB = uniqueTipLabFirst[[2]], nThreads = numThreads)

      attributes(res) <- NULL
      if(!is.null(names(trees))){
        names(res) <- names(trees)
      }
    }
  } else {
    stop('incorrect type, available types are "real", "complex" or "tipLabel" ')
  }

  if(singleTree){
    return(res[[1]])
  } else {
    return(res)
  }

}


#' Align various types of coefficient matrices
#'
#' @param coefficientMatrices a list of coefficient matrices of various sizes
#' @return the aligned list of coefficient matrices
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @details
#' Alignment depends on the type of matrix:
#' \describe{
#'   \item{\dQuote{real}}{the smaller matrices columns are prepended with zero columns to align with the max number of columns and the rows are appended with zero rows to match the max number of rows}
#'   \item{\dQuote{complex}}{the smaller vectors are appended with zeroes to match the max length vector}
#'   \item{\dQuote{tipLabel}}{the smaller matrices are appended with zeroes to match the max number of rows and columns}
#' }
#' @examples
#'
#' library(treenomial)
#' library(ape)
#' differentSizeTrees <- c(rtree(2), rmtree(10,10))
#' coeffs <- treeToPoly(differentSizeTrees, numThreads = 0)
#' alignedCoeffs <- alignPoly(coeffs)
#'
#'
#' @export
alignPoly <- function(coefficientMatrices){
  # check input arguments
  if(class(coefficientMatrices) == "list"){
    if(class(coefficientMatrices[[1]]) == "matrix"){
      if(typeof(coefficientMatrices[[1]]) == "double"){
        coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "real")
      } else if(typeof(coefficientMatrices[[1]]) == "complex"){
        if(dim(coefficientMatrices[[1]])[[1]] == 1){
          coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "complex")
        } else {
          coefficientMatrices <- alignCoeffs(coefficientMatrices, type = "tipLabel")
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


  return(coefficientMatrices)
}

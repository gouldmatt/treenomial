#' Convert trees to coefficient matrices
#'
#' Converts trees to tree distinguishing polynomials descibed with coefficient matrices.
#'
#' @param trees an object or list of class "phylo"
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
#' treeToPoly(tree)
#'
#' # complex coefficient vector for the tree
#' treeToPoly(tree, type = "complex")
#'
#' # for a list of trees
#' treeToPoly(rmtree(10, 30))
#'
#' @export
treeToPoly <- function(trees, type = "real") {

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

      res <- coeffMatListCpp(list(wedgeOrder), type = type)
      attributes(res) <- NULL
    } else {
      trees <- lapply(trees, function(x) {
        inds <- unique(rev(x$edge[x$edge[, 1] >= length(x$tip.label), ]))
        ifelse(inds <= length(x$tip.label), "0", "1")
      })

      res <- coeffMatListCpp(trees, type = type)

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

      res <- coeffMatListCpp(list(wedgeOrder), type = type, tipLabA = uniqueTipLab[[1]], tipLabB = uniqueTipLab[[2]])
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

      res <- coeffMatListCpp(wedgeOrders, type = type, tipLabA = uniqueTipLabFirst[[1]], tipLabB = uniqueTipLabFirst[[2]])

      attributes(res) <- NULL
      if(!is.null(names(trees))){
        names(res) <- names(trees)
      }
    }
  } else {
    stop('incorrect type, available types are "real", "complex" or "tipLabel" ')
  }

  return(res)
}


#' Align varios types of coefficient matrices
#'
#' @param coefficientMatrices a list of coefficient matrices of various sizes
#' @return the aligned list of coefficient matrices
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#' library(treenomial)
#' library(ape)
#' differentSizeTrees <- c(rtree(2), rmtree(10,10))
#' coeffs <- treeToPoly(differentSizeTrees)
#' alignedCoeffs <- alignPoly(coeffs)
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

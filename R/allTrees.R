#' Calculate all full unordered m-ary trees up to n tips
#'
#' Return real coefficient matrices, complex coefficient vectors, or phylo objects for all possible unordered full m-ary trees up to n tips.
#' For binary trees (m = 2), the number of trees at each number of tips follows the \href{https://oeis.org/A001190}{Wedderburn-Etherington numbers}.
#' @param m max number of children for each node
#' @param n max number of tips
#' @param type one of:
#' \describe{
#'   \item{\dQuote{real}}{real coefficient polynomials}
#'   \item{\dQuote{complex}}{complex coefficient polynomials (real polynomials with y = 1 + i)}
#'   \item{\dQuote{phylo}}{phylo objects}
#' }
#' @inheritParams treeToPoly
#' @return list of lists containing all the trees in \strong{type} format for each number of tips
#' @note only m = 2 is currently supported
#' @useDynLib treenomial
#' @examples
#' \donttest{
#' library(treenomial)
#' library(ape)
#'
#' # generate coefficient matrices describing the polynomials of all possible
#' # unordered full binary trees up to 10 tips
#'
#' allBinTenRealCoeff <- allTrees(10, type = "phylo")
#' lengths(allBinTenRealCoeff)
#'
#' # phylo type example, plot all 6 tip unordered full binary trees
#'
#' allBinSixPhylo <- allTrees(6, type = "phylo")[[6]]
#' par(mfrow=c(1,6))
#' plots <- lapply(allBinSixPhylo, function(t){
#'   plot.phylo(ladderize(t), direction = "downwards", show.tip.label = FALSE)
#' })
#' }
#' @export
allTrees <- function(n, m = 2, type = "real") {
  # wadNum <- c(1,1,1,2,3,6, 11, 23, 46, 98, 207, 451,
  #             983, 2179, 4850, 10905, 24631, 56011, 127912,
  #             293547, 676157, 1563372, 3626149, 8436379, 19680277,
  #             46026618, 107890609)
  if(m != 2){
    stop("only binary (m = 2) trees are currently supported")
  }

  if(type == "real"){
    allBinaryTreeShapesReal(n)
  } else if(type == "complex"){
    allBinaryTreeShapesComplex(n)
  } else if(type == "phylo"){
    wedgeList <- allBinaryTreeShapesPhylo(n)
    wedgeList <- wedgeList[-1]

    nTipsTrees <- vector("list", length = n)
    nTipsTrees <- lapply(1:n, function(i){ nTipsTrees[[i]] <- list()})

    nTipsTrees[[1]][[1]] <- "leaf"

    for (i in 1:length(wedgeList)) {
      for (j in 1:length(wedgeList[[i]])) {
        wedgeListCurr <- wedgeList[[i]][[j]]
        nTipsTrees[[i+1]][[j]] <- wedge(nTipsTrees[[wedgeListCurr[1]]][[wedgeListCurr[2]]],nTipsTrees[[wedgeListCurr[3]]][[wedgeListCurr[4]]])
      }
    }

    return(nTipsTrees)
  } else {
    stop("invalid type")
  }
}

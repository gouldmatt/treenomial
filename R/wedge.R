#' Performs the wedge operation
#'
#' Calculates the result from the wedge operation on two real coefficient
#' matrices, two complex coefficient vectors or two phylo objects.
#'
#' @param A,B two real coefficient matrices, complex coefficient vectors or phylo objects
#' @return the wedge result in the same form as the arguments
#' @import ape
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#'
#' library(treenomial)
#' library(ape)
#'
#' # wedge two real coefficient matrices
#'
#' leaf <- matrix(c(0,1), nrow = 1, ncol = 2)
#' wedge(leaf, leaf)
#'
#' # wedge two complex coefficient vectors
#'
#' leaf <- as.complex(c(0,1))
#' wedge(leaf, leaf)
#'
#' @export
wedge <- function(A, B) {

  type = class(A[1])

  if (type == "numeric" && class(A) == "matrix") {
    wedgeExport(A,B)
  } else if (type == "complex") {
    as.vector(wedgeExportConv(A,B))
  } else {
    if(all(A == "leaf" & B == "leaf")){
      res <- rtree(2, br = NULL)
    } else if (is.character(A)){
      res <- bind.tree(rtree(2, br = NULL), B, where = 1)
    } else if (is.character(B)){
      res <- bind.tree(rtree(2, br = NULL), A, where = 1)
    } else {

      A$edge.length <- NULL
      B$edge.length <- NULL

      # bind first tree with a cherry
      res <- bind.tree(rtree(2, br = NULL), A, where = 1)

      # bind second tree on to other tip of cherry
      res <- bind.tree(res, B, where = 1)
    }
    res$tip.label <- 1:(length(res$tip.label))
    return(res)
  }
}

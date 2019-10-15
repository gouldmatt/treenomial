#' Performs the wedge operation
#'
#' Calculates the result from the wedge operation on two real coefficient
#' matrices (class "dgCMatrix"), two complex coefficient vectors or two phylo objects.
#'
#'
#' @param A,B Two real coefficient matrices, complex coefficient vectors or phylo objects.
#' @param type The type of objects to wedge: "real", "complex" or "phylo".
#' @return The wedge result.
#' @import Matrix
#' @importFrom pracma conv
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#' require(Matrix)
#' require(ape)
#'
#' # wedge two real coefficient matrices
#' leaf <- sparseMatrix(1, 2, x = 1)
#' wedge(leaf, leaf, type = "real")
#'
#' # wedge two complex coefficient matrices
#' leaf <- as.complex(c(0,1))
#' wedge(leaf, leaf, type = "complex")
#'
#' # wedge two phylo objects
#' cherry <- rtree(2, br = NULL)
#' wedge(cherry, cherry, type = "phylo")
#'
#' @export
wedge <- function(A, B, type = "real") {

    if (type == "real") {
      # find the size of the resulting polynomail matrix
      colsA <- ncol(A)
      colsB <- ncol(B)
      newMatrixSizeCol <- colsB + colsA - 1
      newMatrixSizeRow <- ceiling(newMatrixSizeCol / 2)

      # convert operands from sparse matrix to coordinate list matrix => [row, col, val]
      A <- as.matrix(unname(Matrix::summary(A)))
      B <- as.matrix(unname(Matrix::summary(B)))

      # perform the wedge operation on the elements of A and B filling resPolyMat
      res <- matrix(data = 0, nrow = newMatrixSizeRow, ncol = newMatrixSizeCol)
      wedgeFill(A, B, res)

      # convert result to sparse matrix
      as(res, "sparseMatrix")
    } else if (type == "complex") {
      res <- conv(A, B)
      res[1] <- res[1] + 1i + 1
      res
    } else if (type == "phylo") {

      if(A == "leaf" && B == "leaf"){
        rtree(2, br = NULL)
      } else if (is.character(A)){
        bind.tree(rtree(2, br = NULL), B, where = 1)
      } else if (is.character(B)){
        bind.tree(rtree(2, br = NULL), A, where = 1)
      } else {
        # bind first tree with a cherry
        res <- bind.tree(rtree(2, br = NULL), A, where = 1)

        # bind second tree on to other tip of cherry
        bind.tree(res, B, where = 1)
      }
    } else {
      # error missing wedge type
      stop("missing or incorrect wedge type")
    }
}



tipLabelWedge <- function(coefficientMatrixA, coefficientMatrixB) {

  newCols <- ncol(coefficientMatrixA) + ncol(coefficientMatrixB) - 1
  newRows <- nrow(coefficientMatrixB) + nrow(coefficientMatrixA) - 1

  # convert operands from sparse matrix to coordinate list matrix => [row, col, val]
  nonzeroA <- unname(which(coefficientMatrixA != 0, arr.ind = T))
  nonzeroB <- unname(which(coefficientMatrixB != 0, arr.ind = T))

  coefficientMatrixA <- as.matrix(cbind(nonzeroA,coefficientMatrixA[nonzeroA]))
  coefficientMatrixB <- as.matrix(cbind(nonzeroB,coefficientMatrixB[nonzeroB]))

  # perform the wedge operation on the elements of coefficientMatrixA and coefficientMatrixB filling resPolyMat
  resPolyMat <- matrix(data = 0i, nrow = newRows, ncol = newCols)
  wedgeFillComplex(coefficientMatrixA, coefficientMatrixB, resPolyMat)

  # convert result to sparse matrix
  # resPolyMat <- as(resPolyMat, "sparseMatrix")

  return(resPolyMat)
}

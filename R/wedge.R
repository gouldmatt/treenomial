#' Performs the wedge operation on two coefficient matrices
#'
#' Calculates the result from the wedge operation of two coefficient matrices.
#' @param coefficientMatrixA First tree coefficient matrix
#' @param coefficientMatrixB Second tree coefficient matrix
#' @return The resulting polynomail matrix.
#' @import Matrix
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#' leaf <- sparseMatrix(1, 2, x = 1)
#' wedge(leaf, leaf)
#' [1,] . . 1
#' [2,] 1 . .
#' @export
wedge <- function(coefficientMatrixA, coefficientMatrixB) {

  # find the size of the resulting polynomail matrix
  colsA <- ncol(coefficientMatrixA)
  colsB <- ncol(coefficientMatrixB)
  newMatrixSizeCol <- colsB + colsA - 1
  newMatrixSizeRow <- ceiling(newMatrixSizeCol / 2)

  # convert operands from sparse matrix to coordinate list matrix => [row, col, val]
  coefficientMatrixA <- as.matrix(unname(Matrix::summary(coefficientMatrixA)))
  coefficientMatrixB <- as.matrix(unname(Matrix::summary(coefficientMatrixB)))

  # perform the wedge operation on the elements of coefficientMatrixA and coefficientMatrixB filling resPolyMat
  resPolyMat <- matrix(data = 0, nrow = newMatrixSizeRow, ncol = newMatrixSizeCol)
  wedgeFill(coefficientMatrixA, coefficientMatrixB, resPolyMat)

  # convert result to sparse matrix
  resPolyMat <- as(resPolyMat, "sparseMatrix")

  return(resPolyMat)
}

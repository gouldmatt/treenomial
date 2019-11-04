#' Plots a Julia Set for a tree
#'
#' Finds the Julia Set for the complex polynomial of a tree and plots in a square image.
#'
#' @param tree phylo object
#' @param pixelLength integer of pixels on one side of the image
#' @param center complex number giving the center of the image on the complex plane
#' @param maxZ integer giving the max value for the real and imaginary axis in the image
#' @param asp aspect ratio of the image
#' @param col colours to be used for the image
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @examples
#'
#' library(treenomial)
#' library(ape)
#' treeJuliaSet(stree(5,type = "right"))
#'
#' @export
treeJuliaSet <- function(tree, pixelLength = 1000, center = 0, maxZ = 2, asp = 1, col = c("white", colorRampPalette(c("dodgerblue4", "lightblue"))(98) , "black")) {
  parBackup <- par()

  par(mar = c(1,1,1,1))

  coeff <- treeToPoly(tree, type = "complex")

  res <- juliaSet(coeffs = coeff[[1]], pixelLength = pixelLength, center = center, maxZ = maxZ)

  graphics::image(z = log(res), axes = FALSE, asp = asp, col = col)

  par <- parBackup
}

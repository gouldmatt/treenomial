#' Estimates R0
#'
#' @param tree list of coefficient matrices
#' @import TreeSim
#' @return estimated R0
#' # example goes here
#' @export
estimateR0 <- function(tree, start = 1.5, end = 8.5, step = 1, numSim = 100){

  testRange <- seq(from = start, to = end, by = step)

  # generate a range of b-d trees
  simbds=lapply(testRange, function(x) sim.bd.taxa(n=2*(tree$Nnode + 1), numbsim=numSim,mu=1, lambda=x,complete = TRUE))


  R0=unlist(lapply(testRange, function(x) rep(x,numSim))) # R0 values (lambda/mu)

  psim <- unlist(simbds, recursive = FALSE) # put into one list

  psim=lapply(psim, function(oneTree) prunetosize(oneTree,(tree$Nnode + 1))) # uniformly at random drop

  names(psim) <- R0

  # generate the polynomials
  bdTreesCoeff <- coefficientMatrix(psim, complex = FALSE)
  treeCoeff <- coefficientMatrix(tree)
  allCeoff <- c(bdTreesCoeff,"tree"= list(treeCoeff))

  distanceMatrix <- coefficientDist(coefficientMatrices = allCeoff)
  # find the 3 closest trees
  approxTree1 <- which(distanceMatrix[nrow(distanceMatrix),] == min(distanceMatrix[nrow(distanceMatrix),1:(nrow(distanceMatrix)-1)]))
  distanceTree1 <- distanceMatrix[nrow(distanceMatrix),approxTree1]
  distanceMatrix[nrow(distanceMatrix),approxTree1] <- Inf

  approxTree2 <- which(distanceMatrix[nrow(distanceMatrix),] == min(distanceMatrix[nrow(distanceMatrix),1:(nrow(distanceMatrix)-1)]))
  distanceTree2 <- distanceMatrix[nrow(distanceMatrix),approxTree2]
  distanceMatrix[nrow(distanceMatrix),approxTree2] <- Inf

  approxTree3 <- which(distanceMatrix[nrow(distanceMatrix),] == min(distanceMatrix[nrow(distanceMatrix),1:(nrow(distanceMatrix)-1)]))
  distanceTree3 <- distanceMatrix[nrow(distanceMatrix),approxTree3]
  distanceMatrix[nrow(distanceMatrix),approxTree3] <- Inf

  approxTrees <- list(R0 = as.numeric(names(psim[c(approxTree1,approxTree2,approxTree3)])), distance = c(distanceTree1,distanceTree2,distanceTree3))
  print(approxTrees$R0)
  approxR0 <- sum(approxTrees$R0*(1/approxTrees$distance))/sum(1/approxTrees$distance)

  return(approxR0)
}


# # pruning function
prunetosize <- function(tree,size) {
  nt=length(tree$tip.label)
  if (nt <= size) {return(tree)} else {
    todrop =sample(1:nt, nt-size);
    return(drop.tip(tree, todrop))
  }

}

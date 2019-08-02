#' Estimates R0
#'
#' @param tree list of coefficient matrices
#' @import TreeSim
#' @return estimated R0
#' # example goes here
#' @export
estimateR0 <- function(tree, start = 1.5, end = 8.5, step = 1, numStepSim = 100, numTotalSim = 10){

  approxR0 <- vector(mode = "numeric", length = numTotalSim)

  pb <- startpb(0, numTotalSim)
  on.exit(closepb(pb))

  for(i in 1:numTotalSim){

    testRange <- seq(from = start, to = end, by = step)

    # generate a range of b-d trees
    simbds=lapply(testRange, function(x) sim.bd.taxa(n=2*(tree$Nnode + 1), numbsim=numStepSim,mu=1, lambda=x,complete = TRUE))


    R0=unlist(lapply(testRange, function(x) rep(x,numStepSim))) # R0 values (lambda/mu)

    usim <- unlist(simbds, recursive = FALSE) # put into one list

    ppsim <- lapply(usim, function(oneTree) prunetosize(oneTree,(tree$Nnode + 1))) # uniformly at random drop

    names(ppsim) <- R0

    # generate the polynomials
    bdTreesCoeff <- coefficientMatrix(ppsim, complex = FALSE)
    treeCoeff <- coefficientMatrix(tree, complex = FALSE)
    allCeoff <- c(bdTreesCoeff,"tree"= list(treeCoeff))

    distanceMatrix <- coefficientDist(coefficientMatrices = allCeoff, method = "fager")

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

    approxTrees <- list(R0 = as.numeric(names(ppsim[c(approxTree1,approxTree2,approxTree3)])), distance = c(distanceTree1,distanceTree2,distanceTree3))

    approxR0[i] <- sum(approxTrees$R0*(1/approxTrees$distance))/sum(1/approxTrees$distance)
    print(approxTrees$R0)
    setpb(pb, i)
  }

  return(approxR0)
}

#' Estimates p
#'
#' @param tree list of coefficient matrices
#' @import TreeSim
#' @return estimated p
#' # example goes here
#' @export
estimateP <- function(tree, start = 1.5, end = 8.5, step = 1, numStepSim = 100, numTotalSim = 10){

  approxP <- vector(mode = "numeric", length = numTotalSim)

  pb <- startpb(0, numTotalSim)
  on.exit(closepb(pb))

  for(i in 1:numTotalSim){
    # generate the sample
    biasedTreesSample <- rtreeshape(numStepSim,tip.number=(targetTree$Nnode + 1), model="biased", p=p)

    # generate the polynomials
    biasedTreesCoeff <- coefficientMatrix(biasedTreesSample, complex = TRUE)
    treeCoeff <- coefficientMatrix(tree, complex = TRUE)
    allCeoff <- c(biasedTreesCoeff,"tree"= list(treeCoeff))

    distanceMatrix <- coefficientDist(coefficientMatrices = allCeoff, method = "sumLogDiff")

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

    approxTrees <- list(R0 = as.numeric(names(biasedTreesSample[c(approxTree1,approxTree2,approxTree3)])), distance = c(distanceTree1,distanceTree2,distanceTree3))

    approxP[i] <- sum(approxTrees$R0*(1/approxTrees$distance))/sum(1/approxTrees$distance)
    print(approxTrees$R0)
    setpb(pb, i)
  }

  return(approxP)
}


#' @export
bdTreesDistance <- function(targetTree, r0, sampleSize){
  print("generating sample")
  # simulate sample size number of trees with argument R0 and twice the # of tips as the target tree
  unprunedSampleTrees <- sim.bd.taxa(n= 2*(targetTree$Nnode + 1), numbsim= sampleSize, mu=1, lambda= r0,complete = TRUE)

  # prune the sample trees
  prunedSampleTrees <- lapply(unprunedSampleTrees, function(x) prunetosize(x,(targetTree$Nnode + 1)))

  print("generating coefficients")
  # generate the coefficient matrices for the sample and the target
  pboptions(type = "txt")
  bdTreesCoeff <- coefficientMatrix(prunedSampleTrees, complex = FALSE)
  treeCoeff <- coefficientMatrix(targetTree, complex = FALSE)

  print("finding distances to sample")
  # find the distance between each sample and the target
  distances <- vector(mode = "numeric", length = sampleSize)
  for(i in 1:sampleSize){
    distances[i] <- sumLogDiff(treeCoeff,bdTreesCoeff[[i]])
  }

  # printout and return average
  print(c("R0 value:", r0))
  print(c("average",sum(distances)/sampleSize))
  print(c("std dev",std(distances)))

  return(sum(distances)/sampleSize)
}

#' @export
biasedTreesDistance <- function(p, targetTree, sampleSize){

  if(p < 0 || p > 1){
    return(Inf)
  }

  # convert target tree to phylo
  targetTree <- as.phylo(targetTree[[1]])

  # generate the sample
  biasedTreesSample <- rtreeshape(sampleSize,tip.number=(targetTree$Nnode + 1), model="biased", p=p)

  # convert to phylo type
  biasedTreesSample <- lapply(biasedTreesSample, as.phylo)

  # generate the coeffs
  biasedTreesCoeff <- coefficientMatrix(biasedTreesSample, complex = FALSE)
  treeCoeff <- coefficientMatrix(targetTree, complex = FALSE)

  # compute the distances
  distances <- vector(mode = "numeric", length = sampleSize)
  for(i in 1:sampleSize){
    distances[i] <- sumLogDiff(treeCoeff,biasedTreesCoeff[[i]])
  }

  print(c("p value:", p))
  print(c("average",sum(distances)/sampleSize))
  print(c("std dev",std(distances)))

  return(sum(distances)/sampleSize)
}


# # pruning function
prunetosize <- function(tree,size) {
  nt=length(tree$tip.label)
  if (nt <= size) {return(tree)} else {
    todrop =sample(1:nt, nt-size);
    return(drop.tip(tree, todrop))
  }

}

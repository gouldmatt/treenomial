#' Calculates the distance between two tree lattices
#'
#'
#' @param L1,L2 lattices to calculate distance between
#' @param w the amount to weight absent positions in each lattice in the distance
#' @return distance value
#' @examples
#'
#' library(treenomial)
#' library(ape)
#'
#'
#' L1 <- treeToLattice(rtree(10), numThreads = 0)
#' L2 <- treeToLattice(rtree(10), numThreads = 0)
#' d <- lattDist(L1, L2, w = 1/4)
#'
#'
#' @export
lattDist <- function(L1,L2, w = 1/4){

  L1 <- as.matrix(unname(L1))

  L2 <- as.matrix(unname(L2))


  if(ncol(L1) > ncol(L2)){
    L2 <- cbind(L2,matrix(0, nrow = nrow(L2),ncol = ncol(L1) - ncol(L2)))

  } else if(ncol(L1) < ncol(L2)){
    L1 <- cbind(L1,matrix(0, nrow = nrow(L1),ncol = ncol(L2) - ncol(L1)))

  }


  return(lattDistance(L1,L2, w = w))

}

#' Calculates the distance matrix from a list of tree lattices
#'
#'
#' @param lattices list of tree lattices
#' @param w the amount to weight absent positions in each lattice in the distance
#' @param numThreads number of threads to be used, the default (-1) will use the number of cores in the machine and numThreads = 0 will only use the main thread
#' @return distance matrix calculated from argument tree lattices
#' @examples
#'
#' library(treenomial)
#' library(ape)
#'
#' # lattices for ten trees of 10 tips
#' latts <- treeToLattice(rmtree(10, 10), numThreads = 0)
#'
#' # distance matrix from the list of lattices
#' d <- lattToDistMat(latts, numThreads = 0)
#'
#' @export
lattToDistMat <- function(lattices, w = 1/4, numThreads = -1) {
  if(!is(lattices, "list")){
    stop("input must be a list of tree lattices")
  }

  distMat <- lattDistMat(lattices, w = w, nThreads = numThreads)

  rownames(distMat) <- names(lattices)
  colnames(distMat) <- names(lattices)

  return(distMat)
}


#' Convert trees to lattice form
#'
#' Converts rooted full binary trees to a tree lattice using tree distinguishing polynomials.
#'
#' @param trees a single phylo object or a list of phylo objects
#' @return the resulting lattice or lattices as matrices with columns:
#' \describe{
#'   \item{\dQuote{node}}{the node number}
#'   \item{\dQuote{pl}}{the parent of this node}
#'   \item{\dQuote{lattice}}{the position of this node in the tree lattice}
#'   \item{\dQuote{bl}}{branch length to this node from its parent}
#'   \item{\dQuote{depth}}{depth of this node in the tree}
#' }
#' @param numThreads number of threads to be used, the default (-1) will use the number of cores in the machine and numThreads = 0 will only use the main thread
#' @importFrom ape as.phylo
#' @useDynLib treenomial
#' @importFrom Rcpp sourceCpp
#' @importFrom methods is
#' @examples
#' library(treenomial)
#' library(ape)
#'
#' tree <- rtree(10)
#' lattice <- treeToLattice(tree, numThreads = 0)
#'
#' @export
treeToLattice <- function(trees, numThreads = -1) {


  # check input format
  if(!is(trees,"phylo") && !is(trees,"list") && !is(trees,"multiPhylo")){
    tryCatch({
      trees <- as.phylo(trees)
    }, error = function(e) {
      stop("incorrect input format, trees must be phylo or list of phylo objects")
    })
  } else if (!is(trees,"phylo") && !is(trees,"multiPhylo")){
    tryCatch({
      trees <- lapply(trees, as.phylo)
      singleTree <- FALSE
    }, error = function(e) {
      stop("incorrect input format, trees must be phylo or list of phylo objects")
    })
  }

  singleTree <- is(trees,"phylo")


  if(singleTree) trees <- list(trees)


  trees <- lapply(trees, reorder.phylo)

  wedgeOrdersNodes <- lapply(trees, function(x) {
    as.character(unique(rev(x$edge[x$edge[, 1] >= length(x$tip.label), ])))
  })

  wedgeOrders <- lapply(trees, function(x) {
    inds <- unique(rev(x$edge[x$edge[, 1] >= length(x$tip.label), ]))
    ifelse(inds <= length(x$tip.label), "0", "1")
  })

  tips <- sapply(trees, function(x) {
    xTips <- length(x$tip.label)
    if(xTips < 2){
      stop("Tree with less than two tips found!")
    }
    xTips
  })

  latList <- lapply(trees, function(x) {
    allocateLattice(x, length(x$tip.label))
  })


  latList <- latticeList(wedgeOrders, wedgeOrdersNodes, latList, tips, nThreads = numThreads)


  attributes(latList) <- NULL


  for (i in 1:length(latList)){
    # latList[[i]][,3] <- lattPositions[[i]]
    colnames(latList[[i]]) <- c("node","pl","bl","depth",rep(" ",ncol(latList[[i]])-4))
  }

  if(!is.null(names(trees))){
    names(latList) <- names(trees)
  }

  if(singleTree) latList <- latList[[1]]

  return(latList)

}



allocateLattice <- function(tree, tips){
  tree <- reorder.phylo(tree,"postorder")

  nodes <- 2*tips-1 # number of nodes
  ind <- 1:nodes
  pos <- rep(0,nodes)

  # parent list
  ans <- Ancestors(tree,type = "parent")
  pl <- cbind(ind,ans,pos)
  # root
  r <- pl[pl[,2]==0,1]

  # find the node depths
  tr<-tree
  tr$edge.length<-rep(1,nodes)
  dpt <- as.data.frame(cbind(nodeHeights(tr),tr$edge))
  # dpt <- arrange(dpt,V4)
  dpt <- dpt[order(dpt$V4),]

  dpt <- dpt[,2]
  dptn <- length(dpt)

  if(r<=dptn){
    dpt <- c(dpt[1:(r-1)],0,dpt[r:dptn])
  } else {
    dpt <- c(dpt[1:(r-1)],0)
  }

  h <- max(dpt) # height of the tree

  # branch length
  bl <- as.data.frame(cbind(tree$edge,tree$edge.length))
  # bl <- arrange(bl,V2)
  bl <- bl[order(bl$V2),]

  bl <- bl[,3]
  bln <- length(bl)

  if(r<=bln){
    bl <- c(bl[1:(r-1)],0,bl[r:bln])
  } else {
    bl <- c(bl[1:(r-1)],0)
  }

  # preallocating the lattice
  pl <- pl[,2]
  L <- as.data.frame(cbind(ind,pl,pos,bl,dpt))
  L[L[,2]==0,3] <- 1

  BinPos <- matrix(0L, nrow = nodes, ncol = h+1)
  BinL <- as.data.frame(cbind(ind,pl,bl,dpt,BinPos))

  BinL[BinL[,2]==0,5:(5+h)] <- c(1,rep(0,h))

  return(as.matrix(unname(BinL)))
}



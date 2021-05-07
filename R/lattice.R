
#' @export
lattDist <- function(L1,L2){

  L1 <- as.matrix(unname(L1))

  L2 <- as.matrix(unname(L2))

  return(lattDistance(L1,L2))

}
#' @export
lattToDistMat <- function(lattices, numThreads = -1) {
  if(!is(lattices, "list")){
    stop("input must be a list of tree lattices")
  }

  distMat <- lattDistMat(lattices, nThreads = numThreads)

  rownames(distMat) <- names(lattices)
  colnames(distMat) <- names(lattices)

  return(distMat)
}


#' @export
latticize <- function(trees, numThreads = -1) {


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


  lattPositions <- latticeList(wedgeOrders, wedgeOrdersNodes, latList, tips, nThreads = numThreads)


  attributes(lattPositions) <- NULL


  for (i in 1:length(latList)){
    latList[[i]][,3] <- lattPositions[[i]]
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

  # preallocating the ladderized form
  pl <- pl[,2]
  L <- as.data.frame(cbind(ind,pl,pos,bl,dpt))
  L[L[,2]==0,3] <- 1
  return(as.matrix(unname(L)))
}

########## From phangorn treeManipulation.R
Ancestors <- function (x, node, type = c("all", "parent")) {
  parents <- x$edge[, 1]
  child <- x$edge[, 2]
  pvector <- integer(max(x$edge))
  pvector[child] <- parents
  type <- match.arg(type)
  if (type == "parent")
    return(pvector[node])
  anc <- function(pvector, node) {
    res <- numeric(0)
    repeat {
      anc <- pvector[node]
      if (anc == 0)
        break
      res <- c(res, anc)
      node <- anc
    }
    res
  }
  if (!missing(node) && length(node) == 1)
    return(anc(pvector, node))
  else allAncestors(x)[node]
}

########## From phangorn treeManipulation.R
allAncestors <- function(x) {
  x <- reorder(x, "postorder")
  parents <- x$edge[, 1]
  child <- x$edge[, 2]
  l <- length(parents)
  res <- vector("list", max(x$edge))
  for (i in l:1) {
    pa <- parents[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

########## From phytools utilities.R
nodeHeights <- function (tree, ...) {
  if (hasArg(root.edge))
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge)
    ROOT <- if (!is.null(tree$root.edge))
      tree$root.edge
  else 0
  else ROOT <- 0
  nHeight <- function(tree) {
    tree <- reorder(tree)
    edge <- tree$edge
    el <- tree$edge.length
    res <- numeric(max(tree$edge))
    for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i,1]] + el[i]
    res
  }
  nh <- nHeight(tree)
  return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

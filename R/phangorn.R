################################################################################
#
# Functions in this file are from phangorn (treeManipulation.R) which is released under GPL >= 2.
#
################################################################################
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

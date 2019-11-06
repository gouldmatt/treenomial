
# treenomial

[![CRAN
status](https://www.r-pkg.org/badges/version/treenomial)](https://cran.r-project.org/package=treenomial)

## Overview

The package **treenomial** is an application of polynomials that
uniquely describe trees. It provides tools for tree analysis and
comparison based on polynomials. The core functions are:

  - **`treeToPoly()`**: convert trees to tree distinguishing polynomials
    described with coefficient matrices

  - **`polyToDistMat()`**: construct a distance matrix from multiple
    coefficient matrices using a distance measure

For the mathematical description of the tree defining polynomial see:

[Liu, Pengyu. “A tree distinguishing polynomial.” arXiv preprint
arXiv:1904.03332 (2019).](https://arxiv.org/abs/1904.03332)

## Installation

    install_github("mattgou1d/treenomial")

## Example tree and polynomial

Consider a three tip tree:

``` r
library(ape)
library(treenomial)

threeTipTree <- rtree(3, rooted = T)
plot.phylo(threeTipTree, use.edge.length = F, show.tip.label = F, direction = "downwards")
```

![](man/figures/README-threeTipTree-1.png)<!-- -->

It’s polynomial is x^3+xy+y which can equivalently be described with a
coefficient matrix where the element in the ith row, jth column
represents the y^(i-1) \* x^(j-1) coefficient:

``` r
treeToPoly(threeTipTree)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    1
#> [2,]    1    1    0    0
```

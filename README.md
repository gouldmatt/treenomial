# treenomial

### Install instructions (for now) 

Generate a personal access token from here using repo scope: 
https://github.com/settings/tokens

Copy the token to the clipboard 

Call this function in R replacing “foo” with the token 

```
devtools::install_github(repo = "mattgou1d/treenomial", ref = "master", auth_token = “foo” )
```

### Functions overview

* __`coefficientMatrix`__: convert trees to tree distinguishing polynomials described with coefficient matrices
* __`coefficientDist`__: construct a distance matrix from multiple coefficient matrices
* __`coefficientMds`__: construct 2d/3d mds plots based on tree/coefficient matrix data 
* __`wedge`__: perform the wedge operation on two coefficient matrices 
* __`allTrees`__: find all possible coefficient matrices up to a certain number of tips 

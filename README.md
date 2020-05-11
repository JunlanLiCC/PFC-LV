# PFC-LV
Principal fitted components method with a latent variable for binary responses (PFC-LV)
# Description
Using a continuous latent variable to represent an unobserved response underlying the binary response, a joint model is proposed for dimension reduction in binary regression. The minimal sufficient linear reduction is obtained, and an efficient expectation maximization algorithm is developed for carrying out maximum likelihood estimation.
# Usage
```R
EMquad_esti_para(X, Y, d, r, B=50, iter.max =30, err.tol = 0.0001)
```
> Arguments
* `X`: design matrix.
* `Y`: binary response vector.
* `d`: maximum number of directions to be used in estimating the reduction subspace.
* `r`: degree of polynomial basis function.
* `B`: number of Gauss-Hermite nodes
* `iter.max`: maximum number of iterations within GH-EM algorithm. The default is 30.
* `err.tol`: error threshold to stop the algorithm. The default is 0.0001.

# Examples


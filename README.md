# PFC-LV
Principal fitted components method with a latent variable for binary responses
# Description
We proposed a joint model for dimension reduction in binary regression using a continuous latent variable to represent an unobserved response underlying the binary response. The minimal sufficient linear reduction is obtained, and an efficient expectation maximization algorithm is developed for carrying out maximum likelihood estimation. The algorithm approximates latent variable integration by using Gauss–Hermite quadrature. 

# Usage
```R
EMquad_esti_para(X, Y, d, r, B=50, iter.max =30, err.tol = 0.0001)
```
#### Arguments
* `X`: design matrix
* `Y`: binary response vector
* `d`: maximum number of directions to be used in estimating the reduction subspace
* `r`: degree of polynomial basis function
* `B`: number of Gauss-Hermite nodes
* `iter.max`: maximum number of iterations within GH-EM algorithm (the default is 30)
* `err.tol`: error threshold to stop the algorithm (the default is 0.0001)

#### Value
* `MU`: estimate of μ
* `GAMMA`: estimate of Γ
* `BETA`: estimate of β
* `DELTA`: estimate of Δ
* `L`: approximate value of the log-likelihood
* `numpar`: number of parameters
* `bic`: value of Bayesian information criterion

# Example
```R
## load the required packages ##
library(dr) 
library(gam)
library(ldr)
library(splines)
library(mvtnorm)
library(logitnorm)
library(fastGHQuad)
library(randomForest)

## load the required functions ##
source("PFC_LV_function.R")

set.seed(2020)
B = 50
n = 400
d = 2
p = 10
r = 3
GAMMA <- diag(1,p,d)
DELTA <- diag(1,p,p)
MU <- as.matrix(rep(0,p),dim=c(p,1))
BETA <- matrix(c(1,1,0,1,-1,0),ncol=r,byrow = FALSE)
ZERO <- rep(0,p)

R.true <- solve(DELTA) %*% GAMMA
inv1 <- solve(crossprod(R.true,R.true))
P1 <- R.true %*% inv1 %*% t(R.true)

## 200 data replications ##
C = 200
err_PFC_O <- err_PFC_LV <- rep(0,C)
for(k in 1:C){
  ## generate data from the joint model ##
  a_true<-1
  b_true<-0
  theta_prior <- rnorm(n,0,1)
  Pi_theta <- exp(a_true*theta_prior+b_true)/(1+exp(a_true*theta_prior+b_true))
  X<-matrix(rep(0,n*p),ncol=p)
  Y<-rep(0,n)
  for(i in 1:n){
    X[i,] <- mvrnorm(1,MU + GAMMA %*% BETA %*% (fpi(theta_prior[i],r)),DELTA)
    Y[i] <- rbinom(1,1,Pi_theta[i])
  }

  ## gold standard ##
  fit_gold <- pfc(X, theta_prior, fy = bf(theta_prior,case="poly",degree=r), numdir=d, structure="unstr")
  R.pfc_gold <- solve(fit_gold$Deltahat)%*%fit_gold$Gammahat
  inv_pfc_gold <- solve(crossprod(R.pfc_gold,R.pfc_gold))
  P_gold <- R.pfc_gold %*% inv_pfc_gold %*% t(R.pfc_gold)

  ## calculate the estimation error of PFC-O ##
  err_PFC_O[k] <- sqrt(sum(diag(t(P1-P_gold)%*%(P1-P_gold))))

  ## apply PFC-LV ##
  GHEM_para <- EMquad_esti_para(X, Y, d=d, r=r, B=B, iter.max =30, err.tol = 0.0001)
  R.new <- solve(GHEM_para$DELTA) %*% GHEM_para$GAMMA
  inv2 <- solve(crossprod(R.new,R.new))
  P2<-R.new %*% inv2 %*% t(R.new)
  
  ## calculate the estimation error of PFC-LV ##
  err_PFC_LV[k] <- sqrt(sum(diag(t(P1-P2)%*%(P1-P2))))
}

```

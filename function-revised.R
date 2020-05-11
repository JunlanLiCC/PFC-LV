### PWSVM ###
lwpsvm <- function(x, y, H, lambda) 
{
  require(kernlab)
  n <- length(y)
  p <- ncol(x)
  
  step <- 1/H
  pi.grid <- seq(step, 1-step, by = step)
  
  bar.x <- apply(x, 2, mean)
  cov.x <- cov(x)
  
  # standardization
  temp <- eigen(cov.x)
  D <- diag(sqrt(temp$values))
  V <- temp$vectors
  sd.x <-  V %*% D
  inv.sd.x <- diag(1/sqrt(temp$values)) %*% t(V)
  z <- t(inv.sd.x %*% (t(x) - bar.x))
  
  w <- matrix(0, p, H)
  for (h in 1:(H-1)) {
    alpha <- rep(0, n)
    temp <- ksvm(x = z, y = as.factor(y), type = "C-svc", kernel = "vanilladot", kpar = list(), C = lambda, class.weights = c("1" = 1-pi.grid[h], "-1" = pi.grid[h]))
    alpha[temp@SVindex] <- unlist(temp@alpha)
    w[,h] <- apply(unlist(temp@coef) * z[temp@SVindex,], 2, sum)
  }
  
  psi <- solve(t(sd.x)) %*% w
  Mn <- matrix(0, p, p)
  for (h in 1:H) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
  
  obj <- eigen(Mn)
  obj
}

sqrt_inv <- function(A) 
{
  A <- (A + t(A)) / 2
  eig <- eigen(A)
  eig.vec <- eig$vec
  eig.val <- eig$val
  
  temp <- t(eig.vec) / sqrt(eig$val)
  return(eig.vec %*% temp)
}
### dr ###
dr_dirreg<-function(y, X, d, nslices)
{
  n <- dim(X)[1]
  p <- dim(X)[2] 
  
  #sort y and X
  y.org <- y; X.org <- X
  ord <- order(y.org)
  y <- y.org[ord]; X <- X.org[ord, ]
  
  X <- scale(X, center = TRUE, scale = F)
  sqrtinv.cov <- sqrt_inv(cov(X))
  Z <- X %*% sqrtinv.cov
  
  #equiparition slicing
  slicing <- dr.slices(y, nslices)
  nslices <- slicing$nslices
  slice.sizes <- slicing$slice.sizes	
  slice.indicator <- slicing$slice.indicator
  
  Z.means <- matrix(NA, p, nslices)
  M1 <- matrix(0, p, p)
  for (k in 1 : nslices) 
  {
    Z.k <- Z[slice.indicator == k, , drop = F]
    Z.means[, k] <- colMeans(Z.k)
    M10 <-  t(Z.k) %*% Z.k / slice.sizes[k] - diag(p)
    M1 <- M1 + M10 %*% t(M10) * slice.sizes[k] / n
  }
  M20 <- Z.means %*% (t(Z.means) * slice.sizes / n)
  M2 <- M20 %*% t(M20)
  M3 <- sum(diag(M20)) * M20
  
  M <- M1 + M2 + M3
  B.Z <- eigen(M)$vec[, 1 : d, drop = F]
  
  B.X <- sqrtinv.cov %*% B.Z
  return(list(B.X = B.X, B.Z = B.Z,slice.sizes=slice.sizes))
}

### the polynomial basis function ###
fpi <- function(pii,r=3){
  result<-c()
  for(i in 1:r){
    result[i] <- pii^i
  }
  dim(result) <- c(r,1)
  return(result)
}

### GH-EM ### 
EMquad_esti_para <- function(X, Y, d=2, r=r, B=50, iter.max =30, err.tol = 0.001){
  ###initialize the parameters
  rf.fit<-randomForest(X,as.factor(Y),mtry=2,ntree=500)
  rf.probs <- predict(rf.fit, type = "prob")[,2]
  rf.probs[which(rf.probs>0.975)] <- 0.975
  rf.probs[which(rf.probs<0.025)] <- 0.025
  ynew<-logit(rf.probs)
  fit1 <- pfc(X,ynew, fy = bf(ynew,case="poly",degree=r), numdir=d, structure="unstr")
  MU <- fit1$Muhat - fit1$Gammahat %*% fit1$Betahat %*% attr(bf(ynew,case="poly",degree=r),"scaled:center")
  BETA<-fit1$Betahat
  GAMMA<-fit1$Gammahat
  DELTA<-fit1$Deltahat
  a<-1
  b<-logit(sum(Y)/n)
  
  n <- dim(X)[1]
  p <- dim(X)[2] 
 
  nodes <- gaussHermiteData(B)$x
  quad_weights <- gaussHermiteData(B)$w
  all_nodes <- t(matrix(rep(nodes,n),nrow = B,ncol = n))
  theta <- sqrt(2)*all_nodes
  f_theta <- array(rep(0), c(r,n,B))
  for(i in 1:n){
    for(j in 1:B){
      f_theta[,i,j] <- fpi(theta[i,j],r)
    }
  }

  #iteration between estimation and maximization
  iter <- 1; err2 <- 1; L<-1
  while (err2 >= err.tol & iter <= iter.max){
    #print(iter)
    
    omega <- rep(0,n)
    u <- matrix(rep(0,n*B), nrow = n, ncol = B)
    for(i in 1:n){
      for(j in 1:B){
        u[i,j] <- (1/sqrt(pi))*dmvnorm(X[i,],mean = MU + GAMMA %*% BETA %*% f_theta[ ,i,j],sigma=DELTA,log=FALSE) *
          dbinom(Y[i],1,1/(1+exp(-a*sqrt(2)*nodes[j]-b)),log=FALSE)
        omega[i] <- omega[i] + u[i,j]*quad_weights[j]
      } 
    }
    
    ######## coefficient estimator
    I <- rep(1,B)
    Y_tilde <- kronecker(Y,I)
    I_w <- rep(1,n)
    quad_weights_new<-as.matrix(quad_weights,nrow=1)
    w<- kronecker(quad_weights_new,I_w)
    w_new<-array(as.vector(w),c(n,B))
    wei_matrix<-u*w_new/(omega*n)
    wei <- as.vector(t(wei_matrix))
    theta.new<-as.vector(t(theta))
    glm.fit<-glm(as.vector(Y_tilde) ~ theta.new , family = quasibinomial , weight = wei)
    b<-as.numeric(glm.fit$coefficients[1])
    a<-as.numeric(glm.fit$coefficients[2])
    est_para <- max_para_quad(X, Y, theta,u, omega, quad_weights, f_theta, n, B, d, p, a,b)
    
    ################################ 
    Lnew <-  est_para$L
    err1 <- abs(Lnew - L)
    err2 <- abs((Lnew - L)/L)
    
    MU <- est_para$Muhat
    GAMMA <- est_para$Gammahat
    BETA <- est_para$Betahat
    DELTA <- est_para$Deltahat
    
    iter <- iter + 1
    L<- Lnew
  }
  LV<- rep(0,n)
  fac <- matrix(rep(0,n*B), nrow = n, ncol = B)
  for(i in 1:n){
    for(j in 1:B){
      fac[i,j] <- (1/sqrt(pi))*dmvnorm(X[i,],mean = MU + GAMMA %*% BETA %*% f_theta[ ,i,j],sigma=DELTA,log=FALSE) *
        dbinom(Y[i],1,1/(1+exp(-a*sqrt(2)*nodes[j]-b)),log=FALSE) * quad_weights[j]
    } 
    LV[i] <-log(sum(fac[i,]))
  }
  L_value<-sum(LV)
  
  paranum <- 2+p*(p+3)/2+r*d+d*(p-d)
  BIC <- -2*L_value+log(n)*paranum
 
  return(list(MU = MU, GAMMA = GAMMA, BETA  = BETA, DELTA = DELTA, L=L_value, numpar=paranum,bic=BIC))
}

### M-step ### 
max_para_quad <- function(X, Y, theta, u, omega, quad_weights, f_theta, n, B, d, p,a,b)
{
  "%^%"<-function(M, pow) 
  { 
    if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
    eigenM = eigen(M) 
    return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))  
  }
  
  Xc <- scale(X, center=TRUE, scale=FALSE)
  I <- rep(1,B)
  X_tilde<- kronecker(I,Xc)
  
  f_sum <- rep(0,r)
  for(i in 1:n){
    for(j in 1:B){
      f_sum <- f_sum + f_theta[,i,j] * quad_weights[j] * u[i,j] / omega[i]
    }
  }
  mean_f <- f_sum/n
  
  f_center <- array(rep(0), c(r,n,B))
  f_theta_weigh <- array(rep(0), c(r,n,B))
  X_tilde_weigh <- matrix(rep(0),ncol = p,nrow = n*B)
  for(i in 1:n){
    for(j in 1:B){
      X_tilde_weigh[i+(j-1)*n,] <- sqrt(quad_weights[j] * B * u[i,j]/omega[i]) * X_tilde[i+(j-1)*n,]
      f_center[,i,j] <- f_theta[,i,j] - mean_f
      f_theta_weigh[,i,j] <- sqrt(quad_weights[j] * B * u[i,j]/omega[i]) * f_center[,i,j]
    }
  }
 
  fc<-c()

  fc=t(array(as.vector(f_theta_weigh),c(r,n*B))) 
  s <- svd(fc)
  U <- s$u
  D <- diag(s$d)
  V <- s$v
  P_FX <- U %*% t(U) %*% X_tilde_weigh
  Sigmahat <- cov(X_tilde_weigh)
  Sigmahat_fit <- cov(P_FX)
  Sigmahat_res <- Sigmahat - Sigmahat_fit
  
  sqrt_Sigmahat_res <- Sigmahat_res%^%0.5 
  Inv_Sqrt_Sigmahat_res <- solve(sqrt_Sigmahat_res)
  lf_matrix <- Inv_Sqrt_Sigmahat_res %*% Sigmahat_fit %*% Inv_Sqrt_Sigmahat_res
  V_evalues <- eigen(lf_matrix, symmetric=TRUE)$values
  evalues <- V_evalues[1:d]
  
  Vhat <- eigen(lf_matrix, symmetric=TRUE)$vectors
  Vhat_d <- matrix(Vhat[,1:d], ncol=d)
  Gammahat <- (Sigmahat_res%^%0.5)%*%Vhat_d%*%solve((t(Vhat_d)%*%Sigmahat_res%*%Vhat_d)%^%0.5)  
  
  Khat<-diag(0, p) 
  if (d < min(ncol(fc),p)) {diag(Khat)[(d+1):p]<- V_evalues[(d+1):p]}
  Deltahat <- sqrt_Sigmahat_res%*%Vhat%*%(diag(p)+Khat)%*%t(Vhat)%*%sqrt_Sigmahat_res
  Betahat <- ((t(Vhat_d)%*%Sigmahat_res%*%Vhat_d)%^%0.5)%*%t(Vhat_d)%*%solve(Sigmahat_res%^%0.5)%*%t(X_tilde_weigh)%*% U %*% (D%^%(-1)) %*% t(V)
  Muhat <- apply(X, 2, mean) - Gammahat %*% Betahat %*% mean_f
  pi_theta <- 1/(1+exp(-a*theta-b))
  temp1 <- -(n/2)*log(det(Sigmahat_res)) -(n*p/2)
  temp2 <- -(n/2)*sum(log(1 + V_evalues[(d+1):p]))
  sum_temp3 <- 0
  for(i in 1:n){
    for(j in 1:B){
      if(exp(a*theta[i,j]+b)==Inf||exp(a*theta[i,j]+b)==-Inf){
        sum_temp3<-sum_temp3 + (Y[i] * (a*theta[i,j]+b)-(a*theta[i,j]+b)) * u[i,j] * quad_weights[j]/omega[i]
      } else{
        sum_temp3<-sum_temp3 + (Y[i] * (a*theta[i,j]+b)-log(1+exp(a*theta[i,j]+b))) * u[i,j] * quad_weights[j]/omega[i]
      }
    }
  }
  temp3 <- (1/2)*sum_temp3
  loglik <- temp1 + temp2 + temp3
  return(list(Muhat=Muhat, Betahat=Betahat,Gammahat=Gammahat, Deltahat=Deltahat , L=loglik))
}


### nonlinear logit ###
max_para_quad_1 <- function(X, Y, theta, u, omega, quad_weights, f_theta, n, B, d, p,h_theta)
{
  "%^%"<-function(M, pow)
  {
    if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
    eigenM = eigen(M)
    return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))
  }
  
  Xc <- scale(X, center=TRUE, scale=FALSE)
  I <- rep(1,B)
  X_tilde<- kronecker(I,Xc)
  
  f_sum <- rep(0,r)
  for(i in 1:n){
    for(j in 1:B){
      f_sum <- f_sum + f_theta[,i,j] * quad_weights[j] * u[i,j] / omega[i]
    }
  }
  mean_f <- f_sum/n
  
  f_center <- array(rep(0), c(r,n,B))
  f_theta_weigh <- array(rep(0), c(r,n,B))
  X_tilde_weigh <- matrix(rep(0),ncol = p,nrow = n*B)
  for(i in 1:n){
    for(j in 1:B){
      X_tilde_weigh[i+(j-1)*n,] <- sqrt(quad_weights[j] * B * u[i,j]/omega[i]) * X_tilde[i+(j-1)*n,]
      f_center[,i,j] <- f_theta[,i,j] - mean_f
      f_theta_weigh[,i,j] <- sqrt(quad_weights[j] * B * u[i,j]/omega[i]) * f_center[,i,j]
    }
  }
  
  # fc is the weighted centered f.
  fc<-c()
  fc=t(array(as.vector(f_theta_weigh),c(r,n*B)))
  s <- svd(fc)
  U <- s$u
  D <- diag(s$d)
  V <- s$v
  P_FX <- U %*% t(U) %*% X_tilde_weigh
  Sigmahat <- cov(X_tilde_weigh)
  Sigmahat_fit <- cov(P_FX)
  Sigmahat_res <- Sigmahat - Sigmahat_fit
  
  sqrt_Sigmahat_res <- Sigmahat_res%^%0.5
  Inv_Sqrt_Sigmahat_res <- solve(sqrt_Sigmahat_res)
  lf_matrix <- Inv_Sqrt_Sigmahat_res %*% Sigmahat_fit %*% Inv_Sqrt_Sigmahat_res
  V_evalues <- eigen(lf_matrix, symmetric=TRUE)$values
  evalues <- V_evalues[1:d]
  
  Vhat <- eigen(lf_matrix, symmetric=TRUE)$vectors
  Vhat_d <- matrix(Vhat[,1:d], ncol=d)
  Gammahat <- (Sigmahat_res%^%0.5)%*%Vhat_d%*%solve((t(Vhat_d)%*%Sigmahat_res%*%Vhat_d)%^%0.5)
  
  Khat<-diag(0, p)
  if (d < min(ncol(fc),p)) {diag(Khat)[(d+1):p]<- V_evalues[(d+1):p]}
  Deltahat <- sqrt_Sigmahat_res%*%Vhat%*%(diag(p)+Khat)%*%t(Vhat)%*%sqrt_Sigmahat_res
  Betahat <- ((t(Vhat_d)%*%Sigmahat_res%*%Vhat_d)%^%0.5)%*%t(Vhat_d)%*%solve(Sigmahat_res%^%0.5)%*%t(X_tilde_weigh)%*% U %*% (D%^%(-1)) %*% t(V)
  Muhat <- apply(X, 2, mean) - Gammahat %*% Betahat %*% mean_f
  pi_theta <- exp(h_theta)/(1+exp(h_theta))
  temp1 <- -(n/2)*log(det(Sigmahat_res)) -(n*p/2)
  temp2 <- 0
  if (d < min(ncol(fc),p)) temp2 <- -(n/2)*sum(log(1 + V_evalues[(d+1):p]))
  sum_temp3 <- 0
  for(i in 1:n){
    for(j in 1:B){
      if(exp(h_theta[i,j])==Inf||exp(h_theta[i,j])==-Inf){
        sum_temp3<-sum_temp3 + (Y[i] * (h_theta[i,j])-(h_theta[i,j])) * u[i,j] * quad_weights[j]/omega[i]
      } else{
        sum_temp3<-sum_temp3 + (Y[i] * (h_theta[i,j])-log(1+exp(h_theta[i,j]))) * u[i,j] * quad_weights[j]/omega[i]
      }
    }
  }
  temp3 <- (1/2)*sum_temp3
  loglik <- temp1 + temp2 +temp3
  return(list(Muhat=Muhat, Betahat=Betahat,Gammahat=Gammahat, Deltahat=Deltahat , L=loglik))
}

EMquad_esti_para_1 <- function(X, Y, d=2,r=r,B=50,iter.max =35, err.tol = 0.001){
  
  ###initialize the parameters
  rf.fit<-randomForest(X,as.factor(Y),mtry=2,ntree=500)
  rf.probs <- predict(rf.fit, type = "prob")[,2]
  rf.probs[which(rf.probs>0.975)] <- 0.975
  rf.probs[which(rf.probs<0.025)] <- 0.025
  ynew<-logit(rf.probs)
  fit1 <- pfc(X,ynew, fy = bf(ynew,case="poly",degree=r), numdir=d, structure="unstr")
  MU <- fit1$Muhat - fit1$Gammahat %*% fit1$Betahat %*% attr(bf(ynew,case="poly",degree=r),"scaled:center")
  BETA<-fit1$Betahat
  GAMMA<-fit1$Gammahat
  DELTA<-fit1$Deltahat
  a<-1
  b<-logit(sum(Y)/n)
  h_theta <- sqrt(2)*all_nodes
  
  n <- dim(X)[1]
  p <- dim(X)[2] 
  
  nodes <- gaussHermiteData(B)$x
  quad_weights <- gaussHermiteData(B)$w
  all_nodes <- t(matrix(rep(nodes,n),nrow = B,ncol = n))
  theta <- sqrt(2)*all_nodes
  f_theta <- array(rep(0), c(r,n,B))
  for(i in 1:n){
    for(j in 1:B){
      f_theta[,i,j] <- fpi(theta[i,j],r)
    }
  }
  
  #iteration between estimation and maximization
  iter <- 1; err2 <- 1; L<-1
  while (err2 >= err.tol & iter <= iter.max){
    print(iter)
    
    omega <- rep(0,n)
    u <- matrix(rep(0,n*B), nrow = n, ncol = B)
    for(i in 1:n){
      for(j in 1:B){
        u[i,j] <- (1/sqrt(pi))*dmvnorm(X[i,],mean = MU + GAMMA %*% BETA %*% f_theta[ ,i,j],sigma=DELTA,log=FALSE) *
          dbinom(Y[i],1,1/(1+exp(-h_theta[i,j])),log=FALSE)
        omega[i] <- omega[i] + u[i,j]*quad_weights[j]
      }
    }
    
    ######## coefficient estimator
    I <- rep(1,B)
    Y_tilde <- kronecker(Y,I)
    I_w <- rep(1,n)
    quad_weights_new<-as.matrix(quad_weights,nrow=1)
    w<- kronecker(quad_weights_new,I_w)
    w_new<-array(as.vector(w),c(n,B))
    wei_matrix<-u*w_new/(omega*n)
    wei <- as.vector(t(wei_matrix))
    theta.new<-as.vector(t(theta))
    gam.fit<-gam(as.vector(Y_tilde) ~ s(theta.new,4), family = quasibinomial, weight = wei)
    h_theta <- matrix(predict(gam.fit, newdata=data.frame(theta.new)),byrow=TRUE,ncol=B)
    est_para <- max_para_quad_1(X, Y, theta,u, omega, quad_weights, f_theta, n, B, d, p, Y_tilde,h_theta)
    
    ################################
    #print(est_para)
    Lnew <-  est_para$L
    err1 <- abs(Lnew - L)
    err2 <- abs((Lnew - L)/L)
    
    MU <- est_para$Muhat
    GAMMA <- est_para$Gammahat
    BETA <- est_para$Betahat
    DELTA <- est_para$Deltahat
    
    iter <- iter + 1
    L<- Lnew
  }
  LV<- rep(0,n)
  fac <- matrix(rep(0,n*B), nrow = n, ncol = B)
  for(i in 1:n){
    for(j in 1:B){
      fac[i,j] <- (1/sqrt(pi))*dmvnorm(X[i,],mean = MU + GAMMA %*% BETA %*% f_theta[ ,i,j],sigma=DELTA,log=FALSE) *
        dbinom(Y[i],1,1/(1+exp(-h_theta[i,j])),log=FALSE) * quad_weights[j]
    }
    LV[i] <-log(sum(fac[i,]))
  }
  L_value<-sum(LV)
  
  paranum <- p*(p+3)/2+r*d+d*(p-d)
  BIC <- -2*L_value+log(n)*paranum
  
  return(list(MU = MU, GAMMA = GAMMA, BETA  = BETA, DELTA = DELTA,L=L_value, numpar=paranum,bic=BIC))
}


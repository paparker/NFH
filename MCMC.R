FH_Fit <- function(Y, X, S2, sig2b=1000, iter=1000, burn=500){
  n <- length(Y)
  p <- ncol(X)
  r <- n
  tau2 <- 1
  eta1 <- rnorm(n)
  tau2Out <- rep(NA, iter)
  beta1Out <- matrix(NA, nrow=iter, ncol=p)
  eta1Out <- matrix(NA, nrow=iter, ncol=r)
  XP <- cbind(X, diag(n))
  thetaOut <- xbOut <- matrix(NA, nrow=n, ncol=iter)
  ll <- rep(NA, iter)
  XtX <- t(X)%*%Diagonal(x=1/S2)%*%X
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## Fixed effects
    meanBeta <- t(X)%*%Diagonal(x=1/S2)%*%(Y-eta1)
    Ub <- chol(forceSymmetric(XtX + 1/sig2b*Diagonal(p)))
    b <- rnorm(p)
    beta1 <- beta1Out[i,] <- backsolve(Ub, backsolve(Ub, meanBeta, transpose=T) + b)
    
    
    ## Random effects
    var <- solve(Diagonal(x=1/S2+(1/tau2)))
    meanEta <- var%*%Diagonal(x=1/S2)%*%(Y-X%*%beta1)
    eta1 <- eta1Out[i,] <- rnorm(n, mean=as.numeric(meanEta), sd=sqrt(diag(var)))
    
    
    theta <- thetaOut[,i] <- X%*%beta1 + eta1
    xbOut[,i] <- theta - eta1
    
    tau2 <- tau2Out[i] <- 1/rgamma(1,
                                   0.1 + n/2,
                                   0.1 + t(theta-X%*%beta1)%*%(theta-X%*%beta1)/2)
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(S2), log=T))
    
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(S2), log=T))
  return(list(Preds=thetaOut[,-c(1:burn)], Tau2=tau2Out[-c(1:burn)], Beta=beta1Out[-c(1:burn),], XB=xbOut[,-c(1:burn)], DIC=DIC))
}



FH_RNN_Fit <- function(Y, X, S2, nh=200, iter=1000, burn=500){
  n <- length(Y)
  p <- ncol(X)
  r <- n
  tau2 <- sig2b <- 1
  eta1 <- rnorm(n)
  tau2Out <- sig2bOut <- rep(NA, iter)
  beta1Out <- matrix(NA, nrow=iter, ncol=nh)
  eta1Out <- matrix(NA, nrow=iter, ncol=r)
  HL <- plogis(X%*%matrix(rnorm(p*nh, sd=1), nrow=p))
  XtX <- t(HL)%*%Diagonal(x=1/S2)%*%HL
  thetaOut <- xbOut <- matrix(NA, nrow=n, ncol=iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## Fixed effects
    meanBeta <- t(HL)%*%Diagonal(x=1/S2)%*%(Y-eta1)
    Ub <- chol(forceSymmetric(XtX + 1/sig2b*Diagonal(nh)))
    b <- rnorm(nh)
    beta1 <- beta1Out[i,] <- backsolve(Ub, backsolve(Ub, meanBeta, transpose=T) + b)

    ## Random effects
    var <- solve(Diagonal(x=1/S2+(1/tau2)))
    meanEta <- var%*%Diagonal(x=1/S2)%*%(Y-HL%*%beta1)
    eta1 <- eta1Out[i,] <- rnorm(n, mean=as.numeric(meanEta), sd=sqrt(diag(var)))

    
    theta <- thetaOut[,i] <- HL%*%beta1 + eta1
    xbOut[,i] <- theta - eta1
    
    tau2 <- tau2Out[i] <- 1/rgamma(1,
                                   0.1 + n/2,
                                   0.1 + t(theta-HL%*%beta1)%*%(theta-HL%*%beta1)/2)
    
    sig2b <- sig2bOut[i] <- 1/rgamma(1,
                                   20 + nh/2,
                                   8 + sum(beta1^2)/2)
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(S2), log=T))
    
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(S2), log=T))
  return(list(Preds=thetaOut[,-c(1:burn)], Tau2=tau2Out[-c(1:burn)], sig2b=sig2bOut[-c(1:burn)], Beta=beta1Out[-c(1:burn),], XB=xbOut[,-c(1:burn)], DIC=DIC))
}

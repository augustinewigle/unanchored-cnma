# Arm-based Anchored model - full multivariate, no decomposition - slower
model{
  
  for(i in 1:nstudy) {
    
    eps[1,i] <- 0 # no heterogeneity in arm 1 because it contains reference treatment
    
    eps[2:(narm[i]),i] ~ dmnorm.vcov(rep(0,(narm[i]-1)), sigma^2*Sigma[1:(narm[i]-1), 1:(narm[i]-1),i]) # remaining arms are MVN
    
    for(j in 1:narm[i]) {
      
      logitp[i,j] <- alpha[i] + V1[j,1:(ncomp-1),i]%*%d1 + eps[j,i]
      p[i,j] <- ilogit(logitp[i,j])
      r[i,j] ~ dbinom(p[i,j], n[i,j])
      
      # deviance contribution
      #rhat[i,j] <- p[i,j]*n[i,j]
      #dev[i,j] <- 2*(r[i,j]*(log(r[i,j])-log(rhat[i,j])) + (n[i,j]-r[i,j])*(log(n[i,j]-r[i,j])-log(n[i,j]-rhat[i,j])))
      
    }
    
    #dev2[i] <- sum(dev[i,1:narm[i]])
    
  }
  #totresdev <- sum(dev2[])
  
  # priors
  for(i in 1:nstudy) {
    
    alpha[i] ~ dnorm(0, 0.001)
    
  }
  for(j in 1:(ncomp-1)) {
    
    d1[j] ~ dnorm(0, 0.001)#I(-5,)
    
  }
  sigma ~ dunif(0,2)
  
}
# Unanchored contrast-level model
model {
  
  for(i in 1:nstudy) {
    
    ystar[1:(narm[i]-1),i] ~ dmnorm.vcov(mu[1:(narm[i]-1),i], Sstar[1:(narm[i]-1), 1:(narm[i]-1),i])
    
    mu[1:(narm[i]-1),i] <- Ustar[1:(narm[i]-1),1:narm[i],i]%*%V[1:narm[i],1:ncomp,i]%*%d + eps[1:(narm[i]-1),i]
    
    eps[1:(narm[i]-1),i] ~ dmnorm.vcov(rep(0, narm[i]-1), sigma^2*Sigma[1:(narm[i]-1), 1:(narm[i]-1), i])
    
  }
  
  # Priors
  for(j in 1:ncomp) {
    
    d[j] ~ dnorm(0, 0.001)
    
  }
  sigma ~ dunif(0,2)
  
}
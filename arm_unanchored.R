# Arm-based Unanchored model

model {
  
  for(i in 1:nstudy) {
    
    logitp[i,1:narm[i]] <- alpha[i] + V[1:narm[i],1:ncomp,i]%*%d + eps[1:narm[i],i]
    
    eps[1:narm[i],i] ~ dmnorm.vcov(rep(0,narm[i]), sigma^2*Sigma[1:narm[i], 1:narm[i],i])
    
    for(j in 1:narm[i]) {
      
      p[i,j] <- ilogit(logitp[i,j])
      r[i,j] ~ dbinom(p[i,j], n[i,j])
      
    }
    
  }
  
  # priors
  for(i in 1:nstudy) {
    
    alpha[i] ~ dnorm(0, 0.001)
    
  }
  for(j in 1:ncomp) {
    
    d[j] ~ dnorm(0, 0.001)
    
  }
  sigma ~ dunif(0,2)
  
}
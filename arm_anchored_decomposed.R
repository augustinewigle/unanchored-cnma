# Arm-based Anchored model - decomposed
model{
  
  for(i in 1:nstudy) {
    
    for(j in 1:narm[i]) {
      
      logitp[i,j] <- alpha[i] + V1[j,1:(ncomp-1),i]%*%d1 + eps[i,j] 
      p[i,j] <- ilogit(logitp[i,j])
      r[i,j] ~ dbinom(p[i,j], n[i,j])
      
    }
    
    eps[i,1] <- 0 # no heterogeneity in arm 1 because arm 1 contains reference treatment
    md[i,1] <- 0
    taud[i,1] <- tau
    
    for(j in 2:narm[i]) {
      
      eps[i,j] ~ dnorm(md[i,j], taud[i,j])
      md[i,j] <- sum(eps[i,1:(j-1)])/(j-1) 
      taud[i,j]<- tau*(2*(j-1)/(j)) 
      
    }
    
  }
  
  tau <- 1/sigma^2
  
  # priors
  for(i in 1:nstudy) {
    
    alpha[i] ~ dnorm(0, 0.001)
    
  }
  for(j in 1:(ncomp-1)) {
    
    d1[j] ~ dnorm(0, 0.001)
    
  }
  sigma ~ dunif(0,2)
  
}
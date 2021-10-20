# Unanchored contrast-level model
model {
  
  for(i in 1:nstudy) {
    
    ystar[1:(narm[i]-1),i] ~ dmnorm.vcov(mu[1:(narm[i]-1),i], Sstar[1:(narm[i]-1), 1:(narm[i]-1),i])
    
    for(j in 1:(narm[i]-1)) {
      
      mu[j,i] <- Ustar[j,1:narm[i],i]%*%V[1:narm[i],1:ncomp,i]%*%d + eps[j,i]
      
      eps[j,i] ~ dnorm(0+md[i,j], taud[i,j])
      
      taud[i,j] <- 2*(j)/(sigma^2*(j+1))
      
    }
    
    md[i,1] <- 0
    
    for(j in 2:(narm[i]-1)) {
      
      md[i,j] <- sum(eps[1:(j-1), i])/(j)
      
    }
    
  }
  
  # Priors
  for(j in 1:ncomp) {
    
    d[j] ~ dnorm(0, 0.001)
    
  }
  sigma ~ dunif(0,2)
  
}
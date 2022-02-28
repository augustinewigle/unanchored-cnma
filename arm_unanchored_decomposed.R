# Arm - based unanchored model

model {
  
  for(i in 1:nstudy) {
    
    for(j in 1:narm[i]) {
      
      logitp[i,j] <- alpha[i] + V[j,1:ncomp,i]%*%d + eps[i,j]
      eps[i,j] ~ dnorm(0 + md[i,j], taud[i,j])
      p[i,j] <- ilogit(logitp[i,j])
      r[i,j] ~ dbinom(p[i,j], n[i,j])
      taud[i,j] <- 2*j/(sigma^2*(j+1))
      
    }
    
    md[i,1] <- 0
    
    for(j in 2:narm[i]) {
      
      md[i,j] <- sum(eps[i,1:(j-1)])/j
      
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
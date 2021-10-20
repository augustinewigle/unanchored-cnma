#' Prepare data for running through JAGS
#' @param data_obj the dataframe
#' @param studycol name of the column with study ids
#' @param components vector of names of components (must have columns in dataframe). Order of components in this vector determines
#' order in the code (the first treatment will be treatment 1, etc.)
#' @return A list of lists.
#' arm_un - a list of the data to be passed to JAGS if the unanchored arm-based model is used
#' contrast_un - a list of the data to be passed to JAGS if the unanchored contrast-based model is used
#' jagsdata - the data matrices created from the input data object - mostly for debugging

make_jags_data <- function(data_obj, studycol, components) {
  
  data_obj$study <- data_obj[,studycol]
  
  nstudy <- length(unique(data_obj$study))
  ncomp <- length(components)
  
  # Temp variable which makes it easier to create the matrices in Jags data format
  temp <- data_obj %>% arrange(study) %>% group_by(study) %>% mutate(arm = row_number())%>%
    ungroup() %>% gather("variable", "value", -study, -arm)  %>% spread(arm, value)
  
  jagsdata <- list()
  
  for(v in unique(temp$variable)) {
    
    jagsdata[[v]] <- as.matrix(temp %>% filter(variable == v) %>% select(-study, -variable))
    
  }
  
  narm <- numeric(nstudy)
  maxna <- ncol(jagsdata[[1]])
  
  for(i in 1:nstudy) {
    
    narm[i] <- maxna - sum(is.na(jagsdata[[1]][i,]))
    
  }
  
  # Make Sigma

  Sigma <- array(dim = c(maxna, maxna, nstudy))
  for(i in 1:nstudy) {

    for(j in 1:narm[i]) {

      for(k in 1:narm[i]) {

        Sigma[j,k,i] <- ifelse(j == k, 1, 0.5)

      }

    }

  }
  
  # Make V and V1
  V <- array(dim = c(maxna, ncomp, nstudy))
  V1 <- array(dim = c(maxna, ncomp-1, nstudy))
  for(i in 1:nstudy) {
    
    for(j in 1:narm[i]) {
      
      for(k in 1:ncomp) {
        
        V[j,k,i] <- ifelse(jagsdata[[components[k]]][i,j]>1,1,jagsdata[[components[k]]][i,j]) # make sure its all just 1 and zeros
        
      }
      
      V1[j,1:(ncomp-1),i] <- V[j,-1,i]
      
    }
    
  }
  
  # Make Ustar
  
  Ustar <- array(dim = c(maxna-1, maxna, nstudy))
  
  for(i in 1:nstudy) {
    
    for(k in 1:(narm[i]-1)) {# comparison {
      
      for(j in 1:narm[i]) {
        
        if(j == 1) {
          
          Ustar[k,j,i] <- -1
          
        } else if (k == j-1) {
          
          Ustar[k,j,i] <- 1
          
        } else {
          
          Ustar[k,j,i] <- 0
          
        }
        
      }
      
    }
  }
  
  # Make ystar and Sstar
  
  ystar <- matrix(nrow = maxna-1, ncol = nstudy)
  Sstar <- array(dim = c(maxna-1, maxna-1, nstudy))
  
  for(i in 1:nstudy) {
    
    for(j in 1:(maxna-1)) {
      
      newr <- ifelse(jagsdata$r[i,j+1] == 0, jagsdata$r[i,j+1] + 0.5, jagsdata$r[i,j+1])
      newr1 <- ifelse(jagsdata$r[i,1] == 0, jagsdata$r[i,1] + 0.5, jagsdata$r[i,1])
      newn <- ifelse(jagsdata$r[i,j+1] == 0, jagsdata$n[i,j+1] + 0.5, jagsdata$n[i,j+1])
      newn1 <- ifelse(jagsdata$r[i,1] == 0, jagsdata$n[i,1] + 0.5, jagsdata$n[i,1])
      
      ystar[j,i] <- log((newr/(newn-newr))/(newr1/(newn1-newr1)))
      
      # For Sstar
      for(k in 1:(narm[i]-1)) {
        
        if(j == k) {
          
          Sstar[j,k,i] <- 1/newr + 1/(newn-newr) + 1/newr1 +  1/(newn1-newr1)
          
        } else {
          
          var_arm1 <- 1/newr1*(1-newr1/newn1)
          Sstar[j,k,i] <- var_arm1
          Sstar[k,j,i] <- Sstar[j,k,i]
          
        }
        
      }
      
    }
    
  }
  
  
  return(list(arm_un = list(narm = narm,
                            nstudy = nstudy,
                            ncomp = ncomp,
                            V = V,
                            Sigma = Sigma,
                            n = jagsdata$n,
                            r = jagsdata$r),
              contrast_un = list(narm = narm,
                                 nstudy = nstudy,
                                 ncomp = ncomp,
                                 V = V,
                                 Sigma = Sigma,
                                 Ustar = Ustar,
                                 ystar = ystar,
                                 Sstar = Sstar),
              jagsdata = jagsdata))
  
  
}
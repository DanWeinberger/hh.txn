longiniSim.m = function( max.hh.size, N.households, RR.community, RR.sar, baseline.community.prob, baseline.HH.prob,sizes){

  K = max.hh.size  # Maximum household size = # of columns 
  m = array(0,dim=c(K+1,K,2))
  logit_P1 <- rep(NA,2)
  logit_P2 <- rep(NA,2)
  P1 <- rep(NA,2)
  P2 <- rep(NA,2)
  B <- rep(NA,2)
  Q <- rep(NA,2)
  
  beta1=log(set.RR.community)
  delta1=log(set.RR.sar)
  
  for(v in 1:2){ #vaccine gro\up
    
    for(k in 1:max.hh.size){
      ds[1:(k+1),k,v] ~ dmulti(m[1:(k+1), k,v], n.hh[k,v])
    }
    
    #Take baseline probs and convert to logit
    alpha1 <- log(baseline.community.prob/(1-baseline.community.prob))
    alpha2 <- log(baseline.HH.prob/(1-baseline.HH.prob))
    
    
    logit_P1[v] <- alpha1 + (v==2)*beta1
    logit_P2[v] <- alpha2 + (v==2)*delta1

    P1[v] <- exp(logit_P1[v])/(1+ exp(logit_P1[v]))
    P2[v] <- exp(logit_P2[v])/(1+ exp(logit_P2[v]))
    
    #Could define B=1-P1; Q=1-P2; then have P1 and P2 be a function of covariates e.g., proportion vaccinated)
    B[v] <- 1 - P1[v]
    Q[v] <- 1 - P2[v]
    
    m[1,1,v] = B[v]   # Probability of 0 of 1 HH member infected 
    m[2,1,v] = 1-B[v] # Probability of 1 of 1 HH member infected 
    
    
    for (k in 2:max.hh.size){
      m[1,k,v] = B[v]^k # Probability everyone in HH escapes infection from the community 
      
      for (j in 1:(k-1)){
        m[j+1,k,v] = choose(k,j)*m[j+1,j, v]*(B[v]^(k-j))*Q[v]^(j*(k-j)) # Probabilty j out of k HH members infected
      }
      
      m[k+1,k,v] = 1-sum(m[1:k,k,v]) #Probability that everyone in HH infected
    }
  }
  
  sizes <- rep(50, max.hh.size)
  simN <- array(NA, dim=dim(m))
    for(k in 1: max.hh.size){
      for(v in 1:2){
        simN[,k,v] <- rmultinom(n=1, size=sizes[k], prob=m[,k,v])
      }
    }

  return(simN)
}

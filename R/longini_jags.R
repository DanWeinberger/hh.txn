longini_jags  <- "
model{
  
  for(k in 1:max.hh.size){
    ds[,k] ~ dmulti(m[1:(k+1), k], n.hh[k])
  }
  
  logit_B ~ dnorm(0, 1e-4)
  logit_Q ~ dnorm(0, 1e-4)

  B <- exp(logit_B)/(1+ exp(logit_B))
  Q <- exp(logit_Q)/(1+ exp(logit_Q))

  
  m[1,1] = B   # Probability of 0 of 1 HH member infected 
  m[2,1] = 1-B # Probability of 1 of 1 HH member infected 
  
  for (k in 2:max.hh.size){
    m[1,k] = B^k # Probability everyone in HH escapes infection from the community 
    for (j in 1:k){
      
      m[j+1,k] = choose_kj_mat[k,j]*m[j+1,j]*(B^(k-j))*Q^(j*(k-j)) # Probabilty j out of k HH members infected
    }
      #m[k+1:max.hh.size,k] <- 0
     m[k+1,k] = 1-sum(m[,k])
  }
}
"


longiniLL = function(pars,data){
  B = pars[1]
  Q = pars[2]
  
  K = dim(data)[2]   # Maximum household size = # of columns 
  m = matrix(0,K+1,K)
  
  m[1,1] = B   # Probability of 0 of 1 HH member infected 
  m[2,1] = 1-B # Probability of 1 of 1 HH member infected 
  
  for (k in 2:K){
    m[1,k] = B^k # Probability everyone in HH escapes infection from the community 
    for (j in 1:k){
      m[j+1,k] = choose(k,j)*m[j+1,j]*(B^(k-j))*Q^(j*(k-j)) # Probabilty j out of k HH members infected
    }
    m[k+1,k] = 1-sum(m[,k])
  }
  
  llikl = 0
  for (k in 1:K){
    for (j in 0:k){
       llikl = llikl + data[j+1,k]*log(m[j+1,k])
     }
   }
  
 # llikl = sum(data*log(m)) #this doesn't work bc only non 0s in ~upper half
  return(-llikl)
}
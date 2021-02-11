#Pitzer AJE 2007  Winbugs code: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7110150/
  model_string<-"
model{
{
  p[1] <- 0
  tp[1] <- log(1 - p[1])          #transform of p
  
  for( d in 2 : 14 ) {              #d=day of illness
    p[d] ~ dbeta(1,1)I(0,1)  #likelihood of txn on each day of illness d
    tp[d] <- log(1 - p[d])    #prob of escaping on day d
  }
  
  for( i in 1 : 2013 ) { #i=2013 contacts 
    nY[i] ~ dbern(q[i]) 
  }
  
  for( i in 1 : 2013 ) {   #i=2013 contacts 
    for( t in 1 : 38 ) {    #t= 38 ___ ???? time from sympton onset to hospitalization for index
      for( j in 1 : 5 ) {    #Index case j--there are multiple index cases
        te[i , t , j] <- tp[1 + T[i , t , j]] #Ti,j(t) is the day of illness of index case j at time t
      }
      loge[i , t] <- sum(te[i , t , ]) #The probability that subject i escapes infection at time t over all of j
    }
    for( t in 2 : 38 ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1]
    }
    ploge[i , 1] <- 0
  }
  
  ##120 infected contacts, 1,893 uninfected contacts
  
  for( i in 1 : 120 ) {
    for( t in 1 : 38 ) { #up to 38 day-why 38? Paper says 19 days of follow up
      inf[i , t] <- (exp(ploge[i , t])) * (1 - exp(loge[i , t]))  #Equation 3 parts?
      g[i , t] <- v[40 - t] #g is the probability density function for the incubation period
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ]) #Infected contacts  #inprod=sum(x1*x2) Equation 3
  }
  
  for( i in 121 : 2013 ) {
    q[i] <- exp(ploge[i , 38]) #Uninfected contacts; equation 2
  }
  for( t in 1 : 38 ) {
    incb[t] <- v[40 - t]
  }
  a ~ dgamma(.001,.001) #part of Incubation period?
  b ~ dgamma(.001,.001) #part of Incubation period?
  
  #v is the density function for the incubation period
  for( t in 2 : 42 ) {
    v[t] <- exp((a - 1) * log(b * (t - 1.5)) - b * (t - 1.5) - loggam(a)) * b
  }
}

"
  
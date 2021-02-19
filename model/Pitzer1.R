#Pitzer AJE 2007  Winbugs code: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7110150/
  model_string<-"
model{

  p[1] <- 0
  tp[1] <- log(1 - p[1])          #transform of p
  
  for( d in 2 : 14 ) {              #d=day of illness
    p[d] ~ dbeta(1,1)I(0,1)  #prior on txn on each day of illness d
    tp[d] <- log(1 - p[d])    #prob of escaping on day d
  }
  
  for( i in 1 : N ) { #i=N contacts 
    nY[i] ~ dbern(q[i]) 
  }
  
  for( i in 1 : N ) {   #i=N contacts 
    for( t in 1 : n.times ) {    #t= n.times ___ ???? time from sympton onset to hospitalization for index
      for( j in 1 : max.index ) {    #Index case j--there are multiple index cases
        te[i , t , j] <- tp[1 + T[i , t , j]] #Ti,j(t) is the day of illness of index case j at time t
      }
      loge[i , t] <- sum(te[i , t , ]) #The probability that subject i escapes infection at time t over all of j
    }
    for( t in 2 : n.times ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1] #cumulative over all prev times
    }
    ploge[i , 1] <- 0
  }
  
  ##N.infected infected contacts, 1,893 uninfected contacts
  
  for( i in 1 : N.infected ) {
    for( t in 1 : n.times ) { #up to n.times day-why n.times? Paper says 19 days of follow up
      inf[i , t] <- (exp(ploge[i , t])) * (1 - exp(loge[i , t]))  #Equation 3 parts?
      g[i , t] <- v[40 - t] #g is the probability density function for the incubation period
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ]) #Infected contacts  #inprod=sum(x1*x2) Equation 3
  }
  
  for( i in (N.infected+1) : N ) {
    q[i] <- exp(ploge[i , n.times]) #Uninfected contacts; equation 2
  }
  for( t in 1 : n.times ) {
    incb[t] <- v[40 - t]
  }
  a ~ dgamma(.001,.001) #part of Incubation period?
  b ~ dgamma(.001,.001) #part of Incubation period?
  
  #v is the density function for the incubation period
  for( t in 2 : 42 ) {
    v[t] <- exp((a - 1) * log(b * (t - 1.5)) - b * (t - 1.5) - loggam(a)) * b
  }
    v[1] <- 0


   for( i in 1 : 81 ) {
    Q[i] <- sum(expos[i , ]) #Likelihood of all possible incubation periods for 81 contacts with known exposure period
  }
  
  for( i in 1 : 81 ) {
    obs[i] ~ dbern(Q[i]) #obs is the infection status of cases with known incubation period (=1)
  }

  for( i in 1 : 81 ) {
    for( j in 1 : 41 ) {
      expos[i , j] <- v[incub[i , j] + 1] / sum(v[]) #Likelihood incubation period for i is incub[i,j] days
    }
  }

}

"
  
library(rjags)

model_string <- "
model{
  tp0[1] <- 0
  tp1[1] <- 0
  tpu[1] <- 0
  p[1] ~ dbeta(1,1)I(0,1)
  tp0[2] <- -p[1]
  tp1[2] <- -p[1]*exp(S[1])
  tpu[2] <- -p[1]*exp(S[2])
  
  for( d in 1 : 4 ) {
    P[d] ~ dnorm(0.0,1.0E-6)
  }
  for( d in 1 : 2 ) {
    S[d] ~ dnorm(0.0,1.0E-6)
  }
  for( d in 1 : 3 ) { 
    p[d+1] <- p[1]*exp(P[d])
    tp0[d+2] <- -p[d+1]
    tp1[d+2] <- -p[d+1]*exp(S[1])
    tpu[d+2] <- -p[d+1]*exp(S[2])
  }
  for( d in 4 : 7 ) {
    p[d+1] <- p[d-3]*exp(P[4])
    tp0[d+2] <- -p[d+1]
    tp1[d+2] <- -p[d+1]*exp(S[1])
    tpu[d+2] <- -p[d+1]*exp(S[2])
  }
  
  for( i in 1 : 2014 ) {
    nY[i] ~ dbern(q[i])
  }
  
  for( i in 1 : 934 ) {
    for( t in 1 : 38 ) {
      for( j in 1 : 5 ) {
        te[i , t , j] <- tp0[1+T[i , t , j]]
      }
      loge[i , t] <- sum(te[i , t , ]) 
    }
    for( t in 2 : 38 ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1]
    }
    ploge[i , 1] <- 0
  }
  for( i in 66 : 934 ) {
    q[i] <- exp(ploge[i , 38])
  }
  for( i in 1 : 65 ) {
    for( t in 1 : 38 ) {
      g[i,t] <- incb[t]/inprod(exps[i , ],incb[ ])
      exps[i,t] <- step(sum(T[i , t , ]) - 1)
      inf[i , t] <- exp(ploge[i , t]) * (1 - exp(loge[i , t]))
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ])
  }
  
  for( i in 935 : 1850 ) {
    for( t in 1 : 38 ) {
      for( j in 1 : 5 ) {
        te[i , t , j] <- tp1[1+T[i , t , j]]
      }
      loge[i , t] <- sum(te[i , t , ]) 
    }
    for( t in 2 : 38 ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1]
    }
    ploge[i , 1] <- 0
  }
  for( i in 980 : 1850 ) {
    q[i] <- exp(ploge[i , 38])
  }
  for( i in 935 : 979 ) {
    for( t in 1 : 38 ) {
      g[i,t] <- incb[t]/inprod(exps[i , ],incb[ ])
      exps[i,t] <- step(sum(T[i , t , ]) - 1)
      inf[i , t] <- exp(ploge[i , t]) * (1 - exp(loge[i , t]))
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ])
  }
  
  for( i in 1851 : 1859 ) {
    for( t in 1 : 38 ) {
      for( j in 1 : 5 ) {
        te[i , t , j] <- tpu[1+T[i , t , j]]
      }
      loge[i , t] <- sum(te[i , t , ]) 
    }
    for( t in 2 : 38 ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1]
    }
    ploge[i , 1] <- 0
  }
  for( i in 1854 : 1859 ) {
    q[i] <- exp(ploge[i , 38])
  }
  for( i in 1851 : 1853 ) {
    for( t in 1 : 38 ) {
      g[i,t] <- incb[t]/inprod(exps[i , ],incb[ ])
      exps[i,t] <- step(sum(T[i , t , ]) - 1)
      inf[i , t] <- exp(ploge[i , t]) * (1 - exp(loge[i , t]))
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ])
  }
  
  for( i in 1860 : 2014 ) {
    for( t in 1 : 38 ) {
      for( j in 1 : 5 ) {
        te[i , t , j] <- tpu[1+T[i , t , j]]
      }
      loge[i , t] <- sum(te[i , t , ]) 
    }
    for( t in 2 : 38 ) {
      ploge[i , t] <- ploge[i , t - 1] + loge[i , t-1]
    }
    ploge[i , 1] <- 0
  }
  for( i in 1868 : 2014 ) {
    q[i] <- exp(ploge[i , 38])
  }
  for( i in 1860 : 1867 ) {
    for( t in 1 : 38 ) {
      g[i,t] <- incb[t]/inprod(exps[i , ],incb[ ])
      exps[i,t] <- step(sum(T[i , t , ]) - 1)
      inf[i , t] <- exp(ploge[i , t]) * (1 - exp(loge[i , t]))
    }
    q[i] <- 1 - inprod(g[i , ],inf[i , ])
  }
  
  for( t in 1 : 38 ) {
    incb[t] <- v[40 - t]
  }
  a ~ dgamma(.001,.001)
  b ~ dgamma(.001,.001)
  for( t in 2 : 42 ) {
    v[t] <- exp((a - 1) * log(b * (t - 1.5)) - b * (t - 1.5) - loggam(a)) * b
  }
  for( i in 1 : 81 ) {
    Q[i] <- sum(expos[i , ])
  }
  for( i in 1 : 81 ) {
    obs[i] ~ dbern(Q[i])
  }
  for( i in 1 : 81 ) {
    for( j in 1 : 41 ) {
      expos[i , j] <- v[incub[i , j] + 1] / sum(v[])
    }
  }
  v[1] <- 0
}
"

##############################################################
#Model Fitting
##############################################################
inits1=list(".RNG.seed"=c(123), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(456), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(789), ".RNG.name"='base::Wichmann-Hill')

#data.sim <- gen_sim_data()


T <- data.list$T
nY <- data.list$nY
obs <- data.list$obs
incub <- data.list$incub

##############################################
#Model Organization
##############################################
model_spec<-textConnection(model_string)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('T'=T,
                                 'nY'=nY,
                                 'obs'=obs,
                                 'incub'=incub),
                       n.adapt=10000, 
                       n.chains=3)











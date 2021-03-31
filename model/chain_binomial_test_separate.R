chain_binomial_jags <- "
model{

### Part 1: d_jk(t) definition -- this first part is computed out of JAGS
selecter ~ runif(1,10000) #select random draw 
d[,,,] <- d.iter[round(selecter),,,,]
#### Part 2: Likelihood definition

for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){ 

    y[i,j] ~ dbern(q[i,j]) #(Person j in HH i: 1= uninfected; 0= infected)

    for(t in 1:time.study.HH[i]){
        log.p[i,j,(N.hh.members[i]+1),t] <- delta0 + beta1*vaxdose1[i,j] +  beta2*vaxdose2[i,j] #outside HH infection
      for(k in 1:N.hh.members[i]){ #k is different from j
        log.prob_uninf[i,j,k,t] <- step(abs(j-k)-0.5)*log((1-p[i,j,k,t])^(d[i,j,k,t]))
        log.p[i,j,k,t] <- alpha0 + beta1*vaxdose1[i,j] +  beta2*vaxdose2[i,j] + kappa1*vaxdose1[i,k] +  kappa2*vaxdose2[i,k]
      }

      for(k in 1:(N.hh.members[i] + 1)){
        p[i,j,k,] <- exp(log.p[i,j,k,])
        log.prob_inf_timej[i,j,k] <- step(abs(j-k)-0.5)*(1-log((1-p[i,j,k,day.exposed[i,j]])^(d[i,j,k,day.exposed[i,j]]))
      }
    }
    
    prob_uninf_to_timej[i,j] <- exp(sum(log.prob_uninf[i,j,,1:(day.exposed[i,j]-1)]))
    prob_uninf[i,j] <- exp(sum(log.prob_uninf[i,j,,1:time.study.HH[i]]))
    prob_inf_timej[i,j] <- exp(sum(log.prob_inf_timej[i,j,]))
    prob_inf[i,j] <- prob_uninf_to_timej[i,j]* prob_inf_timej[i,j]
    q[i,j] <- prob_uninf[i,j]*prob_inf[i,j]
  }  
}

delta0 ~dnorm(0,1e-4) 
alpha0 ~dnorm(0,1e-4)

beta1 ~ dnorm(0,1e-4) 
beta2 ~ dnorm(0,1e-4) 
kappa1 ~ dnorm(0,1e-4)
kappa2 ~ dnorm(0,1e-4)

}
}
"


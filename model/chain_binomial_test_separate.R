chain_binomial_jags <- "
model{

#select random draw from the delay distributions
selecter ~ dunif(1,10000) 
selecter1 <- round(selecter)

#### Part 2: Likelihood definition

for(i in 1:N.HH){ 
k.last[i] <- N.hh.members[i] + 1
  for(j in 1:N.hh.members[i]){ 

    #y[i,j] ~ dbern(q[i,j]) #(Person j in HH i: 1= uninfected; 0= infected)

    for(t in 1:time.study.HH[i, selecter1]){
        log.p[i,j,(N.hh.members[i]+1),t] <- delta0 + beta1*vaxdose1[i,j]  #outside HH infection
      for(k in 1:N.hh.members[i]){ #k is different from j
        d[i,j,k,t] <- step(abs(j-k)-0.5)*step( t > day.infectious[i,k,selecter1]+0.5)*step( t < min(day.exposed[i,j,selecter1], day.infectious.end[i,k,selecter1]) + 0.5)
        log.prob_uninf[i,j,k,t] <- step(abs(j-k)-0.5)*log((1-p[i,j,k,t])^(d[i,j,k,t]))
        log.p[i,j,k,t] <- alpha0 + beta1*vaxdose1[i,j]  + kappa1*vaxdose1[i,k] 
      }
       
      d[i,j,k.last[i],t] <- step( t > day.infectious[i,k.last[i],selecter1]+0.5)*step( t < min(day.exposed[i,j,selecter1], day.infectious.end[i,k.last[i],selecter1]) + 0.5)

      for(k in 1:(N.hh.members[i] + 1)){
        p[i,j,k,] <- exp(log.p[i,j,k,])
        log.prob_inf_timej[i,j,k] <- step(abs(j-k)-0.5)*log(1-(1-p[i,j,k,day.exposed[i,j,selecter1]])^(d[i,j,k,day.exposed[i,j,selecter1]]))
      }
    }
    
    prob_uninf_to_timej[i,j] <- exp(sum(log.prob_uninf[i,j,,1:(day.exposed[i,j]-1)]))
    prob_uninf[i,j] <- exp(sum(log.prob_uninf[i,j,,1:time.study.HH[i,selecter1]]))
    prob_inf_timej[i,j] <- exp(sum(log.prob_inf_timej[i,j,]))
    prob_inf[i,j] <- prob_uninf_to_timej[i,j]* prob_inf_timej[i,j]
    like[i,j] <- y[i,j]*prob_uninf[i,j] + (1-y[i,j])*prob_inf[i,j]
  }  
}

delta0 ~dnorm(0,1e-4) 
alpha0 ~dnorm(0,1e-4)

beta1 ~ dnorm(0,1e-4) 
beta2 ~ dnorm(0,1e-4) 
kappa1 ~ dnorm(0,1e-4)
kappa2 ~ dnorm(0,1e-4)

}
"


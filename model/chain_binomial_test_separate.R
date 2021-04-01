chain_binomial_jags <- "
model{

#select random draw from the delay distributions
#selecter ~ dunif(1,10000) 
#selecter1 ~ dround(selecter)

#selecter1 <- 1

#### Part 2: Likelihood definition

for(i in 1:N.HH){ 
#k.last[i] <- N.hh.members[i] + 1
  for(j in 1:N.hh.members[i]){ 

    y[i,j] ~ dbern(like[i,j]) #(Person j in HH i: 1= uninfected; 0= infected)

    for(t in 1:time.study.HH[i, selecter1]){
        log.p[i,j,(N.hh.members[i]+1),t] <- delta0 + beta1*vaxdose1[i,j]  #outside HH infection
      for(k in 1:N.hh.members[i]){ #k is different from j
        d[i,j,k,t] <- step(abs(j-k)-0.5)*step( t > day.infectious[i,k,selecter1]+0.5)*step( t < min(day.exposed[i,j,selecter1], day.infectious.end[i,k,selecter1]) + 0.5)
        log.prob_uninf[i,j,k,t] <- step(abs(j-k)-0.5)*log((1-p[i,j,k,t])^(d[i,j,k,t]))
        log.p[i,j,k,t] <- alpha0 + beta1*vaxdose1[i,j]  + kappa1*vaxdose1[i,k] 
        p[i,j,k,t] <- exp(log.p[i,j,k,t])

      }
       
      d[i,j,(N.hh.members[i]+1),t] <- 1 #step( t > day.infectious[i,(N.hh.members[i] + 1),selecter1]+0.5)*step( t < min(day.exposed[i,j,selecter1], day.infectious.end[i,(N.hh.members[i] + 1),selecter1]) + 0.5)
      p[i,j,(N.hh.members[i]+1),t] <- exp(log.p[i,j,(N.hh.members[i]+1),t])
      log.prob_uninf[i,j,(N.hh.members[i]+1),t] <- log((1-p[i,j,(N.hh.members[i]+1),t])^(d[i,j,(N.hh.members[i]+1),t]))


    }
    
    for(k in 1:(N.hh.members[i])){
        log.prob_inf_timej[i,j,k] <- step(abs(j-k)-0.5)*log(1-(1-p[i,j,k,day.exposed[i,j,selecter1]])^(d[i,j,k,day.exposed[i,j,selecter1]]))
    }
    
    log.prob_inf_timej[i,j,N.hh.members[i]+1] <- log(1-(1-p[i,j,(N.hh.members[i]+1),day.exposed[i,j,selecter1]])^(d[i,j,(N.hh.members[i]+1),day.exposed[i,j,selecter1]]))
    prob_uninf_to_timej[i,j] <- exp(sum(log.prob_uninf[i,j,1:N.hh.members[i],1:(day.exposed[i,j,selecter1]-1)]))
    prob_uninf[i,j] <- exp(sum(log.prob_uninf[i,j,1:N.hh.members[i],1:time.study.HH[i,selecter1]]))
    prob_inf_timej[i,j] <- exp(sum(log.prob_inf_timej[i,j,1:(N.hh.members[i]+1)]))
    prob_inf[i,j] <- prob_uninf_to_timej[i,j]* prob_inf_timej[i,j]
    like[i,j] <- y2[i,j]*prob_uninf[i,j] + (1-y2[i,j])*prob_inf[i,j]+ 1e-6  
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


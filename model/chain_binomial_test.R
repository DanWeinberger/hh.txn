chain_binomial_jags <- "
model{

### Part 1: d_jk(t) definition -- this first part can be computed out of JAGS

for(i in 1:N.HH){
  for(j in 1:N.hh.members[i]){
    #day.matrix=day of test for the person

    day.infectious[i,j] <- day.matrix[i,j] - infect.dist[i,j]  #infectious prior to test (t_onset in the document)
    day.exposed[i,j] <- day.infectious[i,j] - expose.dist[i,j]  #latent period (t_infect in the document)
    day.infectious.end[i,j] <- day.infectious[i,j] + end.inf[i,j]  #how long infectious after onset of infectiousness

    infect.dist[i,j] ~ dgamma(sh1, ra1) #duration infectiousness
    expose.dist[i,j] ~ dgamma(sh2, ra2)
    end.inf[i,j] ~ dgamma(sh3,ra3)

    for(t in 1:time.study.HH[i]){
       d[i,j,(N.hh.members[i]+1),t] <- 1 #placeholder for exogenous

      for(k in 1:N.hh.members[i]){
        #### step(abs(j-k)-0.5) ensures don't count the current j in k

        d[i,j,k,t] <- step(abs(j-k)-0.5)*step( t > day.infectious[i,k]+0.5)*step( t < min(day.exposed[i,j], day.infectious.end[i,k]) + 0.5)
     
       }
    }
  }
}


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


# Hyperpriors for the latent distributions
# parameterized by mode (m) and standard deviation (sd):
sh1 <- 1 + m1 * ra1
ra1 <- ( m1 + sqrt( m1^2 + 4*sd1^2 ) ) / ( 2 * sd1^2 )
m1 ~ dunif(2,6)
sd1 ~ dunif(2,3) 

sh2 <- 1 + m2 * ra2
ra2 <- ( m2 + sqrt( m2^2 + 4*sd2^2 ) ) / ( 2 * sd2^2 )
m2 ~ dunif(2,6)
sd2 ~ dunif(2,3) 

ra3 <- ( m3 + sqrt( m3^2 + 4*sd3^2 ) ) / ( 2 * sd3^2 )
sh3 <- 1 + m3 * ra3
m3 ~ dunif(2,6) #days
sd3 ~ dunif(2,3) #SD on days

#mu1 ~dnorm(0,1e-4)
#mu2 ~ dnorm(0,1e-4)
#tau1 ~ dgamma(0.001,0.001)
#tau2 ~ dgamma(0.001,0.001)


delta0 ~dnorm(0,1e-4) 
alpha0 ~dnorm(0,1e-4)

beta1 ~ dnorm(0,1e-4) 
beta2 ~ dnorm(0,1e-4) 
kappa1 ~ dnorm(0,1e-4)
kappa2 ~ dnorm(0,1e-4)

}
}
"


zelner_jags2 <- "
model{
for(i in 1:N.HH){

    y[i] ~ dbinom(prob.hh[i],N.hh.members[i]) #y=N infected in HH i
    
    S0[i] = S_sum[i,1] # Initial number of susceptibles in each household
    Sf[i] = S_sum[i,tmax[i]] # Final number of susceptibles in each household
    ninfect[i] = S0[i] - Sf[i] # number infected in each household 

    for(t in 1:tmax[i]){
      log_prob_no_inf_t[i,t] <- -1* S_sum[i,t]*(beta * I_sum[i,t] + alpha)
    }
    
    for(j in 1:N.hh.members[i]){
       prob_inf[i,j] <- infected_matrix[i,j] * S_sum[i,day.matrix[i,j]]*(beta * I_sum[i,day.matrix[i,j]] + alpha) #for infected people only
    }
    
    prob_no_inf[i] <- exp(sum(log_prob_no_inf_t[i,1:tmax[i]])) #P no infections over all time intervals
    prob_inf.hh[i] <- exp(sum(log(prob_inf[i,])))
   
    prob.hh[i] <- prob_inf.hh[i] * prob_no_inf[i] 
}

##### Attaching latent piece

for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){
    
    day.infectious[i,j] <- day.matrix[i,j] - infect.dist[i,j]  #infectious prior to test
    day.exposed[i,j] <- day.infectious[i,j] - expose.dist[i,j]  #latent period
    day.infectious.end[i,j] <- day.infectious[i,j] + end.inf[i,j]  #how long infectious after test
    
    infect.dist[i,j] ~ dgamma(sh1, ra1) #duration infectiusness
    expose.dist[i,j] ~dgamma(sh2, ra2)
    end.inf[i,j] ~dgamma(sh3,ra3)
    
    for(t in 1:tmax[i]){ #how can we do this without looping over all t?
      #need to add 1 to some of these?
      I[i,j,t] <- step(t - day.infectious[i,j]) * (1-step(t-day.infectious.end[i,j])) * infected_matrix[i,j] #they are infected and during infectious period
      E[i,j,t] <- step(t - day.exposed[i,j]) * (1-step(t-day.infectious[i,j])*infected_matrix[i,j])
      S[i,j,t] <- step(-0.5 + (1-infected_matrix[i,j]) + (infected_matrix[i,j]*(1-step(t-day.exposed[i,j])) )) #if they do not have recorded infection OR if they haven't yet been exposed by time t
    }
  }
}

for(i in 1:N.HH){
  for(t in 1:tmax[i]){ 
    I_sum[i,t] <- sum(I[i,,t])
    E_sum[i,t] <- sum(E[i,,t])
    S_sum[i,t] <- sum(S[i,,t])
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




alpha ~ dnorm(0,1e-4)
beta ~ dnorm(0,1e-4)

}
"

# for(i in 1:N.hh){
#     S[i,1] <- N.hh.members[i] - 1
#     E[i,1] <- 0
#     I[i,1] <- 1
#     R[i,1] <- 0
#   
#     for(t in 2:42){
#       S[i,t] <- -S[i,(t-1)]*(beta * I[i,(t-1)] + alpha)
#       E[i,t] <- E[i, (t-1)]+ S[i,(t-1)]*(beta * I[i,(t-1)] + alpha) - E[i, (t-1)]*delta
#       I[i,t] <- I[i,(t-1)] + S[i,(t-1)]*(beta * I[i,(t-1)] + alpha) - I[i,(t-1)] * epsilon
#       R[i,t] <- I[i,(t-1)] * epsilon
#   }
# }

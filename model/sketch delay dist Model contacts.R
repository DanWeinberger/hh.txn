zelner_jags2 <- "
model{

for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){
    
    y[i,j] ~ dbern(q[i,j]) #y=N infected in HH i

    #day.matrix=day of test for the person
    day.infectious[i,j] <- day.matrix[i,j] - infect.dist[i,j]  #infectious prior to test
    day.exposed[i,j] <- day.infectious[i,j] - expose.dist[i,j]  #latent period
    day.infectious.end[i,j] <- day.infectious[i,j] + end.inf[i,j]  #how long infectious after test
    
    index.case[i,j] <- (1-y[i,j]) * step((min(day.exposed[i,]) - day.exposed[i,j] + 0.5) #indicator for whether this person is an index

    infect.dist[i,j] ~ dgamma(sh1, ra1) #duration infectiousness
    expose.dist[i,j] ~dgamma(sh2, ra2)
    end.inf[i,j] ~dgamma(sh3,ra3)
    
    alpha[i,j] <- delta[1] + delta[2]*vax[i,j]  #could add in things like HH-level random intercept, or intercept to link matched households
    beta[i,j] <- epsilon[1] + epsilon[2]*vax[i,j]
    
   ##NOTE THESE ARE PROBABLY NOT RIGHT--PROB NEED TO FLIP SOME OF THEM AROUND
    prob_no_inf_uninf[i,j] <- exp(sum(log_prob_no_inf_t[i,1:day.exposed[i,j]])) #P no infections over all time intervals
    
    prob_no_inf_inf_person[i,j] <- exp(sum(log_prob_no_inf_t[i,1:(day.exposed[i,j]-1) )) ])) #P no infections prior to infection
    
    prob_inf[i,j] <- infected_matrix[i,j] * (1-exp(log_prob_no_inf[i,j,day.exposed[i,j]])) #prob infected on dayfor infected people only
    
    q[i,j] <- (prob_no_inf_inf_person[i,j]) * (1 - prob_inf[i,j]) * infected_matrix[i,j] + #prob for infection at time t and not before
      (prob_no_inf_uninf[i,j]) * (1- infected_matrix[i,j] ) + 1e-6 #prob for uninfected peopel    
    
    #for(t in 1:tmax[i]){ #how can we do this without looping over all t?
     # I[i,j,t] <- step(t - day.infectious[i,j]) * (1-step(t-day.infectious.end[i,j])) * infected_matrix[i,j] #they are infected and during infectious period
    #  log_prob_no_inf_t[i,j,t] <-  log(1 - beta[i,j] *(sum_I[i,(t-1)])  + alpha[i,j]) )

    #}
  
  }
}
##Loop over the HH index cases
for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){ #j is the contact
    for(m in 1:N.hh.members[i]){ #m is the index
      
      dur.inf[i,j,m] <- (1-y[i,j])*(day.infectious.end[i,m] - day.infectious[i,m]) #0 for uninfected people
      
      dur.inf.contact[i,j,m] <- step(day.exposed[i,m] - day.infectious[i,j] ) * (
                                (1-y[i,j])* step(day.infectious.end[i,m] - day.infectious[i,j] +0.5) *(day.infectious[i,m] - day.infectious[i,j]) + # if contact is infected before end of index infectious period 
                                (1-y[i,j])* step(day.infectious[i,j] - day.infectious.end[i,m]) *(day.infectious[i,m] - day.infectious.end[i,m])  # if contact is infected AFTER end of index infectious period 
                        )
      #step(abs(j-m)-0.5) ensures don't count the current j in m
    log.prob.uninf.contact[i,j,m]= step(abs(j-m)-0.5) *log(
          y[i,j]* (1-p_inf[i,j,m])^dur.inf[i,j,m] + #For uninfected person j
          (1-y[i,j])* index.case[i,j] * 0  + #For INDEX case in the household --just community rate--add outside loop
          (1-y[i,j])*(1 - index.case[i,j])) * (1-p_inf[i,j,m])^(dur.inf.contact[i,j,m])   #For non-index case in the household
 
    log.prob.uninf.contact.day.inf[i,j,m] =step(abs(j-m)-0.5) * log(  (
       (1-y[i,j]) * (1-p_inf[i,j,m])^(step(day.exposed[i,j] -  day.infectious[i,m] )*(step(day.infectious.end[i,m] - day.exposed[i,j] ) ) )
        )
      )
    }
    
    prob.uninf[i,j] <- exp(sum(log.prob.uninf.contact[i,j,]))
    prob.inf.day.inf[i,j] <- (1-sum(log.prob.uninf.contact.day.inf[i,j,m]))
    q[i,j] <- (1-index.case[i,j] * (1-y[i,j]) * (1 - alpha[i,j]) * prob.uninf[i,j] * prob.inf[i,j] + #for contacts
              index.case[i,j] * alpha[i,j] + #for index case
              y[i,j] *(1 - alpha[i,j]) * prob.uninf[i,j]                           #for uninfected person
  }
}

#How many Infected people are there in the HH at each time point?
for(i in 1:N.HH){
  for(t in 1:tmax[i]){ 
    I_sum[i,t] <- sum(I[i,,t])
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

for(k in 1:2){
delta[k] ~dnorm(0,1e-4)
epsilon[k] ~dnorm(0,1e-4)
}


}
"


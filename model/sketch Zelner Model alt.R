#Note this doesn't make sense because we aren't taking into account observed time of PCR
zelner_jags2 <- "
model{
for(i in 1:N.HH){

   for(t in 1:tmax[i]){
      log_prob_no_inf_t[i,t] <- -1 * S_sum[i,t]*(beta * I_sum[i,t] + alpha)
   }
    
     for(j in 1:N.hh.members[i]){
      y[i,j] ~ dbern(total_prob[i,j]) #y=infected matrix
     
     ##NOTE THESE ARE PROBABLY NOT RIGHT--PROB NEED TO FLIP SOME OF THEM AROUND
      prob_no_inf[i,j] <- exp(sum(log_prob_no_inf_t[i,1:day.matrix[i,j]])) #P no infections over all time intervals
      prob_inf[i,j] <- infected_matrix[i,j] * S_sum[i,day.matrix[i,j]]*(beta * I_sum[i,day.matrix[i,j]] + alpha) #for infected people only
      
      total_prob[i,j] <- (1 - prob_no_inf[i,j]) * prob_inf[i,j]
    }
    

}

##### Attaching latent piece

for(i in 1:N.hh){
    Dur.Inf[i] ~ dgamma(sh1, ra1) #duration infectiusness
    Dur.Latent[i] ~dgamma(sh2, ra2) #Duration latent

    delta[i] <- 1/Dur.Inf[i]
    epsilon[i] <- 1/Dur.Latent[i]
    
    S[i,1] <- N.hh.members[i] - 1
    E[i,1] <- 0
    I[i,1] <- 1
    R[i,1] <- 0
    NewI[i,1] <- 0
    CumNewI[i,1] <- 0

    for(t in 2:42){
      dS[i,t] <- -S[i,(t-1)]*(beta * I[i,(t-1)] + alpha)
      dE[i,t] <- S[i,(t-1)]*(beta * I[i,(t-1)] + alpha) - E[i, (t-1)]*delta[i]
      dI[i,t] <- S[i,(t-1)]*(beta * I[i,(t-1)] + alpha) - I[i,(t-1)] * epsilon[i]
      dR[i,t] <- I[i,(t-1)] * epsilon
      
      S[i,t] <- dS[i,t] + S[i,(t-1)]
      E[i,t] <- dE[i,t] + E[i,(t-1)]
      I[i,t] <- dI[i,t] + I[i,(t-1)]
      R[i,t] <- dR[i,t] + R[i,(t-1)]
      NewI[i,t] <- S[i,(t-1)]*(beta * I[i,(t-1)] + alpha)
      CumNewI[i,t] <- NewI[i,(t-1)] + NewI[i,t]
    }
    SecondaryI[i] <- sum(NewI[i,])
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



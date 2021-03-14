#Note this doesn't make sense because we aren't taking into account observed time of PCR
zelner_jags2 <- "
model{
for(i in 1:N.HH){
     for(j in 1:N.hh.members[i]){
     
        for(t in 1:tmax[i]){
          log_prob_no_inf_t[i,j,t] <- log( 1 - (beta[i,j] * sum_I[i,t] + alpha[i,j]) )
        }
        
        y[i,j] ~ dbern(q[i,j]) #y=UNinfected matrix
     
     ##NOTE THESE ARE PROBABLY NOT RIGHT--PROB NEED TO FLIP SOME OF THEM AROUND
      prob_no_inf_uninf[i,j] <- exp(sum(log_prob_no_inf_t[i,j,1:day.matrix[i,j]])) #P no infections over all time intervals
      
      prob_no_inf_inf_person[i,j] <- exp(sum(log_prob_no_inf_t[i,j,1:((day.matrix[i,j] - 1)) ])) #P no infections prior to infection
      
      prob_inf[i,j] <- infected_matrix[i,j] * (beta[i,j] * sum_I[i,day.matrix[i,j]] + alpha[i,j]) #for infected people only

    q[i,j] <- (1 - prob_no_inf_inf_person[i,j]) * (1 - prob_inf[i,j]) * infected_matrix[i,j] + #prob for infection at time t and not before
                    (1 - prob_no_inf_uninf[i,j]) * (1- infected_matrix[i,j] ) + 1e-6 #prob for uninfected peopel    
}
    

}

##### Attaching latent piece

for(i in 1:N.HH){
     for(j in 1:N.hh.members[i]){

 
     alpha[i,j] <- exp(log.alpha[i,j])
     beta[i,j]  <- exp(log.beta[i,j])
     log.alpha[i,j] ~ dnorm(mu1,tau1)
     log.beta[i,j] ~ dnorm(mu2,tau2)
     
    Dur.Inf[i,j] ~ dgamma(sh1, ra1) #duration infectiusness
    Dur.Latent[i,j] ~dgamma(sh2, ra2) #Duration latent

    delta[i,j] <- 1/Dur.Inf[i,j]
    epsilon[i,j] <- 1/Dur.Latent[i,j]
    
    ##In this formulation, S,E,I,R are defined by individual; sum_S, sum_E, sum_I are defined by household
    S[i,j,1] <- (N.hh.members[i] - 1)/N.hh.members[i]
    E[i,j,1] <- 0
    I[i,j,1] <- 1/N.hh.members[i]
    R[i,j,1] <- 0
    NewI[i,j,1] <- 0
    CumNewI[i,j,1] <- 0
    dS[i,j,1] <-0
    dE[i,j,1] <-0
    dI[i,j,1] <-0
    dR[i,j,1] <-0

    for(t in 2:43){
      dS[i,j,t] <- (-sum_S[i,j,(t-1)]*(beta[i,j] * sum_I[i,j,(t-1)] + alpha[i,j]) )/N.hh.members[i]
      dE[i,j,t] <- (sum_S[i,j,(t-1)]*(beta[i,j] * sum_I[i,j,(t-1)] + alpha[i,j]) - sum_E[i, j,(t-1)]*delta[i,j])/N.hh.members[i]
      dI[i,j,t] <- (sum_S[i,j,(t-1)]*(beta[i,j] * sum_I[i,(t-1)] + alpha[i,j]) - sum_I[i,(t-1)] * epsilon[i,j])/N.hh.members[i]
      dR[i,j,t] <- (sum_I[i,j,(t-1)] * epsilon[i,j])/N.hh.members[i]
      
      S[i,j,t] <- dS[i,j,t] + S[i,j,(t-1)]
      E[i,j,t] <- dE[i,j,t] + E[i,j,(t-1)]
      I[i,j,t] <- dI[i,j,t] + I[i,j,(t-1)]
      R[i,j,t] <- dR[i,j,t] + R[i,j,(t-1)]
      NewI[i,j,t] <- S[i,j,(t-1)]*(beta[i,j] * I[i,j,(t-1)] + alpha[i,j])
      CumNewI[i,j,t] <- NewI[i,j,(t-1)] + NewI[i,j,t]
    }
    SecondaryI[i,,] <- sum(NewI[i,,])
     }
     
   for(t in 1:43){
     sum_S[i,t] <- sum(S[i,,t])
     sum_E[i,t] <- sum(E[i,,t])
     sum_I[i,t] <- sum(I[i,,t])
     sum_R[i,t] <- sum(R[i,,t])

   }

}



mu1 ~dnorm(0,1e-4)
mu2 ~dnorm(0,1e-4)

tau1 ~dgamma(0.001, 0.001) #Check prior!
tau2 ~dgamma(0.001, 0.001)

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




}
"



#Note this doesn't make sense because we aren't taking into account observed time of PCR
zelner_jags2 <- "
model{
for(i in 1:N.HH){

   for(t in 1:tmax[i]){
      log_prob_no_inf_t[i,t] <- log( 1 - (beta * I[i,t] + alpha) )
   }
    
     for(j in 1:N.hh.members[i]){
        y[i,j] ~ dbern(q[i,j]) #y=UNinfected matrix
     
     ##NOTE THESE ARE PROBABLY NOT RIGHT--PROB NEED TO FLIP SOME OF THEM AROUND
      prob_no_inf_uninf[i,j] <- exp(sum(log_prob_no_inf_t[i,1:day.matrix[i,j]])) #P no infections over all time intervals
      
      prob_no_inf_inf_person[i,j] <- exp(sum(log_prob_no_inf_t[i,1:((day.matrix[i,j] - 1)) ])) #P no infections prior to infection
      
      prob_inf[i,j] <- infected_matrix[i,j] * (beta * I[i,day.matrix[i,j]] + alpha) #for infected people only

    q[i,j] <- (1 - prob_no_inf_inf_person[i,j]) * (1 - prob_inf[i,j]) * infected_matrix[i,j] + #prob for infection at time t and not before
                    (1 - prob_no_inf_uninf[i,j]) * (1- infected_matrix[i,j] ) + 1e-6 #prob for uninfected peopel    
}
    

}

##### Attaching latent piece

for(i in 1:N.HH){
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
    dS[i,1] <-0
    dE[i,1] <-0
    dI[i,1] <-0
    dR[i,1] <-0

#NOTE in Zelner, dS, dE, dI,, dR are whole people; as it is here, it is fractional
    for(t in 2:43){
      dS[i,t] <- -S[i,(t-1)]*(beta * I[i,(t-1)] + alpha)
      dE[i,t] <- S[i,(t-1)]*(beta * I[i,(t-1)] + alpha) - E[i, (t-1)]*delta[i]
      dI[i,t] <- E[i, (t-1)]*delta[i] - I[i,(t-1)] * epsilon[i]
      dR[i,t] <- I[i,(t-1)] * epsilon[i]
      
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

 for(i in 1:N.HH){
      for(j in 1:N.hh.members[i]){

     alpha[i,j] <- exp(log.alpha[i,j])
     beta[i,j]  <- exp(log.beta[i,j])
     log.alpha[i,j] ~ dnorm(mu1,tau1)
     log.beta[i,j] ~ dnorm(mu2,tau2)
     
     mu1 ~dnorm(0,1e-4)
    mu2 ~dnorm(0,1e-4)
    
    tau1 ~dgamma(0.001, 0.001) #Check prior!
    tau2 ~dgamma(0.001, 0.001)

      }
 }
 


"


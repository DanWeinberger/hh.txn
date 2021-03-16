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

 
 #Check this hyperprior structure..will result in different values being drawn for each person inHH
     alpha[i,j] <- exp(log.alpha[i,j])
     beta[i,j]  <- exp(log.beta[i,j])
     log.alpha[i,j] ~ dnorm(mu1,tau1)
     log.beta[i,j] ~ dnorm(mu2,tau2)
     
    Dur.Inf[i,j] ~ dgamma(sh1, ra1) #duration infectiusness
    Dur.Latent[i,j] ~dgamma(sh2, ra2) #Duration latent

    delta[i,j] <- 1/Dur.Inf[i,j]
    epsilon[i,j] <- 1/Dur.Latent[i,j]
    
    ##In this formulation, S,E,I,R are defined by individual; sum_S, sum_E, sum_I are defined by household
    S[i,j,1] <- 1-index.case[i,j] #input 0/1
    E[i,j,1] <- 0
    I[i,j,1] <- index.case[i,j]
    R[i,j,1] <- 0


    for(t in 2:(tmax[i]+1)){
      #Define current state
      
      p.S[i,j,t] <- S[i,j,(t-1)] * (1-(beta[i,j] * (sum_I[i,(t-1)]-I[i,j,(t-1)]) + alpha[i,j])) #Prob remain in S
      p.E[i,j,t] <- E[i,j,(t-1)] * (1-epsilon[i,j]) #Prob stay in E
      p.I[i,j,t] <- I[i,j,(t-1)] * (1-delta[i,j])   #Prob stay infected
      
      S.Draw[i,j,t] ~ dbern(p.S[i,j,t]) 
      E.Draw[i,j,t] ~ dbern(p.E[i,j,t]) 
      I.Draw[i,j,t] ~ dbern(p.I[i,j,t]) 

      S[i,j,t] <- S[i,j,(t-1)] * S.Draw[i,j,t] #previously S, draw a 0 or 1
      
      E[i,j,t] <- S[i,j,(t-1)]*(1-S.Draw[i,j,t]) + #previously S, current draw S=0 
                  E[i,j,(t-1)]*E.Draw[i,j,t]  #Previously E, stay E?
      
      I[i,j,t] <- E[i,j,(t-1)]*(1-E.Draw[i,j,t]) + #previously E, current draw E=0 
                  I[i,j,(t-1)]*I.Draw[i,j,t]  #Previously I, stay I? 
      
      R[i,j,t] <- R[i,j,t-1] + I[i,j,(t-1)]*(1-I.Draw[i,j,t]) #previously I, current draw I=0 

    }
     }
     
   for(t in 1:tmax[i]){
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



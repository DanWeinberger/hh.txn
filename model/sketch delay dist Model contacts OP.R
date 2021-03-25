zelner_jags <- "
model{
for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){

    y[i,j] ~ dbern(q[i,j]) #(Person j in HH i: 1= uninfected; 0= infected)

    #day.matrix=day of test for the person

    day.infectious[i,j] <- day.matrix[i,j] - infect.dist[i,j]  #infectious prior to test (t_onset in the document)
    day.exposed[i,j] <- day.infectious[i,j] - expose.dist[i,j]  #latent period (t_infect in the document)
    day.infectious.end[i,j] <- day.infectious[i,j] + end.inf[i,j]  #how long infectious after onset of infectiousness

    infect.dist[i,j] ~ dgamma(sh1, ra1) #duration infectiousness
    expose.dist[i,j] ~ dgamma(sh2, ra2)
    end.inf[i,j] ~ dgamma(sh3,ra3)

    alpha[i,j] <- delta[1] + delta[2]*vax[i,j]  ### Probability of infection from the community (b_ij in the document)
    #could add in things like HH-level random intercept, or intercept to link matched households

    }
  }

##Loop over the HH index cases
for(i in 1:N.HH){ 
  for(j in 1:N.hh.members[i]){ #j is the contact
    for(m in 1:N.hh.members[i]){ #m is the index

      dur.inf[i,j,m] <- y2[i,j]*(day.infectious.end[i,m] - day.infectious[i,m]) # Duration  of infectious period for other infected members (=0 for uninfected people)
      # It corresponds to d in the document

      dur.inf.contact[i,j,m] <- (1-y2[i,j])*step(day.exposed[i,j]-day.infectious[i,m])*
      (min(day.exposed[i,j]-day.infectious[i,m],day.infectious.end[i,m] - day.infectious[i,m])) ### It corresponds to d_inf in the document
      
      #step(abs(j-m)-0.5) ensures don't count the current j in m
      log.prob.uninf.contact[i,j,m]= step(abs(j-m)-0.5) *log(
      y2[i,j]* (1-p_inf[i,j,m])^dur.inf[i,j,m] +     #For uninfected person j
      step(( day.exposed[i,j])-min(day.exposed[i,1:N.hh.members[i]])-0.5)*(1-y2[i,j])* (1-p_inf[i,j,m])^(dur.inf.contact[i,j,m])) #For non-index case in the household
      ### This quantity represents the probability to escape infection from the HH for uninfected and infected (excluding the index) people respectively

      log.prob.uninf.contact.day.inf[i,j,m] =step(abs(j-m)-0.5) * log((
      (1-y2[i,j]) * (1-p_inf[i,j,m])^(step(day.exposed[i,j] - day.infectious[i,m] )*(step(day.infectious.end[i,m] - day.exposed[i,j] ) ) )
      )) #### This quantity represents the probability of NOT escaping infection from the HH for infected contacts (excluding the index)

      p_inf[i,j,m] <- exp(log.p_inf[i,j,m])
      log.p_inf[i,j,m] <- beta[1,i,j] + beta[2,i,j]*vax[i,j] + kappa[1,i,m] + kappa[2,i,m]*vax[i,m]
      }   #### Probability of infection from HH member m (p in the document)

    prob.uninf[i,j] <- exp(sum(log.prob.uninf.contact[i,j,1:N.hh.members[i]]))   ### Take the sum over all HH members and then exponentiate to go back to original scale
    prob.inf.day.inf[i,j] <- (1-exp(sum(log.prob.uninf.contact.day.inf[i,j,1:N.hh.members[i]]))+alpha[i,j])    ### Take the sum over all HH members and then exponentiate to go back to original scale

    q[i,j] <- (1-y2[i,j]) * (1 - alpha[i,j])^(step((day.exposed[i,j]) -min(day.exposed[i,1:N.hh.members[i]]) - 0.5)) * prob.uninf[i,j] * prob.inf.day.inf[i,j] + #for infected contacts (including the index)
    y2[i,j] *(1 - alpha[i,j]) * prob.uninf[i,j]+ 1e-6               #for uninfected person
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

mu1 ~dnorm(0,1e-4)
mu2 ~ dnorm(0,1e-4)
tau1 ~ dgamma(0.001,0.001)
tau2 ~ dgamma(0.001,0.001)

  for(k in 1:2){
    delta[k] ~dnorm(0,1e-4)
    epsilon[k] ~dnorm(0,1e-4)

    for(i in 1:N.HH){ 
      for(j in 1:N.hh.members[i]){ #j is the contact
        beta[k,i,j] ~ dnorm(mu1, tau1) 
        kappa[k,i,j]  ~dnorm(mu2, tau2)
        }
    }
  }
}
"


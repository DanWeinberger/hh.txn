model_string <- "

model{
  p[1] = 0 #Probability of transmission when NOT exposed to index case
  tp[1] <- log(1 - p[1]) #Log probability of NO transmission on day d
  
  for (d in 2:14) {
    p[d] ~ dbeta(1, 1)I(0, 1) #Uninformative Beta prior for probability of transmission on day d-1 of index case illness
    tp[d] <- log(1 - p[d])
  }
  
  for (i in 1:N_contacts) {
    #Index for number of contacts in our dataset
    nY[i] ~ dbern(q[i]) #nY is our data on whether the contact was NOT infected (0=infected, 1=not infected),
  }                       #which is Bernoulli distributed with probability q (defined later)

  for (i in 1:N_contacts) {
    for (t in 1:N_times) {
      #Time window over which we track exposure (infected contacts symptomatic at t=39.5)
      for (j in 1:N_index_cases) {
        #Number of index cases to whom contact I might have been exposed
        te[i , t , j] <-
          tp[1 + T[i , t , j]] #Log likl i escaped infection from j at time t (function of day of illness of j at t)
      } 
      
      loge[i , t] <-
        sum(te[i , t , ]) #Log likelihood i escaped infection from all contact at time t
    }
    for (t in 2:N_times) {
      ploge[i , t] <-
        ploge[i , t - 1] + loge[i , t - 1] #Log likl i escaped infection from all contact prior to time t
    }
    
    ploge[i , 1] <-
      0 #Log likl of escaping infection prior to t=1 (Note the order you define things doesn’t matter)
  }
  
  for (i in 1: N_infected) {
    #Index for infected contacts
    for (t in 1:N_times) {
      inf[infected_indices[i] , t] <-
        (exp(ploge[infected_indices[i] , t])) * (1 - exp(loge[infected_indices[i] , t])) #Likelihood i was infected at time t (and escaped prior to t)
      g[infected_indices[i] , t] <-
        v[40 - t] #Likelihood of incubation period being 39.5 - t days
    }
    
    q[infected_indices[i]] <-
      1 - inprod(g[infected_indices[i] , ], inf[infected_indices[i] , ]) #Likelihood infected contact escaped infection
  }
  
  for (j in 1:N_uninfected) {
    q[uninfected_indices[j]] <-
      exp(ploge[uninfected_indices[j] , N_times]) #Likelihood uninfected contact escaped infection
  }
  
  a ~ dgamma(.001, .001) #Uninformative priors for hyper-parameters of gamma-distributed incubation period
  b ~ dgamma(.001, .001)
  
  for (t in 2:42) {
    v[t] <-
      exp((a - 1) * log(b * (t - 1.5)) - b * (t - 1.5) - loggam(a)) * b #Incubation period is gamma distributed with
  } #hyper-parameters a and b
  
  for (i in 1:81) {
    Q[i] <-
      sum(expos[i , ]) #Likelihood of all possible incubation periods for 81 contacts with known exposure period
  }
  
  for (i in 1:81) {
    obs[i] ~ dbern(Q[i])     #obs is the infection status of cases with known incubation period (=1)
  }
  for (i in 1:81) {
    for (j in 1:41) {
      expos[i , j] <-
        v[incub[i , j] + 1] / sum(v[]) #Likelihood incubation period for i is incub[i,j] days
    }
  }
  v[1] <- 0
}
"
#Data inputs: T, nY, obs, incub, d
#Parameters estimated: p (probability of transmission on day d of illness), a &amp; b
#(posterior estimates for hyper-parameters of gamma-distributed incubation period)

#prob_transmission_function <- function(){


#Notes: we have date of test, need a distribution to account for date of infection
##38 days max is from old paper--might want to have earlier

#Dirichlet for incubation?


model_string <- "

model{
{
  d <- T+1

    p[1:N_contacts,1,] <- 0 #Probability of transmission when NOT exposed to index case
    log_p_uninf_j[1:N_contacts,1,] <- log(1 - p[1:N_contacts,1,]) #Log probability of NO transmission on day d
    cum_log_p_uninf[1:N_contacts , 1] <- 0 #Log likl of escaping infection prior to t=1 (Note the order you define things doesn’t matter)

  for(i in 1: N_contacts){

         nY[i] ~ dbern(q[i]) #nY is our data on whether the contact was NOT infected (0=infected, 1=not infected),
        
         q[i] <- (1-nY[i])*exp(log_p_uninf[i , N_times[i]]) + #Likelihood uninfected contact escaped infection
                 nY[i]*(1 - inprod(g[i , ],inf[i , ])) #Likelihood infected contact escaped infection

  
   for(t in 1: N_times[i]){
   
      inf[i , t] <- (exp(cum_log_p_uninf[i , t])) * (1 - exp(log_p_uninf[i , t])) #Likelihood i was infected at time t (and escaped prior to t)
      g[i , t] <- v[40 - t] #Likelihood of incubation period being 39.5 - t days
 
    for(j in 1:N.indexes[i]){
      #Define prob infection for person i, from contact j at time t
      ## Is a function of the baseline prob for delay for index j and
      ## the vaccine status of the index and the contact at time t
      ## Is d= T+1?
      
      logit_p[i,t,j] <- (alpha[d[t,i,j]] + #Baseline prob probability of transmission on day d-1 of index case illness
                 beta[1]*vax1.index[t,i,j] +  beta[2]*vax2.index[t,i,j] + #Effect of vaccination of the index
                 beta[3]*vax1.contact[t,i] +   beta[4]*vax2.contact[t,i]    ##effect of vaccinaion of the contact       
      )
      
      p[i,t,j] <- exp(logit_p[i,t,j])/exp(logit_p[i,t,j] + 1)  #Inverse logit
       
    }
   }
  
    for( t in 2: N_times[i] ) {
      log_p_uninf_j[i,t, j] <- log(1 - p[i,t,j]) #log(Prob NOT infected)
      log_p_uninf[i , t] <- sum(log_p_uninf_j[i , t ,1:N.contacts[i] ]) #Log likelihood i escaped infection from all contact at time t
      cum_log_p_uninf[i , t] <- cum_log_p_uninf[i , t - 1] + log_p_uninf[i , t-1] #Log likl i escaped infection from all contact prior to time t
    }
  }
  
  a ~ dgamma(.001,.001) #Uninformative priors for hyper-parameters of gamma-distributed incubation period
  b ~ dgamma(.001,.001)
  
  for(d in in 1:14){
    alpha[d] ~dnorm(0,1e-4)
  }
  for(k in 1:4){
  beta[k] ~dnorm(0,1e-4)
  }

  #v[t] ~ #Incubation period


}"
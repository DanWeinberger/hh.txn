#Data inputs: T, nY, obs, incub, d
#Parameters estimated: p (probability of transmission on day d of illness), a &amp; b
#(posterior estimates for hyper-parameters of gamma-distributed incubation period)

#prob_transmission_function <- function(){


#Notes: we have date of test, need a distribution to account for date of infection
## Or does this cancel out because we have date of test for both index and the contact, assuming they have same delay distribution from infection to test
##38 days max is from old paper--might want to have earlier

##Need to add in community transmission rate
## We will miss the 

#Refs:
#LancetID study covid-19 Li https://www.thelancet.com/cms/10.1016/S1473-3099(20)30981-6/attachment/aa38eb0c-395b-4a9a-93de-5977ce4b6a1b/mmc2.pdf

model_string <- "
model{

  for(i in 1: N_contacts){
  
    Uninf.state[i] ~ dbern(q[i]) #Uninf.state is our data on whether the contact was NOT infected (1=uinfected, 0=infected),
    
    q[i] <- Uninf.state2[i]*exp(cum_log_p_uninf[i , N_times[i]]) + #Likelihood uninfected contact escaped infection
            (1- Uninf.state2[i])*(1 - inprod(g[i , ],inf[i , ]))   #Likelihood infected contact escaped infection
  
    for(t in 1: N_times[i]){
    
      inf[i , t] <- exp(cum_log_p_uninf_lag1[i , t]) * (1 - exp(log_p_uninf[i , t])) #Likelihood i was infected at time t (and escaped prior to t)
      g[i , t] <- v[t] #Likelihood of incubation period being  t days--CHECK THIS
    
          #Define prob infection for person i, from contact j at time t
          ## Is a function of the baseline prob for delay for index j and
          ## the vaccine status of the index and the contact at time t
        for(j in 1:N.indexes[i]){
          logit_p[i,t,j] <- (alpha[(T[i,t,j]+1) ] 
            #just baseline
            #+ #Baseline prob probability of transmission on day d-1 of index case illness
            #beta[1]*vax1.index[t,i,j] +  beta[2]*vax2.index[t,i,j] + #Effect of vaccination of the index
            #beta[3]*vax1.contact[t,i] +   beta[4]*vax2.contact[t,i]    ##effect of vaccinaion of the contact       
           )
      
            #When t=1, p[i,t,j]=0
            p[i,t,j] <- (t>1)*exp(logit_p[i,t,j])/(exp(logit_p[i,t,j]) + 1)  #Inverse logit
            log_p_uninf_j[i,t, j] <- log(1 - p[i,t,j]) #log(Prob NOT infected)
        }
    log_p_uninf[i,t] <- sum(log_p_uninf_j[i , t ,1:N.indexes[i] ]) #Log likelihood i escaped infection from all contact at time t
    }
  
    cum_log_p_uninf[i , 1] <- 0
    cum_log_p_uninf_lag1[i , 1] <- 0
  
    for(t in 2:N_times[i]) {
      cum_log_p_uninf[i , t] <- cum_log_p_uninf[i , t - 1] + log_p_uninf[i , t-1] #Log likl i escaped infection from all contact prior to time t
      cum_log_p_uninf_lag1[i , t] <- cum_log_p_uninf[i , t-1] 
    }
  }
  
  
  for(d in 1:max.t){
    alpha[d] ~dnorm(0,1e-4)
  }
  
  for(k in 1:4){
    beta[k] ~dnorm(0,1e-4)
  }

  ###########################


  
}"
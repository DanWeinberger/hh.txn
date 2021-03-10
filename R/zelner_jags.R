zelner_jags <- "
model{

  #beta,S,I,tinfect
  for (t in 1:tmax){
    for(j in 1:hh){
      I[t,j] ~ dmulti(p[t,j], hh.size[j])
    }
  }
  #logit_p ~ dnorm(0,1e-4)
  beta ~ dnorm(0,1e-4)
  #p <- exp(logit_p)/(1+exp(logit_p))


  for (j in 1:hh){ # j represents each household
    for (t in 1:tmax){ #  i in equation
    ##### First chunk of the equation in Lab 7 #####
    # log likelihood contribution of uninfected people in household j
      logit_p[t,j] = logit_p[t,j] - Sf[j]*beta*S[t,j]*I[t,j]*dt 
    }

    ##### Second chunk of the equation in Lab 7 #####

   # if (ninfect[j]>0){  # if the household has at least one case
      for (k in 1:ninfect[j]){ # for each infected person (k) in household j
       # if (tinfect[k,j]>0){
          for (t in 1:(tinfect[k,j]-1)){
          #log-likelihood for escaping prior to time of infection
            logit_p[t,j] = logit_p[t,j] - beta*S[t,j]*I[t,j]*dt 
          }
        # log likelihood they did NOT escape at time of infection
        # (Because we are using discrete time SIR model, we have to use CDF
        #  rather than PDF.)
        logit_p[tinfect[k,j],j] = logit_p[tinfect[k,j],j] + log(1-exp(-beta*S[tinfect[k,j],j]*I[tinfect[k,j],j]*dt)) 
        }
      }
    #}
  #}
  for(t in 1:tmax){
    for(j in 1:hh){
    p[t,j] <- exp(logit_p[t,j])/(exp(logit_p[t,j]) + 1)  #Inverse logit
    }
  }
}
"
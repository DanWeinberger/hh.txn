################################################################################
# EMD 538 Quantitative Methods for Infectious Disease Epidemiology (FALL 2020) #
#                                                                              #
#          Lab: 7 (Data augmentation and Hidden Markov Models)                  #
#         Date: Friday, October 16, 2020                                       #
#     Coded by: Kayoko Shioda & Maile Phillips (maile.phillips@yale.edu)       #
################################################################################

#------------------------------------------------------------------------------#
# Create a function hhsirLL()
#------------------------------------------------------------------------------#

# Equation can be found in a Word document for Lab 7 and Slide 15-16 in Lecture 7 

hhsirLL = function(beta,S,I,tinfect){ # "beta" is a parameter we want to estimate. S, I, & tinfect are the data we have.
  hh = dim(S)[2]   # number of households   = number of columns in the matrix S (50 HHs in Lab 7)
  tmax = dim(S)[1] # duration of follow-up  = number of rows in the matrix S (504 hours in Lab 7)
  dt = 1/24 # time step (in days, i.e. one-24th day=hour, because the questions asks for the rate of transmission "per day", but our data are "per hour")

  S0 = S[1,] # Initial number of susceptibles in each household
  Sf = S[tmax,] # Final number of susceptibles in each household
  ninfect = S0 - Sf # number infected in each household 

  logLL = matrix(0,tmax,hh)
  
  for (j in 1:hh){ # j represents each household
    for (t in 1:tmax){ #  i in equation
      
      ##### First chunk of the equation in Lab 7 #####
      
      # log likelihood contribution of uninfected people in household j
      logLL[t,j] = logLL[t,j] - Sf[j]*beta*S[t,j]*I[t,j]*dt 
    }
    
    ##### Second chunk of the equation in Lab 7 #####
    
    if (ninfect[j]>0){ # if the household has at least one case
      for (k in 1:ninfect[j]){ # for each infected person (k) in household j
        if (tinfect[k,j]>0){
          for (t in 1:(tinfect[k,j]-1)){
            # log-likelihood for escaping prior to time of infection
            logLL[t,j] = logLL[t,j] - beta*S[t,j]*I[t,j]*dt 
          }
          # log likelihood they did NOT escape at time of infection
          # (Because we are using discrete time SIR model, we have to use CDF
          #  rather than PDF.)
          logLL[t,j] = logLL[t,j] + log(1-exp(-beta*S[tinfect[k,j],j]*I[tinfect[k,j],j]*dt)) 
        }
      }
    }
  }
  return(-sum(logLL))
}

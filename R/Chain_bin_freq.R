library(boot)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
library(stats4)

source('./R/simulate_data.R')
source('./R/delay_dist_sim.R')
source('./R/Chain_bin_lik.R')
source('./R/data_manipulation.R')

#### Likelihood definition for all HH
model.run <- function(N.HH){
###### Set random values for the parameters
  alpha0=0.4
  delta0=0.7
  beta=0.2
  kappa=0.5
  params <- c(alpha0,delta0,beta,kappa)

  #Generate the data and store as a data frame
  sim.data.ls <- pblapply(1:N.HH,gen.hh)
  #This is like the data we would get from KSM
  sim.data.df <- do.call('rbind.data.frame', sim.data.ls)
  ##
  Y <- delay.gen(input_df = sim.data.df)[[1]]
  X <- delay.gen(input_df = sim.data.df)[[2]]
  chain_bin_lik(params,Y,X)
  #optim(params,chain_bin_lik,Y=Y,X=X)
}


#test1 <- df4[,c('hhID', 'ID','ID_b', 't.index','Y')]

##1. Run delay.gen() to create X and Y for a single set of random delay dist value
##2. Run mle() or optim() to estimate parameters
##3 Repeat 1-2 N times
##4 average of all N

#Need Y to represent each HH and time, from 0 to censor time; 
#Check that simulated data truncated to 21 days -> this goes from 1 to ncol(infect.status) so correct




#mod.data <- pbreplicate(10000,delay.gen(sim.data.df))

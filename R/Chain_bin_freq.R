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

#Generate the synthetic data and store as a data frame
N.HH <- 500
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), prob.trans.day=(1-0.968),irr.vax=0.2)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)


#### Likelihood definition for all HH
model.run <- function(ds){
###### Set random values for the parameters
  alpha0= log(1-0.968)
  delta0= log(1-0.9995)
  beta= log(0.2)
  kappa= log(1)
  params <- c(alpha0,delta0,beta,kappa)

  ##
  LatentData <-  delay.gen(input_df = ds)
  Y <- LatentData$Y
  X <- LatentData$X
  #chain_bin_lik(params,Y,X)
  optim(params,chain_bin_lik,Y=Y,X=X, method='BFGS')
 # mle(chain_bin_lik, start=params)
}

mod.data <- pbreplicate(2,model.run(ds=sim.data.df), simplify=F)

parms <- sapply(mod.data,'[[','par')

hist(parms[1,])
hist(parms[2,])
hist(parms[3,])
hist(parms[4,])

#test1 <- df4[,c('hhID', 'ID','ID_b', 't.index','Y')]

##1. Run delay.gen() to create X and Y for a single set of random delay dist value
##2. Run mle() or optim() to estimate parameters
##3 Repeat 1-2 N times
##4 average of all N

#Need Y to represent each HH and time, from 0 to censor time; 
#Check that simulated data truncated to 21 days -> this goes from 1 to ncol(infect.status) so correct




#mod.data <- pbreplicate(10000,delay.gen(sim.data.df))

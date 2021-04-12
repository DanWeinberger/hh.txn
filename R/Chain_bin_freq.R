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
N.HH <- 5000
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), prob.trans.day=(1-0.968),irr.vax1=0.5,irr.vax2=1)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)

#### Likelihood definition for all HH
model.run <- function(ds){
###### Set random values for the parameters
  alpha0=0
  delta0= 0
  beta= 0
  kappa= 0
  params <- c(alpha0,delta0,beta,kappa)

  ##
  LatentData <-  delay.gen(input_df = ds)
  Y <- LatentData$Y
  X <- LatentData$X
  #chain_bin_lik(params,Y,X)
  optim(params,chain_bin_lik,Y=Y,X=X, method='BFGS')
 # mle(chain_bin_lik, start=params)
}

mod.data <- pbreplicate(1,model.run(ds=sim.data.df), simplify=F)

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

############################################
#clock the pieces of the likelihood function
#############################################

#### Define logit_p = X*params; need  to add as.matrix
ptm <- proc.time()
logit_p <- as.vector(as.matrix(X[,1:4]) %*% params) ## Added as.matrix
proc.time() - ptm

### Go back to p (probability  of transmission) with inverse logit: 
ptm <- proc.time()
q <- 1 - 1/(1 + exp(-logit_p))
proc.time() - ptm

##Pi needs to be a single value by ID/hhID/time point; Y should be same length
ptm <- proc.time()
  data_tab <- cbind.data.frame(log.q=log(q),X)
  data_t = data.table(data_tab)
  ans = data_t[,list(A = sum(log.q)), by = 'ID,hhID,t.index']
  pi=ans$A
proc.time() - ptm

#Old, slow version
# ptm <- proc.time()
# pi <- 1 - exp(aggregate(log(q), by=list(X$ID, X$hhID, X$t.index ), FUN=sum)$x )
# proc.time() - ptm

### Likelihood definition (for the moment no log-lik, so there is just a product over all HH members and time steps):
ptm <- proc.time()
ll= sum(dbinom(x=Y,size=1,prob = pi,log = TRUE),na.rm = TRUE)
proc.time() - ptm



#mod.data <- pbreplicate(10000,delay.gen(sim.data.df))

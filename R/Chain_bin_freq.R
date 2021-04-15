library(boot)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
library(stats4)
library(dplyr)

source('./R/simulate_data.R')
source('./R/delay_dist_sim.R')
source('./R/Chain_bin_lik2.R')
source('./R/data_manipulation_df.R')

#Generate the synthetic data and store as a data frame
N.HH <- 20000
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), prob.trans.day=(1-0.968),irr.vax1=0.5,irr.vax2=1)


#How many infections per household
hh.inf <- sapply(sim.data.ls, function(x) sum(x$infected))

sim.data.ls.pos <- sim.data.ls[hh.inf>0]

sim.data.ls.pos.sub <- sim.data.ls.pos[1:1000]


sim.data.df <- do.call('rbind.data.frame', sim.data.ls.pos.sub)

LatentData <-  delay.gen(input_df = sim.data.df)

params <- c(alpha0=0,delta0=0,beta=0,kappa=0)
ptm <- proc.time()
mod1 <- optim(params,chain_bin_lik,Y=LatentData$Y,X=LatentData$X, method='BFGS')
proc.time() - ptm

parms <- sapply(mod.data,'[[','par')
parms










############################################
#clock the pieces of the likelihood function
#############################################

alpha0=0
delta0= 0
beta= 0
kappa= 0
params <- c(alpha0,delta0,beta,kappa)

##
ptm <- proc.time()
LatentData <-  delay.gen(sim.data.df)
proc.time() - ptm

X.ls <- LatentData$X.ls
Y.ls <- LatentData$Y.ls

nrowX <- sapply(X.ls, nrow)


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

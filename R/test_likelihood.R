library(reshape2)
source("./R/Chain_bin_lik.R")
###################
### Set the true parameter values
alpha0_true= log(0.3)
delta0_true=log(0.08)
beta_true= log(0.5)
kappa_true= log(0.09)
params_true <- c(alpha0_true,delta0_true,beta_true,kappa_true) 

### Define the functions used in this file
###Simulate  people
gen.hh.test <- function(idN,prop.vax1=0.5, prop.vax2=0.5){
  HH.size <- min(1+ rpois(n=1,1.5),5) #cap at 5
  df1 <- as.data.frame(matrix(NA, nrow=HH.size, ncol=2))
  names(df1) <- c('ID', 'hhID')
  df1$hhID <- idN
  df1$ID <- 1:HH.size
  df1$vax1dose <- rbinom(nrow(df1), 1,prop.vax1)
  df1$vax2dose <- rbinom(nrow(df1), 1,prop.vax2)*df1$vax1dose
  return(df1)
}


data.manipulation.test <- function(input_df){
  sim.data.df.spl <- split(input_df, input_df$hhID)
  df1 <- do.call('rbind.data.frame',sim.data.df.spl)
  df1$max.time.hh <- 15

  ### Create a  copy of df1 and call it df2
  df2 <- df1

  ### Rename all of the variables on df2 as HHID_b, ID_b, Vax_dose1b
  colnames(df2) <- paste(colnames(df2), 'b', sep='_')

  ### Merge df1 and df2 by HHID; name this df3 
  df3 <- merge(df1,df2,by.y='hhID_b' , by.x='hhID')
  df3$expand.t <- NA
  df3$expand.t <- df3$max.time.hh #censor uninfected person at max.time.hh
  #Add an index for each person
  df3$rowN <- 1:nrow(df3)
  df.t <- df3[rep(df3$rowN, times=df3$expand.t),] #expands the df to number of time points

  #Add an index for each person
  df.t$t.index <- unlist(lapply(df3$expand.t, function(x){ 
    z <- 1:x 
    return(z)
  }))

  df3 <- df.t

  #Create elements for design matrix
  df3$delta0 <- 0
  df3$delta0[(df3$ID == df3$ID_b)] <- 1 #exogenous/community risk 

  df3$alpha0 <- 0
  df3$alpha0[df3$ID != df3$ID_b ] <- 1

  df3$vax1 <- df3$vax1dose
  df3$vax2 <- df3$vax1dose_b
  df3$vax2[df3$alpha0==0] <- 0

  #Design matrix 
  X <- df3[c('alpha0','delta0','vax1','vax2','ID','hhID','t.index')]

  ### Generate the probabilities
  #### Define logit_p = X*params; need  to add as.matrix
  logit_p <- as.vector(as.matrix(X[,1:4]) %*% params_true) 
  ### Go back to p (probability  of transmission) with inverse logit: 
  q <- 1 - 1/(1 + exp(-logit_p))
  ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
  data_tab <- cbind.data.frame(log.q=log(q),X)
  data_t = data.table(data_tab)
  ans = data_t[,list(A = sum(log.q)), by = 'ID,hhID,t.index']
  ans <- setorder(ans, t.index)
  pi= 1- exp(ans$A)
  ### Generate the infections
  ans$infect.status <- rbinom(n=nrow(ans),size=1,prob = pi)
  time.infection <- ans[ , .SD[which.min(t.index)], by = list(infect.status,hhID,ID)][infect.status==1]$t.index
  ####  Just keep the infections at the first time step
  try1 <- ans[ , .SD[which.min(t.index)], by = list(infect.status,hhID,ID)][infect.status==1]
  try0 <- ans[ans$infect.status==0]
  
  try2  = rbind(try0,try1)
  data_try2 = data.table(try2)
  ans.final = data_try2[,list(A = mean(infect.status)), by = 'ID,hhID,t.index']
  ans.final <- setorder(ans.final, t.index)
  data_tab <- cbind.data.frame(ans.final[,c('A','ID','hhID','t.index')])
  data_t = data.table(data_tab)
  
  #### Generate the Ys
  Y.df = data_t[,list(A = mean(A)), by = 'ID,hhID,t.index']
  Y <-  Y.df$A
  
  out.list=list('Y'=Y, 'X'=X)  

  return(out.list)}

######  Now compute the likelihood using the "observed data":


#Generate the synthetic data and store as a data frame
N.HH <- 1000
sim.data.ls <- pblapply(1:N.HH, gen.hh.test)

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
  LatentData <-  data.manipulation.test(input_df = ds)
  Y <- LatentData$Y
  X <- LatentData$X
  #chain_bin_lik(params,Y,X)
  optim(params,chain_bin_lik,Y=Y,X=X, method='BFGS',hessian = T)
  # mle(chain_bin_lik, start=params)
}

mod.data <- pbreplicate(1,model.run(ds=sim.data.df), simplify=F)


### Check parameters
parms <- sapply(mod.data,'[[','par')

###  Compute SE

hessian <- mod.data[[1]]$hessian
OI<-solve(hessian)
se<-sqrt(diag(OI))


### Compare with true values
parms
se
params_true


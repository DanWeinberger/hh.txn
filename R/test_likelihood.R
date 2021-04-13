library(reshape2)
library(pbapply)
library(data.table)
source("./R/Chain_bin_lik.R")

###################
### Set the true parameter values
alpha0_true= log(0.01)
delta0_true= log(0.08)
beta_true= log(0.09)
kappa_true= log(0.07)
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
  
  #Create elements for design matrix
  df3$delta0 <- 0
  df3$delta0[(df3$ID == df3$ID_b)] <- 1 #exogenous/community risk 

  df3$alpha0 <- 0
  df3$alpha0[df3$ID != df3$ID_b ] <- 1

  df3$vax1 <- df3$vax1dose
  df3$vax2 <- df3$vax1dose_b
  df3$vax2[df3$alpha0==0] <- 0


  
  #Design matrix 
  X <- df3[c('alpha0','delta0','vax1','vax2','ID','hhID')]
  
  ### Set the time at t=0  --> in  the loop it starts computing 
  X$t.index <- 1
  infect.stat <- 0
  #X$infect.status <- infect.stat
  N.hh.members <- aggregate( df1$ID, by=list( 'hhID'=df1$hhID), FUN=length)$x
  mat.inf.final <- list()
  mat.inf  <- aggregate( X$t.index, by=list( 'ID'=X$ID, 'hhID'=X$hhID), FUN=mean)
  colnames(mat.inf)[3] <- 'tindex'
  mat.inf$infect.status <- infect.stat
  for(i in 1:N.HH){
    mat.inf.2  <- mat.inf[mat.inf$hhID==i,]
    for (t in 2:(df1$max.time.hh[1])){
      
        logit_p <- as.vector(as.matrix(X[X$hhID==i,][,1:4]) %*% params_true) 
      ### Go back to p (probability  of transmission) with inverse logit: 
        q <- 1 - 1/(1 + exp(-logit_p))
      ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
        data_tab <- cbind.data.frame(log.q=log(q),X[X$hhID==i,])
        data_t = data.table(data_tab)
        ans = data_t[,list(A = sum(log.q)), by = 'ID,hhID,t.index']
        pi= 1- exp(ans$A)
        ### Generate probabilities of becoming infected
        infect.status <- rbinom(n=N.hh.members[i],size=1,prob = pi)
        ### Store results: if infected status==0, place  the new results, otherwise leave 1
        mat.inf.2$infect.status[1:N.hh.members[i]] <- ifelse(mat.inf.2$infect.status[1:N.hh.members[i]]!=1,infect.status,1)
        ### Increase time index if infect.status==1
        mat.inf.2$tindex[1:N.hh.members[i]] <-ifelse(mat.inf.2$infect.status[1:N.hh.members[i]]!=1, (mat.inf.2$tindex[1:N.hh.members[i]]+1),mat.inf.2$tindex[1:N.hh.members[i]])
        
    }
    mat.inf.final[[i]] <- mat.inf.2
    
  }
  
  mat.inf.df <- do.call('rbind.data.frame',mat.inf.final)
  mat.inf.df.final <- mat.inf.df[rep(seq(nrow(mat.inf.df)),times=mat.inf.df$tindex),1:4]
  #Add an index for each person
  mat.inf.df.final$t.index <- NA
  mat.inf.df.final$t.index <- unlist(lapply( mat.inf.df$tindex, function(x){ 
   z <- 1:x 
   return(z)
  }))
  
 
  mat.inf.df.final$infect.status <- ifelse(mat.inf.df.final$t.index<mat.inf.df.final$tindex,0,1)*ifelse(mat.inf.df.final$infect.status==1,1,0)
  
  
  data_try2 = data.table(mat.inf.df.final)
  ans.final = data_try2[,list(A = mean(infect.status)), by = 'ID,hhID,t.index']
  ans.final <- setorder(ans.final, t.index)
  data_tab <- cbind.data.frame(ans.final[,c('A','ID','hhID','t.index')])
  data_t = data.table(data_tab)
  
  #### Generate the Ys
  Y.df = data_t[,list(A = mean(A)), by = 'ID,hhID,t.index']
  Y <-  Y.df$A
  out.list=list('Y'=Y, 'X'=X,'pi'=pi)  

  return(out.list)}

######  Now compute the likelihood using the "observed data":


#Generate the synthetic data and store as a data frame
N.HH <- 1000
sim.data.ls <- pblapply(1:N.HH, gen.hh.test)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)

#### Likelihood definition for all HH
alpha0=0
delta0= 0
beta= 0
kappa= 0
params <- c(alpha0,delta0,beta,kappa)
model.run <- function(ds, params){
  ###### Set random values for the parameters
  ##
  LatentData <-  data.manipulation.test(input_df = ds)
  Y <- LatentData$Y
  X <- LatentData$X
  pi <- LatentData$pi
  #chain_bin_lik(params,Y,X)
  optim(params,chain_bin_lik,Y=Y,X=X, method='BFGS',hessian = T)
  # mle(chain_bin_lik, start=params)
}

mod.data <- pbreplicate(1,model.run(ds=sim.data.df,params=params), simplify=F)


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


### Test the loglikelihood with true parameters: 

mod.data1 <- pbreplicate(1,model.run(ds=sim.data.df,params=params), simplify=F)
parms <- sapply(mod.data1,'[[','par')
mod.data2 <- pbreplicate(1,model.run(ds=sim.data.df,params=params_true), simplify=F)
mod.data3 <- pbreplicate(1,model.run(ds=sim.data.df,params=parms), simplify=F)


mod.data1[[1]]$value
mod.data2[[1]]$value
mod.data3[[1]]$value


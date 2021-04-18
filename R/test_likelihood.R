library(reshape2)
library(pbapply)
library(data.table)
source("./R/Chain_bin_lik.R")
library(dplyr)
library(tidyr)

###################
### Set the true parameter values
alpha0_true= log(0.5)
delta0_true= log(0.5)
beta_true= log(0.3)
kappa_true= log(0.3)
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
  mat.inf.df <-  mat.inf
  X$infect.status <- infect.stat 
  X_ <- X
  X_HH <- list()

  #Y_final <- list()
  X_final <- list()
  X_final[[1]] <- X
  #rep.inf <- 0
  for (t in 2:df1$max.time.hh[1]){
    #print(t)
   mat.inf <- mat.inf.df 
   for(i in 1:N.HH){
        if((length(X)!=0) & (dim(X[X$hhID==i,])[1]!=0)){
        #print(dim(X[X$hhID==i,])[1]!=0)
        mat.inf.2  <- mat.inf[mat.inf$hhID==i,]
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
       # print(paste0('X',X))
        N.hh.memb.X <- aggregate( X$ID, by=list( 'hhID'=X$hhID), FUN=length)$x
        num.dup <- aggregate(list(numdup=rep(1,nrow(X[X$hhID==i,]))), X[X$hhID==i,], length)$numdup
        #X$infect.status[X$hhID==i] <-  rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],each=length(unique(X[X$hhID==i,]$ID)))#*ifelse(unique(X$ID)==mat.inf.2$ID,1,infect.status)
        #rep.inf <- infect.status[1:length(unique(X[X$hhID==i,]$ID))]
        #print(i)
        #print(paste0('length(rep.inf): ',length(rep.inf)))
        #print(paste0('dim(X[X$hhID==i,])[1]: ',dim(X[X$hhID==i,])[1]))
        #X$infect.status[X$hhID==i] <-  rep.inf*ifelse(dim(X[X$hhID==i,])[1]==length(rep.inf),1,rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],each=length(unique(X[X$hhID==i,]$ID))))
        #X$infect.status[X$hhID==i] <-rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],times=sqrt(N.hh.memb.X[i]))#infect.status[1:length(unique(X[X$hhID==i,]$ID))]
        mat.inf.final[[i]] <- mat.inf.2
        
        if(min(mat.inf.2$infect.status)==0){
          jindex <- mat.inf.2$ID[mat.inf.2$infect.status!=1] 
          HHindex <- mat.inf.2$hhID[mat.inf.2$infect.status!=1]
          X_new <- X
          X_ <- list()
          k=0
          for(j in 1:N.hh.members[i]){
            print(j)
            for(k in 1:length(jindex)){
              print(k)
              if(j==jindex[k]){
              print(j)
              print(k)
              X_new <- X[(X$ID==j) & (X$hhID==i),]
              print(X[(X$ID==j) & (X$hhID==i),])
              X_[[k]] <- X_new
            }
            }
          }
          X_ <- do.call('rbind.data.frame',X_)
          X_HH[[i]] <- X_ 
        }else(X_HH[[i]]<-list())
        
     
     X_HH.df  <- do.call('rbind.data.frame',X_HH)
        }else{break}
   }
   
    mat.inf.df  <- do.call('rbind.data.frame',mat.inf.final)
    #X <- X_HH.df
  if(min(mat.inf.df$infect.status)!=1){
    print(mat.inf.df$infect.status)
    print(t)
    X <- X_HH.df
    X_final[[t]] <- X_HH.df 
    }
  }
  
  
  mat.inf.df.final <- mat.inf.df[rep(seq(nrow(mat.inf.df)),times=mat.inf.df$tindex),1:4]
  #Add an index for each person
  mat.inf.df.final$t.index <- NA
  mat.inf.df.final$t.index <- unlist(lapply( mat.inf.df$tindex, function(x){ 
    z <- 1:x 
    return(z)
  }))
  
  data_try2 = data.table(mat.inf.df.final)
  ans.final = data_try2[,list(A = mean(infect.status)), by = 'ID,hhID,t.index']
  ans.final <- setorder(ans.final, t.index)
  data_tab <- cbind.data.frame(ans.final[,c('A','ID','hhID','t.index')])
  data_t = data.table(data_tab)
  
  #### Generate the Ys
  Y.df = data_t[,list(A = mean(A)), by = 'ID,hhID,t.index']
  Y <-  Y.df$A
  X_t.df  <- do.call('rbind.data.frame',X_final)
  #X.fin  <- do.call('rbind.data.frame',X_final)
  
  #X.spl <- split(X.fin, X.fin$hhID)
  #X.spl.2 <- lapply(X.spl,function(x) split(x,x$ID))
  #max.time.inf  <- aggregate( mat.inf.df.final$tindex, by=list( 'ID'=mat.inf.df.final$ID, 'hhID'=mat.inf.df.final$hhID), FUN=max)
  
  #res <- list()
  #for(l in 1:length(X.spl.2[[1]])){
   # for(k in 1:length(X.spl.2[[l]])){
     # print(k)
     # res[[l]] = X.spl.2[[l]][[k]][X.spl.2[[l]][[k]]$t.index <=max.time.inf$x[l],]
  #}
#}
 # library(data.table)
 #
  #setDT(X.fin)[, .SD[.N > max.time.inf$x], by  =t.index]
  

 # do.call(rbind, lapply(split(df, df$y), function(i) if(nrow(i) >= 3) { i }))
  
  #apply(mat.inf.df,2,function(x) max(x$t.index))
  out.list=list('Y'=Y, 'X'=X_t.df,'pi'=pi)  

  return(out.list)}

######  Now compute the likelihood using the "observed data":
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
  mat.inf.df <-  mat.inf
  X$infect.status <- infect.stat 
  X_ <- X
  X_HH <- list()
  
  #Y_final <- list()
  X_final <- list()
  X_final[[1]] <- X

  mat.inf <- mat.inf.df 
  for(i in 1:N.HH){
    if((length(X)!=0) & (dim(X[X$hhID==i,])[1]!=0)){
        #print(dim(X[X$hhID==i,])[1]!=0)
        mat.inf.2  <- mat.inf[mat.inf$hhID==i,]
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
        # print(paste0('X',X))
        N.hh.memb.X <- aggregate( X$ID, by=list( 'hhID'=X$hhID), FUN=length)$x
        num.dup <- aggregate(list(numdup=rep(1,nrow(X[X$hhID==i,]))), X[X$hhID==i,], length)$numdup
        #X$infect.status[X$hhID==i] <-  rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],each=length(unique(X[X$hhID==i,]$ID)))#*ifelse(unique(X$ID)==mat.inf.2$ID,1,infect.status)
        #rep.inf <- infect.status[1:length(unique(X[X$hhID==i,]$ID))]
        #print(i)
        #print(paste0('length(rep.inf): ',length(rep.inf)))
        #print(paste0('dim(X[X$hhID==i,])[1]: ',dim(X[X$hhID==i,])[1]))
        #X$infect.status[X$hhID==i] <-  rep.inf*ifelse(dim(X[X$hhID==i,])[1]==length(rep.inf),1,rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],each=length(unique(X[X$hhID==i,]$ID))))
        #X$infect.status[X$hhID==i] <-rep(infect.status[1:length(unique(X[X$hhID==i,]$ID))],times=sqrt(N.hh.memb.X[i]))#infect.status[1:length(unique(X[X$hhID==i,]$ID))]
        mat.inf.final[[i]] <- mat.inf.2
    }
    
    mat.inf.df  <- do.call('rbind.data.frame',mat.inf.final)
  }
  
  
  mat.inf.df.final <- mat.inf.df[rep(seq(nrow(mat.inf.df)),times=mat.inf.df$tindex),1:4]
  #Add an index for each person
  mat.inf.df.final$t.index <- NA
  mat.inf.df.final$t.index <- unlist(lapply( mat.inf.df$tindex, function(x){ 
    z <- 1:x 
    return(z)
  }))
  
  data_try2 = data.table(mat.inf.df.final)
  ans.final = data_try2[,list(A = mean(infect.status)), by = 'ID,hhID,t.index']
  ans.final <- setorder(ans.final, t.index)
  data_tab <- cbind.data.frame(ans.final[,c('A','ID','hhID','t.index')])
  data_t = data.table(data_tab)
  
  #### Generate the Ys
  Y.df = data_t[,list(A = mean(A)), by = 'ID,hhID,t.index']
  Y <-  Y.df$A
  X_t.df  <- do.call('rbind.data.frame',X_final)
  out.list=list('Y'=Y, 'X'=X,'pi'=pi)  
  
  return(out.list)}


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

ds=sim.data.df
LatentData <-  data.manipulation.test(input_df = ds)
Y <- LatentData$Y
X <- LatentData$X
pi <- LatentData$pi

#Solution with c(0, 0, 0, 0) starting point
mle<-optim(c(0, 0, 0, 0),
           chain_bin_lik,
           Y=Y,
           X=X, 
           method='BFGS',
           hessian = T)
mle$value
mle$par
params_true
sqrt(diag(chol2inv(chol(mle$hessian))))

#Solution with random starting point
mle<-optim(rnorm(n=4),
           chain_bin_lik,
           Y=Y,
           X=X, 
           method='BFGS',
           hessian = T)
mle$value
mle$par
params_true
sqrt(diag(chol2inv(chol(mle$hessian))))

#Solution with params_true starting point
mle<-optim(params_true,
           chain_bin_lik,
           Y=Y,
           X=X, 
           method='BFGS',
           hessian = T)
mle$value
mle$par
params_true
sqrt(diag(chol2inv(chol(mle$hessian))))




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

LatentData <-  data.manipulation.test(input_df = ds)
Y <- LatentData$Y
X <- LatentData$X
pi <- LatentData$pi
#Evaluating the Likelihood
chain_bin_lik(c(0,0,0,0), Y=Y, X=X)   #All Zeros
chain_bin_lik(params_true, Y=Y, X=X)  #True Parameters
chain_bin_lik(parms, Y=Y, X=X)        #Estimated Parameters

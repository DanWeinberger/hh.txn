library(boot)
library(matrixStats)
library(reshape2)
library(data.table)
source('./R/simulate_data.R')
#### Likelihood definition for all HH


###### Set random values for the parameters
alpha0=0.4
delta0=0.7
beta=0.2
kappa=0.5
params <- c(alpha0,delta0,beta,kappa)

#Generate the data and store as a data frame
list_hh <-lapply(1:10,gen.hh)
hh_list <- lapply(list_hh,function(x) x=x$df1)
infect.status <- lapply(list_hh,function(x) x=x$infect.status)
delay.dist <- lapply(hh_list, function(x) delay_dist_sim(x))
N.HH <- length(hh_list)
max.time.hh.final <- list()
y_hh <- list()
delta_hh <- list()
X_hh <- list()

## Select one HH
for(i in 1:N.HH){
  df1 <- hh_list[[i]]
  infect.status_hh <- infect.status[[i]]

  #### Define delay distributions and max.time.hh
  delay.dist.hh <- delay_dist_sim(df1)
  max.time.hh <- delay.dist.hh$max.time.hh
  max.time.hh.final[[i]] <- max.time.hh
  infect.status_hh <- as.matrix(infect.status_hh[,1:max.time.hh]) ## to avoid that there are infected people after time.study.hh -> this would create
  ### issues with y definition
  day.infectious <- delay.dist.hh$day.infectious
  day.exposed <- delay.dist.hh$day.exposed
  day.infectious.end <- delay.dist.hh$day.infectious.end
  
  #### Definition of delta
  delta <- array(NA,dim=c(nrow(df1)+1,nrow(df1)+1,max.time.hh))
  ### delta is equal to 1 for the first row referring to baseline risk contribution
  delta[1,2:(nrow(df1)+1),]<- 1
  for(j in 2:(nrow(df1)+1)){
    for(t in 1:max.time.hh){
      for(k in 2:(nrow(df1)+1)){
        delta[j,k,t] <-ifelse(k!=j,1,0)*ifelse(t > day.infectious[k-1],1,0)*ifelse(t < min(day.exposed[j-1], day.infectious.end[k-1]),1,0)
      }
    }
  }
  delta_hh[[i]] <- delta
  
  #### Define observed data y which depend on max.time.hh and infect.status
  y <- matrix(NA,nrow=max.time.hh,ncol=nrow(df1))
  #infect.status <- infect.status[[i]]
  for(l in 1:nrow(df1)){
    if(max(infect.status_hh[l,])!=0){
    #y[1:min(0,which(infect.status_hh[l,]==1))-1,l]=0
    #y[min(0,which(infect.status_hh[l,]==1)),l]=1
      y[1:min(which(infect.status_hh[l,]==1))-1,l]=0
      y[min(which(infect.status_hh[l,]==1)),l]=1
    }else{
     y[1:nrow(y),l] <- 0    
    }  
  }
  y_hh[[i]] <- y
  ### Manipulated dfs to extract values for the design matrix X
  
  ### Preprocessing the data: creating y, X and call to delay distributions
  
  ### Create a  copy of df1 and call it df2
  df1 <- hh_list[[i]]
  df2 <- copy(df1)
  ### Rename all of the variables on df1 as HHID_b, ID_b, Vax_dose1b
  colnames(df1) <- paste(colnames(df1), 'b', sep='_')
  colnames(df1)[2] <- 'hhID'
  ### Merge df1 and df2 by HHID; name this df3 
  df3 <- merge(df1,df2,by='hhID')
  vax <- (df3[df3$ID!=df3$ID_b,]$vax1dose)
  
  
  #### Maybe there is a more efficient way of doing so
  #index <- unlist(df1$ID)
  #vax <- list()
  #for(m in 1:nrow(df1)){
  #  if(m==index[m]){
  #    index1 = index[-m]
  #    vax[[m]] <- df1[df1$ID %in% index1, ]$vax1dose
  #  }
  #}
  #vax <- unlist(vax)
  
  
  ### Define the design matrix: first two rows refer to the exogenous contribution, last two to the endogenous contribution
  X = matrix(0, nrow=nrow(df1)^2,ncol=length(params))
  if(nrow(df1)==1){
    X[1:nrow(df1),2] <- 1    ## set delta0
    X[1:nrow(df1),3] <- df1$vax1dose
  }else{
    X[(nrow(df1)+1):nrow(X),1] <-1 ## set alpha0
    X[1:nrow(df1),2] <- 1    ## set delta0
    X[1:nrow(df1),3] <- df1$vax1dose
    X[(nrow(df1)+1):nrow(X),3] <-rep(df1$vax1dose,each=(nrow(df1)-1))
    X[(nrow(df1)+1):nrow(X),4] <- vax 
  }
  X_hh[[i]] <- X
}

chain_bin_lik <- function(params){
  ll=0
  for(i in 1:N.HH){
  #### Define logit_p = X*params
    logit_p <- X_hh[[i]] %*%params
  
  ### Go back to p (probability  of transmission) with inverse logit: 
    p <- exp(logit_p)/(exp(logit_p) + 1) 
  
  ### Create matrix (k+1)*(j+1) with 1-p
    j<-matrix(0,nrow=nrow(delta_hh[[i]]),ncol=nrow(delta_hh[[i]]),byrow=T)
    j[,1]<-NA
    diag(j)<-NA
    j[is.na(j) == 0]<-1-p
    k <- flatten(j)
    k[is.na(k)==0]=1-p
    pe_mat <- matrix(k,nrow=nrow(delta_hh[[i]]),ncol=nrow(delta_hh[[i]]),byrow = T)
  
  #### Create delta's: same dimensions as pe_mat. Once delta is created, define for each time step: pi=(1-pe_mat)^(delta)
    pe_d <- matrix(NA,nrow=nrow(delta_hh[[i]]),ncol=nrow(delta_hh[[i]]))
    pi <- list()
    for(t in 1:max.time.hh.final[[i]]){
    ### define p^d
      pe_d <- colSums(log(pe_mat^(delta_hh[[i]][,,t])),na.rm = T)
      pi[[t]] <- 1-exp(pe_d)
    }

  
  ### Stack list together
    pi<-do.call("rbind",pi)
  
  ### Likelihood definition (for the moment no log-lik, so there is just a product over all HH members and time steps):
    ll_piece= sum(dbinom(x=y,size=1,prob = pi[,2:ncol(pi)],log = TRUE),na.rm = TRUE)
    ll = ll + ll_piece
  }
  return(-ll)
}


###  Simulate delay_distributions

delay_dist_sim <- function(df1){
  n.sim <- nrow(df1)
  expose.dist= rgamma(n.sim,4,0.75) #duration latent
  infect.dist= rgamma(n.sim,4,0.75) #duration infectiousness
  end.inf= rgamma(n.sim,4, 0.75) #duration infectiousness
  
  day.infectious <- round(df1$day_index - infect.dist)
  day.exposed <- round(day.infectious - expose.dist)
  day.infectious.end <-  round(day.infectious + end.inf)
  first.day.hh <- min(day.exposed, na.rm=T) - 1
  last.day.hh <- max(day.infectious.end, na.rm=T)
  
  day.infectious <- day.infectious - first.day.hh
  day.exposed <- day.exposed - first.day.hh
  day.infectious.end <- day.infectious.end - first.day.hh
  max.time.hh <- max(day.infectious.end, na.rm=T)
  
  res.list <- list('max.time.hh'=max.time.hh,'day.infectious'=day.infectious,'day.exposed'=day.exposed,'day.infectious.end'=day.infectious.end)
  return(res.list)
}




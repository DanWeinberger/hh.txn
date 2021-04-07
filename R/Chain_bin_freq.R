library(boot)
library(matrixStats)
source('./R/simulate_data.R')

######
alpha0=0.4
delta0=0.7
beta=0.2
kappa=0.5
params <- c(alpha0,delta0,beta,kappa)


#Generate the data and store as a data frame
list_hh <-lapply(1:10,gen.hh)
hh_list <- lapply(list_hh,function(x) x=x$df1)
infect.status <- lapply(list_hh,function(x) x=x$infect.status)
## Time step =1
df1 <- hh_list[[2]]
infect.status <- infect.status[[2]]



delay.dist <- delay_dist_sim(df1)
max.time.hh <- delay.dist$max.time.hh
day.infectious <- delay.dist$day.infectious
day.exposed <- delay.dist$day.exposed
day.infectious.end <- delay.dist$day.infectious.end


y <- matrix(NA,nrow=max.time.hh,ncol=nrow(df1))
y[1:min(which(infect.status[1,]==1))-1,1]=0
y[min(which(infect.status[1,]==1)),1]=1

y[1:min(which(infect.status[2,]==1))-1,2]=0
y[min(which(infect.status[2,]==1)),2]=1





chain_bin_lik <- function(params,df1){
  
### define the design matrix: first two rows refer to the exogenous contribution, last two to the endogenous contribution
  X = matrix(0, nrow=nrow(df1)^2,ncol=length(params))
  X[(nrow(df1)+1):nrow(X),1] <-1 ## set alpha0
  X[1:nrow(df1),2] <- 1    ## set delta0
  X[1:nrow(df1),3] <- df1$vax1dose
  X[(nrow(df1)+1):nrow(X),3] <-rep(df1$vax1dose,each=(nrow(df1)-1))
  X[(nrow(df1)+1):nrow(X),4] <- df1$vax1dose  ### TO BE CHANGED
 


  logit_p <- X %*%params
  ### Go back to p: 
  p <- exp(logit_p)/(exp(logit_p) + 1) 

  ### Create matrix (k+1)*(j+1)

  pe_mat <- matrix(NA,nrow=nrow(df1)+1,ncol=nrow(df1)+1)  ## prob escape infection
  pe_mat[1,2:(nrow(df1)+1)]<- 1-p[1:nrow(df1)]
  for(j in 2:(nrow(df1)+1)){
    for(k in 2:(nrow(df1)+1)){
      if(j<k){
        pe_mat[j,k] <- ifelse(j!=k,1,0)*(1-p[k])
      }
      else{pe_mat[j,k] <- ifelse(j!=k,1,0)*(1-p[j])}
    }
  }


#### Create delta's: same as p 
  delta <- pe_mat
  delta[1,2:(nrow(df1)+1)]<- 1

  pe_d <- matrix(NA,nrow=nrow(df1)+1,ncol=nrow(df1)+1)
  pi <- list()
  for(j in 2:(nrow(df1)+1)){
    for(t in 1:max.time.hh){
      for(k in 2:(nrow(df1)+1)){
        delta[k,j] <-ifelse(k!=j,1,0)*ifelse(t > day.infectious[k-1],1,0)*ifelse(t < min(day.exposed[j-1], day.infectious.end[k-1]),1,0)
    #delta[3,2] <-ifelse(t > day.infectious[1],1,0)*ifelse(t < min(day.exposed[2], day.infectious.end[1]),1,0)
      ### define p^d
        pe_d <- colProds(pe_mat^(delta),na.rm = T)
        pi[[t]] <- 1-pe_d
      }
    }
  }

  pi<-do.call("rbind",pi)
  
  ### Likelihood definition:
  ll=0
  ll = ll+product(colProds(pi[,2:ncol(pi)])^(colProds(y,na.rm = T)))*product((1-colProds(pi[,2:ncol(pi)]))^(1-colProds(y,na.rm = T)))
  
  return(-ll)
}





#Note here, we have all time points represented in the df, so the likelihood is very simple--no need to exponentiate stuff  
chain_bin_lik <- function(params,Y,X){
  
  #### Define logit_p = X*params; need  to add as.matrix
  logit_p <- as.vector(as.matrix(X[,1:4]) %*% params) ## Added as.matrix
  
  ### Go back to p (probability  of transmission) with inverse logit: 
  q <- 1 - exp(logit_p)/(exp(logit_p) + 1) 
  
  ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
  #q.spl <- split(q, paste(X$ID, X$hhID, X$t.index)) ## CHECK
  #q.spl <- split(q, paste(ID, hhID, t)) 
  
  pi <- 1 - aggregate(log(q), by=list(X$ID, X$hhID, X$t.index ), FUN=sum)$x
  
  #pi <-  1 - exp(sapply(q.spl,function(x) sum(log(x)))) 
  
  ### Likelihood definition (for the moment no log-lik, so there is just a product over all HH members and time steps):
  ll= sum(dbinom(x=Y,size=1,prob = pi,log = TRUE),na.rm = TRUE)
  return(-ll)
}
##

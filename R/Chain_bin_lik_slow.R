#Note here, we have all time points represented in the df, so the likelihood is very simple--no need to exponentiate stuff  
chain_bin_lik <- function(params,Y,X){
  
  #### Define logit_p = X*params; need  to add as.matrix
  logit_p <- as.vector(as.matrix(X[,1:4]) %*% params) ## Added as.matrix
  
  ### Go back to p (probability  of transmission) with inverse logit: 
  q <- 1 - 1/(1 + exp(-logit_p))
  
  
  ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
  #pi <- 1 - exp(aggregate(log(q), by=list(X$ID, X$hhID, X$t.index ), FUN=sum)$x )
  data_tab <- cbind.data.frame(log.q=log(q),X)
  data_t = data.table(data_tab)
  ans = data_t[,list(A = sum(log.q)), by = 'ID,hhID,t.index']
  ans <- setorder(ans, t.index,ID,hhID) ### To order based on t.index and not ID
  pi= 1- exp(ans$A)
  
  ### Likelihood definition (for the moment no log-lik, so there is just a product over all HH members and time steps):
  ll= sum(dbinom(x=Y,size=1,prob = pi,log = TRUE),na.rm = TRUE)
  return(-ll)
}
##
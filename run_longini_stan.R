#Note, this seems set up to analyze each household individually

###FUNCTIONS
longini_pdf <- function(B,Q,j,k) {
  
  if (j == 0) {
    return(B^k)
  } else {
    t1 <- choose(k,j)*(B^(k-j))*(Q^(j*(k-j)))
    t2 <- 0
    for (i in 0:(j-1)) {
      t2 <- t2 + longini_pdf(B,Q,i,j)
    }
    
    return(t1*(1-t2))
  }
  
}


N <- 10
tl <- rep(0,N+1)
for (i in 0:N) {
  tl[i+1] <- longini_pdf(0.95,0.8,i,N)
}


longini_pdf_zero_trunc <- function(B,Q,j,k) {
  dl <- longini_pdf(B,Q,j,k)
  zl <- longini_pdf(B,Q,0,k)
  
  return(dl/(1-zl))
}

N <- 10
zl <- rep(0,N)
for (i in 1:N) {
  zl[i] <- longini_pdf_zero_trunc(0.95,0.8,i,N)
}

rlongini_zero <- function(B,Q,k,n=1) {
  zl <- rep(0,k)
  for (i in 1:k) {
    zl[i] <- longini_pdf_zero_trunc(B,Q,i,k)
  }
  
  rv <- rmultinom(n,1,zl)
  n <- apply(rv,2,function(x) which(x == 1))
  
  return(n)
}


#############################################

## Sample household sizes
hh_size <- rpois(3000,5) + 1 ## Poisson-distributed number of contacts
rv <- rep(0,3000)
for (i in 1:3000) {
  rv[i] <- rlongini_zero(0.9,0.9,hh_size[i])
}

d <- data.frame(j = rv, k = hh_size)

## Format as input data--each of these represents an individual household
data_in <- list(N = nrow(d),
                j = d$j,
                k = d$k)
##########################################################



## Run model
m <- stan("./R/longini_stan.stan",
          iter = 2000,
          chains = 2,
          sample_file = "testtraces.csv",
          data = data_in)

library(boot)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
source('./R/simulate_data.R')
source('./R/delay_dist_sim.R')

#### Likelihood definition for all HH


###### Set random values for the parameters
alpha0=0.4
delta0=0.7
beta=0.2
kappa=0.5
params <- c(alpha0,delta0,beta,kappa)

#Generate the data and store as a data frame
sim.data.ls <- pblapply(1:10000,gen.hh)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)

##
delay.gen <- function(input_df){
    sim.data.df.spl <- split(input_df, sim.data.df$hhID)
    
    df1.ls <- lapply(sim.data.df.spl, delay_dist_sim)
    
    df1 <- do.call('rbind.data.frame',df1.ls)
    
    #keep just variables we need
    df1b <- df1[, c('ID','hhID','vax1dose','vax2dose','infected','day.exposed','day.infectious.end','day.infectious','max.time.hh')]
    
    ### Create a  copy of df1 and call it df2
    df2 <- df1b
    
    ### Rename all of the variables on df2 as HHID_b, ID_b, Vax_dose1b
    colnames(df2) <- paste(colnames(df2), 'b', sep='_')
    
    ### Merge df1 and df2 by HHID; name this df3 
    df3 <- merge(df1b,df2,by.y='hhID_b' , by.x='hhID')
    
    df3$expand.t <- NA
    df3$expand.t[df3$infected==0] <- df3$max.time.hh[df3$infected==0] #censor uninfected person at max.time.hh
    df3$expand.t[df3$infected==1] <- df3$day.exposed[df3$infected==1] #censor infected person at (latent) day of exposure
    
    df3$rowN <- 1:nrow(df3)
    df.t <- df3[rep(df3$rowN, times=df3$expand.t),] #expands the df to number of time points
    
    #Add an index for each person
    df.t$t.index <- unlist(lapply(df3$expand.t, function(x){ 
      z <- 1:x 
      return(z)
      }))
    
    # keep rows where the time index is within range of the infectious period for person b;
    # or where ID_a==ID_b (community risk)
    if(df.t$ID == df.t$ID_b){ #exogenous
      df.t$keep <-  T   
    }else{ #within-HH
      df.t$keep <- (df.t$t.index >= df.t$day.infectious_b & df.t$t.index <= df.t$day.infectious.end_b ) 
    }
    
    df4 <- df.t[df.t$keep==1,]
    
    #Create elements for design matrix
    df4$delta0 <- 0
    df4$delta0[(df4$ID == df4$ID_b)] <- 1 #exogenous/community risk 
    
    df4$alpha0 <- 0
    df4$alpha0[df4$ID != df4$ID_b ] <- 1
    
    df4$vax1 <- df4$vax1dose
    
    df4$vax2 <- df4$vax1dose_b
    df4$vax2[df4$alpha0==0] <- 0
    
    #HMM is this right? seems like we could double-count Y=1 for certain contact pairs?
    
    df4$infect.at.timet <- df4$infected
    df4$infect.at.timet[df4$t.index < df4$day.exposed] <- 0
    Y.df <- aggregate( df4$infect.at.timet, by=list( 'ID'=df4$ID, 'hhID'=df4$hhID, 'tindex'=df4$t.index ), FUN=mean)
    
    #Design matrix
    X <- df4[c('alpha0','delta0','vax1','vax2')]
    
    Y <-  Y.df$x

  out.list=list('Y'=Y, 'X'=X)    
    return(out.list)
}



#test1 <- df4[,c('hhID', 'ID','ID_b', 't.index','Y')]




##1. Run delay.gen() to create X and Y for a single set of random delay dist value
##2. Run mle() or optim() to estimate parameters
##3 Repeat 1-2 N times
##4 average of all N

#Need Y to represent each HH and time, from 0 to censor time; 
#Check that simulated data truncated to 21 days


#Note here, we have all time points represented in the df, so the likelihood is very simple--no need to exponentiate stuff  
chain_bin_lik <- function(params, ID, hhID,t){

    #### Define logit_p = X*params
    logit_p <- X %*% params
  
  ### Go back to p (probability  of transmission) with inverse logit: 
    q <- 1 - exp(logit_p)/(exp(logit_p) + 1) 
  
    ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
    q.spl <- split(q, paste(ID, hhID, t))
    pi <-  1 - exp(sapply(q.spl, function(x) sum(log(x))   ))
    
  ### Likelihood definition (for the moment no log-lik, so there is just a product over all HH members and time steps):
    ll= sum(dbinom(x=Y,size=1,prob = pi,log = TRUE),na.rm = TRUE)
    return(-ll)
  }
##

#mod.data <- pbreplicate(10000,delay.gen(sim.data.df))

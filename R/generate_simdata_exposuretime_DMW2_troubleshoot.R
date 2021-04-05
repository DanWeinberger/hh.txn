###Simulate  people

#Let's say we expect 20% of households to have a case; 80% without case
#0.8 = p_no_infect^5; p_no_infect=0.956


# 
gen.hh <- function(idN, CPI=(1-0.989), prob.trans.day=(1-0.956), prop.vax1=0.5, prop.vax2=0.5, irr.vax=1, IRR.comm=1){
  
  
  HH.size <- min(1+ rpois(n=1,1.5),2) #cap at 2
  df1 <- as.data.frame(matrix(NA, nrow=HH.size, ncol=2))
  names(df1) <- c('ID', 'hhID')
  df1$hhID <- idN
  df1$ID <- 1:HH.size
  
  # df1$date.test <- as.Date(NA)
  # df1$date.test[1] <- as.Date('2020-12-20')
  # df1$date.test[2:nrow(df1)] <- df1$date.test[1] + 1 + 21*runif(n=(nrow(df1)-1)) 
  # 
  
  df1$vax1dose <- rbinom(nrow(df1), 1,prop.vax1)
  df1$vax2dose <- rbinom(nrow(df1), 1,prop.vax2)

  #Randomly set 1 person as the index
  ##Might need to modify this so that index is unvaccinated person
  
  rand1 <- runif(nrow(df1))
  
  df1$index <- 1*(rand1==min(rand1))

  #These are time distributions that specify for the person if they are infected, how much time to add before infectious, and ow long infectious
  expose.dist= rgamma(length(df1$ID),4,  0.75) #duration latent
  infect.dist= rgamma(length(df1$ID),4, 0.75) #duration infectiousness

  exposed.status <- matrix(NA, nrow=nrow(df1), ncol=100)
  infect.status <- matrix(NA, nrow=nrow(df1), ncol=100)
  n.infect.prev <- matrix(NA, nrow=nrow(df1), ncol=100)
  
  #exposed.status[,1] <- df1$index*rbinom(nrow(df1),1, CPI*IRR.comm^(df1$vax1dose)) #10% chance that a tested index is positive 
  infect.status[,1] <- 0
  n.infect.prev[,1] <- 0
  exposed.status[,1] <- 1 - rbinom(nrow(df1), 1, prob.uninf.day.comm )  #exponent ensure once you are exposed, you stay in that category
  
  prob.infect.day <- prob.trans.day * irr.vax^df1$vax1dose   #prob of being infected per day, per exposure,
  prob.uninfect.day <- 1 - prob.infect.day
  
  prob.uninf.day.comm <- 1 - CPI * irr.vax^df1$vax1dose 
  
  for( i in 2:ncol(infect.status)){
    day.exposed <- apply(exposed.status,1, function(x) which(x==1)[1])
    day.exposed[is.na(day.exposed)] <- 0
    
    day.infect.start <- apply(infect.status,1, function(x) which(x==1)[1])
    day.infect.start[is.na(day.infect.start)] <- 0
    
    n.infect.prev[,i] <- sum(infect.status[,(i-1)]) #how many people in HH were infectious at previous time?

    infect.status[,i]  <- ((i - day.infect.start) <= infect.dist) * ((i - day.exposed ) > expose.dist)  #You are infectious for days in specified range
    
    exposed.status[,i] <- (1-rbinom(nrow(df1), 1, prob.uninf.day.comm*prob.uninfect.day^n.infect.prev[,i] )) ^ (1- exposed.status[,(i-1)]) #exponent ensure once you are exposed, you stay in that category
    
  } 
  
  ##Assume that get PCR on day 3 of being infectious (2 days asymptomatic transmission)
  df1$day.test <- day.infect.start+3
  df1$date.test <- as.Date('2021-01-01') + df1$day.test
  df1$day_index <- as.numeric(df1$date.test - min(df1$date.test)) 
  df1$infected <- apply(infect.status,1,max)
  return(df1)
  
}
  
   
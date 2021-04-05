###Simulate  people

#Let's say we expect 20% of households to have a case; 80% without case
#0.95 = (p_no_infect^(21*5); p_no_infect=0.9995

#0.85 = (p_no_infect^(5); p_no_infect=0.968  duration 5 days


# 
gen.hh <- function(idN, CPI=(1-0.9995), prob.trans.day=(1-0.968), prop.vax1=0.5, prop.vax2=0.5, irr.vax=1, IRR.comm=1){
  
  
  HH.size <- min(1+ rpois(n=1,1.5),5) #cap at 5
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

  exposed.status <- matrix(NA, nrow=nrow(df1), ncol=200)
  infect.status <- matrix(NA, nrow=nrow(df1), ncol=200)
  n.infect.prev <- matrix(NA, nrow=nrow(df1), ncol=200)
  
  prob.infect.day <- prob.trans.day * irr.vax^df1$vax1dose   #prob of being infected per day, per exposure,
  prob.uninfect.day <- 1 - prob.infect.day
  prob.uninf.day.comm <- 1 - CPI * irr.vax^df1$vax1dose 
  
  #exposed.status[,1] <- df1$index*rbinom(nrow(df1),1, CPI*IRR.comm^(df1$vax1dose)) #10% chance that a tested index is positive 
  infect.status[,1] <- 0
  n.infect.prev[,1] <- 0
  exposed.status[,1] <- 1 - rbinom(nrow(df1), 1, prob.uninf.day.comm )  #exponent ensure once you are exposed, you stay in that category
  
  
  for( i in 2:ncol(infect.status)){
    day.exposed <- apply(exposed.status,1, function(x) which(x==1)[1])
    day.exposed[is.na(day.exposed)] <- 0

    day.expose.start <- apply(exposed.status,1, function(x) which(x==1)[1])
    day.expose.start[is.na(day.expose.start)] <- 0
        
    day.infect.start <- (day.expose.start + round(expose.dist))*(i>day.expose.start & (day.expose.start !=0 ))
    day.infect.end <- day.infect.start + round(infect.dist)*(i>day.expose.start & (day.expose.start !=0 ))
    
        n.infect.prev[,i] <- sum(infect.status[,(i-1)]) #how many people in HH were infectious at previous time?

    infect.status[,i]  <- (i >= day.infect.start) * (i <=day.infect.end ) * (day.infect.start>0)  #You are infectious for days in specified range
    
    exposed.status[,i] <- (1-rbinom(nrow(df1), 1, prob.uninf.day.comm*prob.uninfect.day^n.infect.prev[,i] )) ^ (1- exposed.status[,(i-1)]) #exponent ensure once you are exposed, you stay in that category
    
  } 
  
  ##Assume that get PCR on day 1 of being infectious #(2 days asymptomatic transmission)
  if(sum(day.expose.start)>0){
    earliest.expose.hh <- min(day.expose.start[day.expose.start>0])
  } else{
    earliest.expose.hh <- 1
  }
  latest.expose.hh <- min((earliest.expose.hh+20),ncol(exposed.status))
  
  df1$day.test <- day.infect.start 
  df1$date.test <- as.Date('2021-01-01') # + df1$day.test
  df1$day_index <- as.numeric(df1$date.test - min(df1$date.test)) 
  df1$infected <- apply(exposed.status[,earliest.expose.hh:latest.expose.hh, drop=F],1,max) #Only count infections if they were exposed within 21 days of first exposure in HH
  return(df1)
  
}
  
   
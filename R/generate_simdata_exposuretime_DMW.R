set.seed(1234)
library(reshape2)
###Simulate  people
gen.hh <- function(idN, prop.hh.infected=0.5, prop.vax1=0.5, prop.vax2=0.5){
  HH.size <- min(2+ rpois(n=1,1.5),5) #cap at 5
  df1 <- as.data.frame(matrix(NA, nrow=HH.size, ncol=2))
  names(df1) <- c('ID', 'hhID')
  df1$hhID <- idN
  df1$ID <- 1:HH.size
  
  df1$date.test <- as.Date(NA)
  df1$date.test[1] <- as.Date('2020-12-20')
  df1$date.test[2:nrow(df1)] <- df1$date.test[1] + 1 + 21*runif(n=(nrow(df1)-1)) 
  
  ## For vax1dose and vax2dose, I am already accounting for the two weeks lag before vaccine starts to have effect
  # df1$date.vax1dose <- as.Date(NA)
  # df1$date.vax1dose[1] <- as.Date('2020-12-20')
  # df1$date.vax1dose[2:nrow(df1)] <- df1$date.vax1dose[1] + 14 + 5*runif(n=(nrow(df1)-1)) 
  # index_random1 <- rbinom(0:HH.size,HH.size,prob=0.5)
  # df1$date.vax1dose[-index_random1] <- NA
  
  df1$vax1dose <- rbinom(nrow(df1), 1,prop.vax1)
  df1$vax2dose <- rbinom(nrow(df1), 1,prop.vax2)
  

  #### To account for infected people, we need to account for exposure time
    ## Identify the index:
  day.infectious <- rep(NA,length(df1$ID))
  day.exposed <- rep(NA,length(df1$ID))
  day.infectious.end <- rep(NA,length(df1$ID))

 # df1$day_index <- 21 + round(as.numeric(df1$date.test - min(df1$date.test))) #Why adding 21?
  #time index based on the test date
  df1$test_day <- round(as.numeric(df1$date.test - min(df1$date.test)))
  
  infect.dist= rgamma(length(df1$ID),1.33,1/3.46) #duration infectiousness
  expose.dist= rgamma(length(df1$ID),1.33,1/3.46)
  end.inf= rgamma(length(df1$ID),1.33,1/3.46)
  
  for (j in 1:length(df1$ID)){
  #day.matrix=day of test for the person
    df1$day.infectious[j] <- df1$test_day[j] - infect.dist[j]  #infectious prior to test
    df1$day.exposed[j] <- df1$day.infectious[j] - expose.dist[j]  #latent period
    df1$day.infectious.end[j] <- df1$day.infectious[j] + end.inf[j]  #how long infectious after you become infectious
    ###
  }
  
  ### Define index case:
  df1$index.case <- 0
  df1$index.case[df1$day.exposed==min(df1$day.exposed) ] <- 1 
  

  ## Index can be positive or negative, with lower probability if it is vaccinated
  df1$infected <- 0
  for (j in 1:length(df1$ID)){
    if((df1$index.case[j]==1) & (is.na(df1$date.vax2dose[j]))){
      df1$infected[j] <- rbinom(n=1,size=1, p=prob.hh.txn) 
    }
    else if((df1$index.case[j]==1) & (!is.na(df1$date.vax2dose[j]))){
      df1$infected[j] <- rbinom(n=1,size=1, p=prob.hh.txn*0.2) 
    }
  }
  
  #
  df1$infect.index <- round(df1$day.infectious - min(df1$day.infectious)) + 1
  
  for( i in 1:max(df1$infect.index))
  
  
  
  
  #### For now, only HH members with a positive index within the HH can become infected (no risk of becoming infected)
  #### from the community
  for (j in 1:length(df1$ID)){
    if(df1$infected[df1$index.case==1]==1){
      #print(paste0('part: ',df1$infected[df1$index.case==1]==1))
      if(is.na(df1$date.vax2dose[j])){
      df1$infected[j] <- rbinom(n=1,size=1, p=prob.hh.txn) 
      }
      else if(!is.na(df1$date.vax2dose[j])){
      df1$infected[j] <- rbinom(n=1,size=1, p=prob.hh.txn*0.5) 
    }
    }
  }
  return(df1)
}

#Generate the data and store as a data frame
hh_list <- lapply(1:100,gen.hh )
hh_df <- do.call('rbind.data.frame', hh_list)


#Count time since first test
##Split by HH to do this
#first case has day index set to 21 to allow for infection/latent period to happen earlier Tt1-20
hh_df.spl <- split(hh_df, hh_df$hhID)
hh_df.spl <- lapply(hh_df.spl, function(x){
  x$day_index <- 21 + round(as.numeric(x$date.test - min(x$date.test)))
  x$day_index[x$infected==0] <- 42
  return(x)
})
hh_df <- do.call('rbind.data.frame', hh_df.spl)

##TO ADD TO THIS: Day when hospitalized relative to date of test--then they no longer contribute


hh_df.m <- melt(hh_df[,c('day_index','hhID','ID','infected','vax1dose','vax2dose')], id.vars=c('hhID','ID'))
hh_df.c <- acast(hh_df.m, hhID~ ID ~variable )

day.matrix <- hh_df.c[,,'day_index']
infected_matrix <- hh_df.c[,,'infected']
vax <- hh_df.c[,,'vax1dose']
#vax1dose <- hh_df.c[,,'vax1dose']
#vax2dose <- hh_df.c[,,'vax2dose']
N.hh.members <- apply(infected_matrix,1,function(x) sum(!is.na(x)))
N.HH <- nrow(infected_matrix)

N.cases.hh <- apply(infected_matrix,1,sum, na.rm=T)


set.seed(1234)
library(reshape2)
###Simulate  people
gen.hh <- function(idN, prob.hh.txn=0.5){
  HH.size <- min(2+ rpois(n=1,1.5),5) #cap at 5
  df1 <- as.data.frame(matrix(NA, nrow=HH.size, ncol=2))
  names(df1) <- c('ID', 'hhID')
  df1$hhID <- idN
  df1$ID <- 1:HH.size
  first.date <- '2020-04-01'
  df1$date.test <- as.Date(NA)
  df1$date.test[1] <- as.Date('2020-12-20')
  df1$date.test[2:nrow(df1)] <- df1$date.test[1] + 1 + 21*runif(n=(nrow(df1)-1)) 
  
  ## For vax1dose and vax2dose, I am already accounting for the two weeks lag before vaccine starts to have effect
  df1$date.vax1dose <- as.Date(NA)
  df1$date.vax1dose[1] <- as.Date('2020-12-20')
  df1$date.vax1dose[2:nrow(df1)] <- df1$date.vax1dose[1] + 14 + 5*runif(n=(nrow(df1)-1)) 
  index_random1 <- sample(0:HH.size,HH.size)
  df1$date.vax1dose[-index_random1] <- NA
  
  df1$vax1dose <- 0
  df1$vax1dose[!(is.na(df1$date.vax1dose))]<-1
  
  df1$date.vax2dose <- as.Date(NA)
  df1$date.vax2dose[1:nrow(df1)] <- df1$date.vax1dose[1] + 21 +14+ 2*runif(n=(nrow(df1))) 
  df1$date.vax2dose[is.na(df1$date.vax1dose)]<-NA
  index_random2 <- sample(0:HH.size,HH.size)
  df1$date.vax2dose[-index_random2] <- NA
  df1$vax2dose <- 0
  df1$vax2dose[!(is.na(df1$date.vax2dose))]<-1
  

  df1$infected <- NA
  df1$infected[1] <- 1
  df1$infected[2:nrow(df1)] <- rbinom(n=(nrow(df1)-1),size=1, p=prob.hh.txn)
  
  
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
vax <- hh_df.c[,,'vax2dose']
#vax1dose <- hh_df.c[,,'vax1dose']
#vax2dose <- hh_df.c[,,'vax2dose']
N.hh.members <- apply(infected_matrix,1,function(x) sum(!is.na(x)))
N.HH <- nrow(infected_matrix)

N.cases.hh <- apply(infected_matrix,1,sum, na.rm=T)
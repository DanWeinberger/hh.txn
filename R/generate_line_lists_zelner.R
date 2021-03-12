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
  df1$date.test[1] <- as.Date('2021-01-01')
  df1$date.test[2:nrow(df1)] <- df1$date.test[1] + 1 + 21*runif(n=(nrow(df1)-1)) 
  
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
#first case has day index set to 21 to allow for infection/latent period to happen earlier Tt1-20)
hh_df.spl <- split(hh_df, hh_df$hhID)
hh_df.spl <- lapply(hh_df.spl, function(x){
  x$day_index <- 21 + round(as.numeric(x$date.test - min(x$date.test)))
  x$day_index[x$infected==0] <- 42
  return(x)
})
hh_df <- do.call('rbind.data.frame', hh_df.spl)

##TO ADD TO THIS: Day when hospitalized relative to date of test--then they no longer contribute


hh_df.m <- melt(hh_df[,c('day_index','hhID','ID','infected')], id.vars=c('hhID','ID'))
hh_df.c <- acast(hh_df.m, hhID~ ID ~variable )

day.matrix <- hh_df.c[,,'day_index']
infected_matrix <- hh_df.c[,,'infected']
N.hh.members <- apply(infected_matrix,1,function(x) sum(!is.na(x)))
N.HH <- nrow(infected_matrix)

N.cases.hh <- apply(infected_matrix,1,sum, na.rm=T)



# 
# 
# #####################
# #distribution for 
# for(t in 1:42){ #how can we do this without looping over all t?
#   
#   for(i in 1:N.HH){
#     for(j in 1:N.hh.members[i]){
#       
#       day.infectious[i,j] <- day.matrix[i,j] - rgamma(1, shape=0.75, scale=3) #infectious prior to test
#       day.exposed[i,j] <- day.infectious[i,j] - rgamma(1, shape=0.75, scale=3) #latent period
#       day.infectious.end[i,j] <- day.infectious[i,j] + rgamma(1, shape=0.75, scale=7) #how long infectious after 
#       
#       
#       I[i,j,t] <- (t >= day.infectious[i,j] & t<day.infectious.end[i,j] & infected_matrix[i,j]==1)
#       E[i,j,t] <- (t >= day.exposed[i,j] & t<day.infectious[i,j] & infected_matrix[i,j]==1)
#       S[i,j,t] <- (infected_matrix[i,j]==0 | (infected_matrix[i,j]==1 & t<=day.exposed[i,j]) )
#     }
#     I_sum[i,t] <- sum(I[i,,t])
#     E_sum[i,t] <- sum(E[i,,t])
#     S_sum[i,t] <- sum(S[i,,t])
#     
#   }
# }
# 
# 
# #JAGS code
# model{
#   for(i in 1:N.HH){
#     for(j in 1:N.hh.members[i]){
#       
#       day.infectious[i,j] <- day.matrix[i,j] - infect.dist[i,j] #infectious prior to test
#       day.exposed[i,j] <- day.infectious[i,j] - expose.dist[i,j] #latent period
#       day.infectious.end[i,j] <- day.infectious[i,j] + end.inf[i,j] #how long infectious after 
#       
#       infect.dist[i,j] ~ dgamma(sh1, ra1)
#       expose.dist[i,j] ~dgamma(sh2, ra2)
#       end.inf[i,j] ~dgamma(sh3,ra3)
#       
# for(t in 1:42){ #how can we do this without looping over all t?
#         
#       I[i,j,t] <- (t >= day.infectious[i,j] & t<day.infectious.end[i,j] & infected_matrix[i,j]==1)
#       E[i,j,t] <- (t >= day.exposed[i,j] & t<day.infectious[i,j] & infected_matrix[i,j]==1)
#       S[i,j,t] <- (infected_matrix[i,j]==0 | (infected_matrix[i,j]==1 & t<=day.exposed[i,j]) )
# }
#     }
# for(t in 1:42){ #how can we do this without looping over all t?
#   I_sum[i,t] <- sum(I[i,,t])
#   E_sum[i,t] <- sum(E[i,,t])
#   S_sum[i,t] <- sum(S[i,,t])
# }
#   }
# 
#    # Hyperpriors for the latent distributions
#     # parameterized by mode (m) and standard deviation (sd):
#     sh1 <- 1 + m1 * ra1
#     ra1 <- ( m1 + sqrt( m1^2 + 4*sd1^2 ) ) / ( 2 * sd1^2 )
#     m1 ~ dunif(2,6)
#     sd1 ~ dunif(2,3) 
#     
#     sh2 <- 1 + m2 * ra2
#     ra2 <- ( m2 + sqrt( m2^2 + 4*sd2^2 ) ) / ( 2 * sd2^2 )
#     m2 ~ dunif(2,6)
#     sd2 ~ dunif(2,3) 
#     
#     ra3 <- ( m3 + sqrt( m3^2 + 4*sd3^2 ) ) / ( 2 * sd3^2 )
#     sh3 <- 1 + m3 * ra3
#     m3 ~ dunif(2,6) #days
#     sd3 ~ dunif(2,3) #SD on days
#     
# }    

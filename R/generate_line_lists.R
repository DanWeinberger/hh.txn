library(reshape2)



 df1 <- as.data.frame(matrix(NA, nrow=3, ncol=5))
 names(df1) <- c('contactID', 'indexID', 'date.index.sympt','date.index.isolated', 'infected')
 df1$contactID <- 1
 df1$indexID <- 1:3
 df1$date.index.sympt <- as.Date(c('2020-03-19','2020-03-23', '2020-03-25'))
 df1$date.index.isolated <- as.Date(c('2020-03-22','2020-03-27', '2020-03-27'))
 
 df1$fu.time <- as.numeric(df1$date.index.isolated - df1$date.index.sympt) +1
 df1$infected <-0
 df1$tf <- max(df1$date.index.isolated) - min(df1$date.index.sympt) +1
 df1$date.df <- as.Date(NA) 
 
 
 
#Contact2 infected 
 df2 <- as.data.frame(matrix(NA, nrow=2, ncol=5))
 names(df2) <- c('contactID', 'indexID', 'date.index.sympt','date.index.isolated', 'infected')
 df2$contactID <- 2
 df2$indexID <- 1:2
 df2$date.index.sympt <- as.Date(c('2020-03-19','2020-03-23'))
 df2$date.index.isolated <- as.Date(c('2020-03-22','2020-03-25'))
 df2$date.df <- as.Date(c('2020-03-25'))
 
  
 df2$fu.time <- as.numeric(df2$date.index.isolated - df2$date.index.sympt) + 1
 df2$infected <-1
 df2$tf <- max(df2$date.index.isolated) - min(df2$date.index.sympt) +1
 
  #Combine data frames from all people
 df.all <- rbind.data.frame(df1, df2)

 timepoints <- max(df.all$tf)
 N.index <- max(df.all$indexID)
 N.contact <- max(df.all$contactID) 
 
 T <- array(0, dim=c(N.contact,timepoints, N.index))
 
 for(i in 1:N.contact){
   ds <- df.all[df.all$contactID==i,]
   
   #For uninfected people
   if(max(ds$infected==0)){
     
    for(j in 1:max(ds$indexID) ){
     T[i, 1:ds$fu.time[j],j] <- 1:ds$fu.time[j]
    }
     
     #For infected people
   }else{
     ds$start.obs.index <- ds$tf - as.numeric(ds$date.df - ds$date.index.sympt )
     ds$end.obs.index <- ds$start.obs.index + ds$fu.time -1
     for(j in 1:max(ds$indexID) ){
     T[i,ds$start.obs.index[j]:ds$end.obs.index[j] ,j] <- 1:ds$fu.time[j]
      }
   }
   
   
   
 }
 
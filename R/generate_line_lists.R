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
 
 
 
 ###Simulate More people
 gen.uninfected.contact <- function(idN){
    N.indexes <- min(1+ rpois(n=1,1.5),5) #cap at 5
    df1 <- as.data.frame(matrix(NA, nrow=N.indexes, ncol=2))
    names(df1) <- c('contactID', 'indexID')
    df1$contactID <- idN
    df1$indexID <- 1:N.indexes
    first.date <- '2020-04-01'
    df1$date.index.test <- sample(seq(as.Date(first.date), as.Date(first.date)+14, by="day"), N.indexes)
    df1$date.contact.test <- df1$date.index.test + 21 #max 21 day follow up
   
    df1$fu.time <- as.numeric(df1$date.contact.test - df1$date.index.test) +1
    df1$infected <-0
    df1$date.df <- as.Date(NA) 
    return(df1)
 }
 
 gen.infected.contact <- function(idN){
    N.indexes <- min(1+ rpois(n=1,1.5),5) #cap at 5
    df1 <- as.data.frame(matrix(NA, nrow=N.indexes, ncol=2))
    names(df1) <- c('contactID', 'indexID')
    df1$contactID <- idN
    df1$indexID <- 1:N.indexes
    first.date <- '2020-04-01'
    df1$date.index.test <- sample(seq(as.Date(first.date), as.Date(first.date)+14, by="day"), N.indexes)
    df1$date.contact.test <- df1$date.index.test + rpois(n=1,5) +1
    df1$fu.time <- as.numeric(df1$date.contact.test - df1$date.index.test) +1
    df1$infected <-1
    df1$date.df <- as.Date(df1$date.contact.test)
    return(df1)
 }
 
 #Generate the data
 uninfected_list <- lapply(1:1200,gen.uninfected.contact )
 infected_list <- lapply(1201:2400,gen.infected.contact)

 #Combine the data frames into 1
 all_subjects <- do.call('rbind.data.frame', c(uninfected_list,infected_list)) 

 timepoints <- max(all_subjects$fu.time)
 N.index <- max(all_subjects$indexID)
 N.contact <- max(all_subjects$contactID) 
 
 T <- array(0, dim=c(N.contact,timepoints, N.index))
 vax1.index <- array(0, dim=c(N.contact,timepoints, N.index))
 vax2.index <- array(0, dim=c(N.contact,timepoints, N.index))
 
 vax1.contact <- array(0, dim=c(N.contact,timepoints))
 vax2.contact <- array(0, dim=c(N.contact,timepoints))
 
 for(i in 1:N.contact){
    ds <- all_subjects[all_subjects$contactID==i,]
    
    #For uninfected people
    if(max(ds$infected==0)){
       
       for(j in 1:max(ds$indexID) ){
          T[i, 1:ds$fu.time[j],j] <- 1:ds$fu.time[j]
       }
       
       #For infected people
    }else{
       ds$start.obs.index <- ds$fu.time - as.numeric(ds$date.df - ds$date.index.test )
       ds$end.obs.index <- ds$start.obs.index + ds$fu.time -1
       for(j in 1:max(ds$indexID) ){
          T[i,ds$start.obs.index[j]:ds$end.obs.index[j] ,j] <- 1:ds$fu.time[j]
       }
    }
    
    
    
 }
 
library(lme4)
delay.gen.simplified <- function(input_df){
  sim.data.df.spl <- split(input_df, input_df$hhID)
  
  df1.ls <- lapply(sim.data.df.spl, delay_dist_sim)
  df1 <- do.call('rbind.data.frame',df1.ls)
 # df1 <- df1[1:4,]
 # df1$ID <- c(1,2,3,4)
 # df1$hhID <- c(1,1,1,1)
 # df1$infected <- c(1,1,0,0)
 # df1$day.infectious  <- c(3,6,4,9)
 # df1$day.exposed <- c(2,4,3,6)
 # df1$day.infectious.end <- c(6,9,8,11)
 # df1$max.time.hh <- c(11,11,11,11)
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
  
  df.t$keep <-  T *ifelse(df.t$ID == df.t$ID_b,1,0)  | (df.t$t.index >= df.t$day.infectious_b & df.t$t.index <= df.t$day.infectious.end_b )*ifelse(df.t$infected_b==1,1,0)
  
  
  df4 <- df.t[df.t$keep==1,]
  
  #Create elements for design matrix
  df4$delta0 <- 0
  df4$delta0[(df4$ID == df4$ID_b)] <- 1 #exogenous/community risk 
  
  df4$alpha0 <- 0
  df4$alpha0[df4$ID != df4$ID_b ] <- 1
  
  df4$vax1 <- df4$vax1dose
  
  df4$vax2 <- df4$vax1dose_b
  df4$vax2[df4$delta0==1] <- 0 #person B's vaccine status doesn't influence whether person A is infected
  
  df4$infect.at.timet <- df4$infected
  df4$infect.at.timet[df4$t.index < df4$day.exposed] <- 0
  
  
  df4$atrisk <- 0
  df4$atrisk[(df4$t.index >= df4$day.infectious_b & df4$t.index <= df4$day.infectious.end_b )  & (df4$infected_b==1)] <- 1
  
  df4$m <- 0
  df4$m[(df4$atrisk==1) &  (df4$vax1dose_b!=1)] <- 1
  
  df4$z <- 0
  df4$z[(df4$atrisk==1) &  (df4$vax1dose_b==1)] <- 1
  
  
  data_tab.X <- cbind.data.frame(df4[,c('m','ID','hhID','t.index','vax1dose')])
  data_t.X = data.table(data_tab.X)
  X.df = data_t.X[,list(A = sum(m)), by = 'ID,hhID,t.index,vax1dose']
  colnames(X.df)[ncol(X.df)]  <- 'm'
  
  data_tab.X.z <- cbind.data.frame(df4[,c('z','ID','hhID','t.index','vax1dose')])
  data_t.X.z = data.table(data_tab.X.z)
  X.df.z = data_t.X.z[,list(A = sum(z)), by = 'ID,hhID,t.index,vax1dose']
  colnames(X.df.z)[ncol(X.df.z)]  <- 'z'
  
  
  X.df$z <-  X.df.z$z
  
  #Y.df <- aggregate( df4$infect.at.timet, by=list( 'ID'=df4$ID, 'hhID'=df4$hhID, 'tindex'=df4$t.index ), FUN=mean)
  data_tab.Y <- cbind.data.frame(df4[,c('infect.at.timet','ID','hhID','t.index','vax1dose')])
  data_t.Y = data.table(data_tab.Y)
  Y.df = data_t.Y[,list(A = mean(infect.at.timet)), by = 'ID,hhID,t.index,vax1dose']
  #Y.df = data_t[,list(A = mean(infect.at.timet)), by = 'ID,hhID,t.index,vax1dose,infected_b']
  
  Y.df <- setorder(Y.df, hhID, ID, t.index,vax1dose) ### To order based on t.index and not ID
  

  
  #Y.df$m <- 0
  #Y.df$z <- 0
  #df4.agg.z <- df4[(df4$vax1dose_b==1) &(df4$t.index>=df4$day.infectious_b) & (df4$t.index<=df4$day.infectious.end_b) & (df4$infected_b==1),]
  #df4.agg.m <- df4[(df4$vax1dose_b!=1) &(df4$t.index>=df4$day.infectious_b) & (df4$t.index<=df4$day.infectious.end_b) & (df4$infected_b==1),]
  
  #df4.agg.m <- df4[(df4$vax1dose_b!=1),]
  #agg.m <- aggregate(df4.agg.m$ID_b, by=list('ID'=df4.agg.m$ID, 'hhID'=df4.agg.m$hhID, 'tindex'=df4.agg.m$t.index ), FUN=length)
  
 
# 
# 
#   if(nrow(df4.agg.m)!=0){
#     agg.m <- aggregate(df4.agg.m$ID_b, by=list('ID'=df4.agg.m$ID, 'hhID'=df4.agg.m$hhID, 'tindex'=df4.agg.m$t.index ), FUN=length)
#     #agg2 <- 0
#     Y.df[(Y.df$ID %in% agg.m$ID) & (Y.df$hhID %in% agg.m$hhID) &  (Y.df$t.index %in% agg.m$tindex) & (Y.df$t.index>=Y.df$day.infectious_b) & (Y.df$t.index<= Y.df$day.infectious.end_b),]$m <- agg.m$x
#   }
#   
#   if(nrow(df4.agg.z)!=0){
#     agg.z <- aggregate(df4.agg.z$ID_b, by=list('ID'=df4.agg.z$ID, 'hhID'=df4.agg.z$hhID, 'tindex'=df4.agg.z$t.index ), FUN=length)
#     Y.df[(Y.df$ID %in% agg.z$ID) & (Y.df$hhID %in% agg.z$hhID) &  (Y.df$t.index %in% agg.z$tindex) & (Y.df$t.index>=Y.df$day.infectious_b) & (Y.df$t.index<= Y.df$day.infectious.end_b),]$z <- agg.z$x
#     
#   }

  Y <-  Y.df$A
  
  X.df$Y <- Y
  
  out.list=list( 'X.df'=X.df)  
  return(out.list)
}






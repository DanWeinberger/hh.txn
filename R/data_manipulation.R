delay.gen <- function(input_df){
  sim.data.df.spl <- split(input_df, input_df$hhID)
  
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
  #if(df.t$ID == df.t$ID_b){ #exogenous ## This doesn't work
  # df.t$keep <-  T   
  #}else{ #within-HH
  #  df.t$keep <- (df.t$t.index >= df.t$day.infectious_b & df.t$t.index <= df.t$day.infectious.end_b ) 
  #}
  
  ### CHECK: 
  df.t$keep <-  T *ifelse(df.t$ID == df.t$ID_b,1,0)  | (df.t$t.index >= df.t$day.infectious_b & df.t$t.index <= df.t$day.infectious.end_b )*ifelse(df.t$infected_b==1,1,0)
  
  
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
  
  #Y.df <- aggregate( df4$infect.at.timet, by=list( 'ID'=df4$ID, 'hhID'=df4$hhID, 'tindex'=df4$t.index ), FUN=mean)
  data_tab <- cbind.data.frame(df4[,c('infect.at.timet','ID','hhID','t.index')])
  data_t = data.table(data_tab)
  Y.df = data_t[,list(A = mean(infect.at.timet)), by = 'ID,hhID,t.index']
  Y.df <- setorder(Y.df, t.index,ID,hhID) ### To order based on t.index and not ID
  
  
  #Design matrix 
  X <- df4[c('alpha0','delta0','vax1','vax2','ID','hhID','t.index')]  ###CHECK if OK   
  
  Y <-  Y.df$A
  
  out.list=list('Y'=Y, 'X'=X)  
  return(out.list)
}



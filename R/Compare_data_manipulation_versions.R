library(boot)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
library(stats4)

source('./R/simulate_data.R')
source('./R/delay_dist_sim.R')
source('./R/Chain_bin_lik.R')
source('./R/data_manipulation_test.R')
source('./R/data_manipulation.R')



#Generate the synthetic data and store as a data frame
N.HH <- 5000
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), #Increase CPI to have more cases, 
                        prob.trans.day=(1-0.968),
                        irr.vax1=0.2,irr.vax2=1)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)

### Compare different versions for data manipulation

### Original data manipulation file
ptm <- proc.time()
LatentData <-  delay.gen(input_df = sim.data.df)
Y <- LatentData$Y
X <- LatentData$X
proc.time() - ptm

### Before merging df1b and df2, convert them to data.tables
ptm <- proc.time()
LatentData <-  delay.gen.dtable(input_df = sim.data.df)
Y <- LatentData$Y
X <- LatentData$X
proc.time() - ptm


### Loop through HHs using lapply

ptm <- proc.time()
LatentDataFast <- lapply(sim.data.ls,delay.gen.singleHH)
X <- lapply(LatentDataFast, function(x) x<- x$X)
Y <- lapply(LatentDataFast, function(x) x<- x$Y)
X<-do.call('rbind.data.frame',X)
Y<-do.call('rbind.data.frame',Y)
proc.time() - ptm


#### Explicit for loop through HHs
ptm <- proc.time()
X.list <- list()
Y.list <- list()
for(i in 1:length(sim.data.ls)){
  LatentData <- delay.gen.singleHH(sim.data.ls[[i]])
  X.list[[i]] <- LatentData$X
  Y.list[[i]] <- LatentData$Y
}
X<-do.call('rbind.data.frame',X.list)
Y<-do.call('rbind.data.frame',Y.list)
proc.time() - ptm


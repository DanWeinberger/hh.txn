library(boot)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
library(stats4)
library(lme4)


source('./R/simulate_data.R')
source('./R/delay_dist_sim.R')
source('./R/data_manipulation_simple_model.R')


################ Define glmer function
#Generate the synthetic data and store as a data frame
N.HH <- 50
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), #Increase CPI to have more cases, 
                        prob.trans.day=(1-0.968),
                        irr.vax1=0.2,irr.vax2=1)

#This is like the data we would get from KSM
sim.data.df <- do.call('rbind.data.frame', sim.data.ls)

input_df <- sim.data.df
LatentData <- delay.gen.simplified(sim.data.df)
Y <- LatentData$Y
X <- LatentData$X.df

X <- data.frame(X)
outcome_name <- 'Y'

formula <- as.formula(paste0(outcome_name, '~',
                             'vax1dose+m+z+(1|hhID)'))
mod1 <- glmer(formula = formula, data=X, family=binomial(link="cloglog"),verbose=TRUE)


summary(mod1)




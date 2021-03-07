#################################
gen_sim_data <- function(popN=2400,IRR.vax=1){

set.seed(123)
nprob=IRR.vax*array(c(0.01, 0.01,0.02,0.02,0.04,0.04,0.01,0.01),dim=c(1,8)) #probability of transmission from "normal" index case on day 1-8 of illness
hprob=IRR.vax*0.06*array(1, dim=c(1,8)) #probability of transmission from "highly infectious" index case on day 1-8
popsize=popN #number of contacts exposed
truth=nprob*2/3+hprob/3 #"true" probability of transmission (population avg)

g=dgamma(0.5:29.5,1.33,1/3.46) # %gamma probability distribution of incubation period
g=rev(g) #flip vector (should correspond to flipud in matlab)

hpny <- array(NA,dim=c(1,8))
npny <- array(NA,dim=c(1,8))

hpny[1]=1-hprob[1] #probability contact exposed to highly infectious index case is not infected on day 1
npny[1]=1-nprob[1]

for (i in 2:8){
  hpny[i]=hpny[i-1]*(1-hprob[i]) #probability contact in not infected on day d=2:8
  npny[i]=npny[i-1]*(1-nprob[i])
}

hPy=1-hpny
nPy=1-npny

hcase=(popsize/24)*hPy #Number of contacts exposed to highly infectious index case who are infected on day d=1:8
ncase=(popsize/12)*nPy

totcase=round(hcase+ncase) #Total number of contacts infected on day d=1:8

S=array(0, dim=c(popsize,8)) #Matrix for the number of index cases to whom contact i is exposed who are on day j of illness 
y=array(0, dim=c(popsize,1)) #Infection status of contact
for (s in 1:8){
  S[(popsize/8*(s-1)+1):popsize,s]=1
  y[(popsize/8*(s-1)+1):(popsize/8*(s-1)+totcase[s]),1]=1
}

P=cbind(y,S)
P = P[order(P[,1],decreasing=F),]
#P=sortrows(P,[1]) 
P= apply(P,2,rev) #Infected contacts must be listed first
S=P[,2:9]
nY=1-P[,1] #Probability contact is not infected

T=array(0,dim=c(popsize,40)) #Matrix for number of index cases to whom contact i is exposed at time j (infected contacts become symptomatic at t=30 for uninfected contacts, exposure begins at t=1) 
### CHECK! Matrix should be of size (popsize,30)
inc <- rep(NA, sum(totcase))
inc <- round(rgamma(sum(totcase),1.33,1/3.46)) #incubation period chosen at random (from gamma distribution) for each infected contact

stepc=array(0,dim=c(8,1))
stepc[1]=totcase[1]
for (i in 2:8){
stepc[i]=totcase[i]-totcase[i-1] #number on contacts who were infected on day d=1:8
}
####

i=0
for (j in 1:8) {#contact infected on day j
  for (k in (i+1):(i+(9-j)*stepc[j])){ #infected contact #
    for (m in 1:j){ #accounting for all prior exposure-days
      T[k,30-inc[k]-j+m]=m #contact symptom onset = day 30, incubation period >= 0
    }
  }
  for (n in (j+1):8) { #accounting for all subsequent exposure-days
    for (k in (i+(n-j)*stepc[j]+1):(i+(9-j)*stepc[j])) {#infected contact #
      T[k,30-inc[k]-j+n]=n
    }
  }
  i=i+(9-j)*stepc[j]
}



for (i in (1+sum(totcase)):popsize){ #uninfected contacts
  for (j in 2:9){
    T[i,j]=T[i,j-1]+1
    if (T[i,j] >=sum(S[i,])){
      break
    }
  }
}      
T=T[,1:38]

#These are the other two variables you need to run the WinBugs code ("known" incubation periods and infection status of those contacts)
incub <- rep(NA,81)
incub <- round(rgamma(81,1.33,1/3.46))

obs=array(1, dim=c(81,1))



data.sim<-list('T'=T, 'nY'=nY,'obs'=obs,'incub'=incub)
return(data.sim)
}

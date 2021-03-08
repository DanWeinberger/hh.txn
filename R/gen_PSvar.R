gen_PSvar <- function(Cstrata,Istrata,Ena)
  {
  
  x <- array(1, dim=c(length(Ena[,1]),1))
  
  for (i in 2:length(x)){
    if (Ena[i,2]==Ena[i-1,2]){
      x[i,1]=0
      }
  }
  
  y=which(x!=0)
  
  Inum <- array(NA,dim=c(length(y),1))
  for (i in 1:length(y)-1){
    Inum[i,1]=y[i+1]-y[i];    #number of index cases contact i is exposed to
  }
  
  Inum = replace(Inum, length(Inum), 1)
  kmax=max(Inum)
  
  Ri=Ena[y,2:3]  #contact ID and infection status
  R=array(0,dim=c(length(Ri[,1]),3,kmax))
    for (i in 1:length(Ri[,1])){
      for (k in 1:kmax){
        if (k<=Inum[i,1]){
          R[i,1,k]=Istrata[y[i]+k-1,1]+1 #index case strata
          R[i,2,k]=Ena[y[i]+k-1,6] #duration of exposure
          R[i,3,k]=Ena[y[i]+k-1,7]-Ena[y[i]+k-1,4] #serial interval (max incubation period)
        }
      }
  } 
  
  
  maxinc=39
  T=array(0,dim=c(length(y),maxinc,kmax))
    for (i in 1:nrow(T)){
      for (k in 1:kmax){
        for (t in 2:maxinc){
          if (R[i,3,k]>=(maxinc-t) || R[i,3,k]<0){ #start tracking exposure at maxinc - serial interval (or t=2 for uninfected)
            T[i,t,k]=R[i,1,k] #index case strata (1=female, 2=male)
          }
          if (t==(maxinc-R[i,3,k]+R[i,2,k]) || (R[i,3,k]<0 && t==2+R[i,2,k])){ #stop tracking exposure when t = maxinc - s.i. + dur
            break                                                        #  or t=2+dur for uninfected
          } 
        }
      }
    }

  return(T)
  
}



### Small test for gen_PS function: 

Ena <- read.csv("data/SARSfakedata.csv")
Istrata <- as.matrix(Ena[,'X.index.sex.'],dim=c(length(Ena[,1]),1))
#Cstrata <- Ena[,c('X.contact.ID.','X.contact.inf.status.','X.contact.sex.')]
Cstrata <- Ena[,c('X.contact.ID.','X.contact.age.','X.contact.sex.')] ## ask Ginny

T <- gen_PSvar(Cstrata,Istrata,Ena)





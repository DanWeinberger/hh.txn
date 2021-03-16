
for(i in 1:N.hh){
  index.case[i]  ~dcat(prob.index[i,])
  prob.index[i,] <- Test.Pos[i,] * 1/sum(Test.Pos[i,])
  
  for(j in 1:N.hh.members[i]){
    
    y[i,j] ~ dbern(q[i,j]) #y=UNinfected matrix
    
    ##NOTE THESE ARE PROBABLY NOT RIGHT--PROB NEED TO FLIP SOME OF THEM AROUND
    prob_no_inf_uninf[i,j] <- exp(sum(log_prob_no_inf_t[i,1:day.matrix[i,j]])) #P no infections over all time intervals
    
    prob_no_inf_inf_person[i,j] <- exp(sum(log_prob_no_inf_t[i,1:((day.matrix[i,j] - 1)) ])) #P no infections prior to infection
    
    prob_inf[i,j] <- infected_matrix[i,j] * (beta[i,j] * I[i,j,day.matrix[i,j]] + alpha) #for infected people only
    
    q[i,j] <- (prob_no_inf_inf_person[i,j]) * (1 - prob_inf[i,j]) * infected_matrix[i,j] + #prob for infection at time t and not before
      (prob_no_inf_uninf[i,j]) * (1- infected_matrix[i,j] ) + 1e-6 #prob for uninfected peopel    
    
    Test.Pos[i,j] <- 1 - y[i,j] #does the person test pos at some point?
    
    S[i,j,1] <- 1 * step(abs(index.case[i] - i ) -0.5) #For index case, S=0
    E[i,j,1] <- 0
    I[i,j,1] <- 1,1] <- 1 * (1- step(abs(index.case[i] - i ) -0.5)) #for index case I=1
    R[i,j,1] <- 0

for(t in 2:(tmax[i]+1)){
  #Define current state
  
  p.S[i,j,t] <- S[i,j,(t-1)] * (1-(beta[i,j] * (sum_I[i,(t-1)]-I[i,j,(t-1)]) + alpha[i,j])) #Prob remain in S;
  p.E[i,j,t] <- E[i,j,(t-1)] * (1-epsilon[i,j]) #Prob stay in E
  p.I[i,j,t] <- I[i,j,(t-1)] * (1-delta[i,j])   #Prob stay infected
  
  log_prob_no_inf_t[i,j,t] <- log( 1 - (beta[i,j] * (sum_I[i,(t-1)]-I[i,j,(t-1)])[i,t] + alpha[i,j]) )
  
  
  S.Draw[i,j,t] ~ dbern(p.S[i,j,t])
  E.Draw[i,j,t] ~ dbern(p.E[i,j,t])
  I.Draw[i,j,t] ~ dbern(p.I[i,j,t])
  
  S[i,j,t] <- S[i,j,(t-1)] * S.Draw[i,j,t] #previously S, draw a 0 or 1
  
  E[i,j,t] <- S[i,j,(t-1)]*(1-S.Draw[i,j,t]) + #previously S, current draw S=0
    E[i,j,(t-1)]*E.Draw[i,j,t]  #Previously E, stay E?
  
  I[i,j,t] <- E[i,j,(t-1)]*(1-E.Draw[i,j,t]) + #previously E, current draw E=0
    I[i,j,(t-1)]*I.Draw[i,j,t]  #Previously I, stay I?
  
  R[i,j,t] <- R[i,j,t-1] + I[i,j,(t-1)]*(1-I.Draw[i,j,t]) #previously I, current draw I=0
  
}
  }
}
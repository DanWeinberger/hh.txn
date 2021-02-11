
N.contacts <- 1200
N.infected <- 300
N.times <- 20
nY <- rep(0, N.contacts)
nY[1:N.infected] <- 1

#j=index cases, i= contacts, t= times
#t is defined base don when household starts
T_ijt <- array(NA, dim=c(N.contacts,5,20 ))

#Contact1 has 1 index case and 20 days of follow up
T_ijt[1,1,1:20] <- 1:20

#Contact 2 has 2 index cases and 15 days of follow up, for 2nd index, starts on day 6
T_ijt[2,1,1:15] <- 1:15
T_ijt[2,2,1:10] <- 6:15

#Contact 3 has 3 index cases and 12 days of follow up
T_ijt[3,1,1:12] <- 1:12
T_ijt[4,2,1:5] <- 8:12
T_ijt[5,2,1:3] <- 10:12


###  Simulate delay_distributions

delay_dist_sim <- function(df1){
  n.sim <- nrow(df1)
  
  #Could modify one of these to be time from infection to test
  expose.dist= rgamma(n.sim,4,0.75) #duration latent
  infect.dist= rgamma(n.sim,4,0.75) #duration infectiousness
  #end.inf= rgamma(n.sim,4, 0.75) #duration infectiousness
  
  day.infectious <- round(df1$day_index - infect.dist)
  day.exposed <- round(day.infectious - expose.dist)
  day.infectious.end <-  round(day.infectious + infect.dist)
  first.day.hh <- min(day.exposed, na.rm=T) - 1
  last.day.hh <- max(day.infectious.end, na.rm=T)
  
  df1$day.infectious <- day.infectious - first.day.hh
  df1$day.exposed <- day.exposed - first.day.hh
  df1$day.infectious.end <- day.infectious.end - first.day.hh
  df1$max.time.hh <- max(df1$day.infectious.end, na.rm=T)
  
  #res.list <- list('max.time.hh'=max.time.hh,'day.infectious'=day.infectious,'day.exposed'=day.exposed,'day.infectious.end'=day.infectious.end)
  return(df1)
}




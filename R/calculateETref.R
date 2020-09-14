calculateETref <- function(x) {
  
  wthr.df <- as.data.frame(x)
  
  ## Check date and number format in wthr.df and then use the functions
  wthr.df$TS<-ymd_hms(wthr.df$TS) 
  
  wthr.df[ ,2:ncol(wthr.df)] <- apply(wthr.df[ ,2:ncol(wthr.df)], 2, 
                                      function(x) {as.numeric(as.character(x))})
  
  ## extract observations for EACH DATE 
  
  wthr.df$date <- date(wthr.df$TS)
  
  ## extract observations for EACH TIME 
  
  wthr.df$time <- strftime(wthr.df$TS, format="%H:%M:%S", tz="UTC")
  
  # wthr.df %>% group_by(time) %>% mutate(Tmax = summarise(Value = max(wthr.df$T))) 
  
  ## Group by date
  wthr.df1<-wthr.df %>%
    arrange(date, time) %>%
    group_by(date) %>%
    mutate(Tmax = max(Temp), Tmin = min(Temp)) 
  
  
  ## Calculate the Penman-Monteith Equation
  
  # del <- (4098(0.6108*exp(17.27*t/(t+237.3))))/(t+273.3)^2
  # 
  # gammaa <- 0.63189 # P at 500m msl = 95.52 kPa & gamma = i.e.0.000665*95.52  
  # 
  # es = 0.6108*exp(17.27*t/(t + 273.3))
  # 
  # ea <- rh/100*es
  # 
  # Rng <- ((slr.rad*15*60/1000000)*2.5)
  #  
  # ETref <- ((0.408*del*Rng) + 0.063189*((900/96)*ws*(es-ea))/(t+273))/(del+0.063189*(1+0.34*ws))
  
  wthr.df1$ETref <- NA
  
  for(i in 1:nrow(wthr.df1))
  {
    t <- wthr.df1$Temp[i]
    ws <- wthr.df1$WS[i]
    tmax <- wthr.df1$Tmax[i]
    tmin <- wthr.df1$Tmin[i]
    rh <- wthr.df1$RH[i]
    slr.rad <- wthr.df1$SR[i]
    
    del <- (4098*(0.6108*exp(17.27*t/(t+237.3))))/(t+273.3)^2
    
    gammaa <- 0.63189  # P at 500m msl = 95.52 kPa
    
    es = 0.6108*exp(17.27*t/(t + 273.3))
    
    ea <- rh/100*es
    
    Rng <- ((slr.rad*15*60/1000000)*2.5)
    
    ETref <- ((0.408*del*Rng) + 0.063189*((900/96)*ws*(es-ea))/(t+273))/(del+0.063189*(1+0.34*ws)) ## 96 is the constant
    
    if(ETref == 0.0){
      
      wthr.df1$ETref[i] <- NA
      
    } else {wthr.df1$ETref[i] <- ETref}
  }
  
  wthr.df1$ETref <- na.approx(wthr.df1$ETref)
  
  wthr.df1$ETref<-round(wthr.df1$ETref, 3)
  
  return(wthr.df1)
}
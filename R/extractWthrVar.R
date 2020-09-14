extractWthrVar <- function(x, y) {
  
  raw.df <- y; id <- x
  
  wthr.tmp <- raw.df[raw.df$variable==unq.clm.var[id], ]
  wthr.tmp$ts <- ymd_hms(wthr.tmp$timestamp); 
  wthr.tmp <- wthr.tmp[order(wthr.tmp$ts),]
  wthr.tmp <- wthr.tmp[!duplicated(wthr.tmp[c(1,5)]),]  # 'sensor', 'ts'
  
  sen.cnt <- unique(wthr.tmp$sensor)
  
  
  if (length(sen.cnt) > 2)
  {
    
    print(paste0(unq.clm.var[id], " has data from > 2 Sensors"))
    
    # Calculate number of elements per sensor and select the top two with the highest number of values
    unq.SENSOR.val.ln <- c() # create vector to store the number of unique Values per sector
    unq.SENSOR.valSKW <- c() # create vector to store the kurtosis of values per sector
    
    # store #values per sensor values
    for (i in 1:length(sen.cnt))
    {
      t1 <- wthr.tmp[wthr.tmp$sensor == unique(wthr.tmp$sensor)[i], ]
      unq.SENSOR.val.ln <- c(unq.SENSOR.val.ln, nrow(t1))
      unq.SENSOR.valSKW <- c(unq.SENSOR.valSKW, skewness(t1$value))
    }
    
    # find mode of the length vector: unq.sec.dates.ln and remove sectors less than mode
    SENSOR.valLen <- data.frame(SENSOR = sen.cnt, valLEN = unq.SENSOR.val.ln, valSKW = unq.SENSOR.valSKW)
    SENSOR.valLen <- SENSOR.valLen[order(SENSOR.valLen$valLEN, decreasing = TRUE), ]
    
    SENSOR.valLen <- SENSOR.valLen[SENSOR.valLen$valSKW > -1 & SENSOR.valLen$valSKW < 1, ]
    
    wthr.tmp1 <- wthr.tmp[wthr.tmp$sensor == SENSOR.valLen$SENSOR[1], ]
    wthr.var.DF <- wthr.tmp1
    
  } else if (length(sen.cnt) == 2) {
    
    print(paste0(unq.clm.var[id], " has data from 2 sensors"))
    
    sen1 <- wthr.tmp[wthr.tmp$sensor==sen.cnt[1], ]
    
    sen2 <- wthr.tmp[wthr.tmp$sensor==sen.cnt[2], ]
    
    L1 <- nrow(sen1); L2 <- nrow(sen2)
    
    if(L1 > 4*L2 )
      
    { wthr.var.DF <- sen1 } else if (L2 > 4*L1) { 
      wthr.var.DF <- sen2 } else {
        
        stat <- quantile(wthr.tmp$value)
        stat1 <- quantile(sen1$value); stat2 <- quantile(sen2$value); 
        
        dist.ord <- c(dist(rbind(stat, stat1)), dist(rbind(stat, stat2)))
        
        sen.f <- which.min(dist.ord)
        
        wthr.var.DF <- wthr.tmp[wthr.tmp$sensor==sen.cnt[sen.f], ]}
    } else {
      
    print(paste0(unq.clm.var[id], " has data from only 1 sensor"))
    
    wthr.var.DF <- wthr.tmp[wthr.tmp$sensor==sen.cnt[1], ]
  }
  
  return(wthr.var.DF)
  
}
getFeatures <- function(x) {
  
  a <- Sys.time()
  
  spl.smth <- x
  
  # TS and load-cell cols of raw data
  dd.ts <- spl.smth[9:nrow(spl.smth), 6:ncol(spl.smth)] 
  rownames(dd.ts) <- spl.smth$old_unit[9:nrow(spl.smth)]
  dd.ts <- as.data.frame(t(dd.ts))
  dd.ts$TS <- ymd_hms(rownames(dd.ts))
  #dd.ts$TS <- ymd_hms(sub("X","",rownames(dd.ts)))
  dd.ts <- dd.ts[ ,c(ncol(dd.ts), 1:(ncol(dd.ts)-1))]
  rownames(dd.ts) <- seq(1, nrow(dd.ts))
  
  # TS, weather and ETref cols of raw data
  dd.meta <- spl.smth[1:8, 6:ncol(spl.smth)] 
  rownames(dd.meta) <- rownames(spl.smth)[1:8]
  dd.meta <- as.data.frame(t(dd.meta))
  dd.meta$TS <- dd.ts$TS
  dd.meta$date <- date(dd.meta$TS)
  dd.meta$time <- strftime(dd.meta$TS, format="%H:%M:%S", tz="UTC")
  dd.meta <- dd.meta[ ,c((ncol(dd.meta)-2), (ncol(dd.meta)-1), ncol(dd.meta), 
                         1:(ncol(dd.meta)-3))]
  rownames(dd.meta) <- seq(1, nrow(dd.meta))
  
  
  # dd.xts <- xts(dd.ts.f[,-1], order.by = as.POSIXct(dd.ts.f$TS, format = "%Y-%m-%d %H:%M")) ## 1392 entries
  
  
  ## Function: Standardize between 0 to 1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  ## Create Weather_Variables Dataframe
  weather.DAT <- data.frame(matrix(NA, nrow = length(unique(dd.meta$date)), ncol = 15))
  colnames(weather.DAT) <- c("p.avg_RH","p.avg_T","p.avg_VPD","p.avg_Slr.Rad","p.avg_WS", 
                             "mdn.avg_RH","mdn.avg_T","mdn.avg_VPD","mdn.avg_Slr.Rad","mdn.avg_WS", 
                             "evn.avg_RH","evn.avg_T","evn.avg_VPD","evn.avg_Slr.Rad","evn.avg_WS")
  
  ## Create Input Data
  
  spl.smth.df <- dd.ts[,-1]
  
  # Features to be extracted: "maxET", "curv.maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45",
  # "auc.10-15", "sd.10-15", "auc.07-19", "sd.07-19", "auc.prop.10-15", "auc.prop.07-19",
  # "auc.night","total.auc","cos.sim.index"  
  
  ## Create Feature Set list
  F.mat <- matrix(NA, nrow = length(unique(dd.meta$date)), ncol = 15)
  
  colnames(F.mat) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                       "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                       "auc.prop.07-19", "auc.night", "cos.sim.index")
  
  F.mat.list <- replicate(ncol(spl.smth.df), F.mat, simplify=F) # Make a list for each of the geno
  
  dates <- unique(dd.meta$date)
  
  ## Main Body
  for(j in 1:ncol(spl.smth.df)){
    
    spl.smth.whtr.sec.df <- as.data.frame(cbind(dd.meta, spl.smth.df[,j]))
    
    colnames(spl.smth.whtr.sec.df)[12] <- c("Data")
    
    ## extract observations for EACH DATE 
    spl.smth.whtr.sec.df$date <- date(spl.smth.whtr.sec.df$TS)
    
    ## extract observations for EACH TIME 
    spl.smth.whtr.sec.df$time <- strftime(spl.smth.whtr.sec.df$TS, format="%H:%M:%S", tz="UTC")
    
    # wthr.df %>% group_by(time) %>% mutate(Tmax = summarise(Value = max(wthr.df$T))) 
    
    ## Group by date
    spl.smth.whtr.sec.df1 <- spl.smth.whtr.sec.df %>%
      arrange(date, time) %>%
      group_by(date) %>%
      mutate(Dmax = max(Data), Dmin = min(Data)) 
    
    if(!(max(spl.smth.whtr.sec.df1$Data)==0)) # condition to prevent error columns with all 0 values
    {
      for(i in 1:length(dates)){
        
        ## Create tmp data only for dates[i]
        data.tmp <- as.data.frame(spl.smth.whtr.sec.df1[spl.smth.whtr.sec.df1$date == as.character(dates[i]),])
        
        ## Calculate- "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET" 
        #1: maxET
        maxET <- unique(data.tmp$Dmax) 
        if (maxET > 0.0) {etmax.loc <- min(which.max(data.tmp$Data))
        F.mat.list[[j]][i,1] <- maxET} else {F.mat.list[[j]][i,1] <- 0.0}
        print("Feature #1: over")
        
        #2: slope.maxET-6
        if (maxET > 0.0 & etmax.loc > 5){d.tmp<-data.frame(x = 1:6, y = data.tmp$Data[(etmax.loc-5):etmax.loc])
        mod <- lm(y~x, d.tmp, na.action=na.exclude)
        slope.maxET.6 <- mod$coefficients[2] 
        F.mat.list[[j]][i,2] <- slope.maxET.6} else {F.mat.list[[j]][i,2] <- 0.0}
        print("Feature #2: over")
        
        #3: slope.07-maxET
        et.07.loc <- which(hms(data.tmp$time) == hms("07:00:00"))
        if (maxET > 0.0 & etmax.loc > et.07.loc) {d.tmp1<-data.frame(x1 = 1:(etmax.loc - (et.07.loc-1)), y1 = data.tmp$Data[et.07.loc:etmax.loc])
        mod1 <- lm(y1~x1, d.tmp1, na.action=na.exclude)
        slope.07.maxET <- mod1$coefficients[2] 
        F.mat.list[[j]][i,3] <- slope.07.maxET} else {F.mat.list[[j]][i,3] <- 0.0}
        print("Feature #3: over")
        
        #4: slope.00-07
        et.00.loc <- which(hms(data.tmp$time) == hms("00:00:00"))
        d.tmp2<-data.frame(x2 = 1:((et.07.loc+1) - et.00.loc), y2 = data.tmp$Data[et.00.loc:et.07.loc])
        mod2 <- lm(y2~x2, d.tmp2, na.action=na.exclude)
        slope.00.07 <- mod2$coefficients[2]
        F.mat.list[[j]][i,4] <- slope.00.07
        print("Feature #4: over")
        
        #5: slope.19-23:45
        et.19.loc <- which(hms(data.tmp$time) == hms("19:00:00"))
        et.23_45.loc <- which(hms(data.tmp$time) == hms("23:45:00"))
        d.tmp3<-data.frame(x3 = 1:(et.23_45.loc - (et.19.loc-1)), y3 = data.tmp$Data[et.19.loc:et.23_45.loc])
        mod3 <- lm(y3~x3, d.tmp3, na.action=na.exclude)
        slope.19.2345 <- mod3$coefficients[2]
        F.mat.list[[j]][i,5] <- slope.19.2345
        print("Feature #5: over")
        
        #6: "curvmaxET"
        if (length(unique(data.tmp$Data)) > 1) {tsf <- tsfeatures(data.tmp$Data)
        F.mat.list[[j]][i,6] <- tsf$curvature} else {F.mat.list[[j]][i,6] <- 0.0}
        print("Feature #6: over")
        
        ## Calculate- "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
        # "auc.prop.07-19", "auc.night", "cos.sim.index"
        
        #7: total.auc 
        tot.auc <-sum(data.tmp$Data[1]/2, sum(data.tmp$Data[ 2 : (nrow(data.tmp) - 1)]),
                      data.tmp$Data[nrow(data.tmp)]/2)
        F.mat.list[[j]][i,7] <- tot.auc
        print("Feature #7: over")
        
        #8: auc.10-15 
        y1.10 <- which(hms(data.tmp$time) == hms("10:00:00"))
        y1.15 <- which(hms(data.tmp$time) == hms("15:00:00"))
        auc.10.15 <- sum(data.tmp$Data[y1.10]/2, sum(data.tmp$Data[(y1.10 + 1) : (y1.15 - 1)]),
                         data.tmp$Data[y1.15]/2)
        F.mat.list[[j]][i,8] <- auc.10.15
        print("Feature #8: over")
        
        #9: sd.10-15 
        sd.10.15 <- sd(data.tmp$Data[(y1.10) : (y1.15)])
        F.mat.list[[j]][i,9] <- sd.10.15
        print("Feature #9: over")
        
        #10: auc.prop.10-15
        auc.prop.10.15 <- (auc.10.15/tot.auc)
        F.mat.list[[j]][i,10] <- auc.prop.10.15
        print("Feature #10: over")
        
        #11: auc.07-19
        y1.07 <- which(hms(data.tmp$time) == hms("07:00:00"))
        y1.19 <- which(hms(data.tmp$time) == hms("19:00:00"))
        auc.07.19 <- sum(data.tmp$Data[y1.07]/2, sum(data.tmp$Data[(y1.07 + 1) : (y1.19 - 1)]),
                         data.tmp$Data[y1.19]/2)
        F.mat.list[[j]][i,11] <- auc.07.19
        print("Feature #11: over")
        
        #12: sd.07-19
        sd.07.19 <- sd(data.tmp$Data[(y1.07) : (y1.19)])
        F.mat.list[[j]][i,12] <- sd.07.19
        print("Feature #12: over")
        
        #13: auc.prop.07-19 
        auc.prop.07.19 <- (auc.07.19/tot.auc)
        F.mat.list[[j]][i,13] <- auc.prop.07.19
        print("Feature #13: over")
        
        #14: auc.night 
        auc.night <- tot.auc-auc.07.19
        F.mat.list[[j]][i,14] <- auc.night
        print("Feature #14: over")
        
        #15: cos.sim.index
        A = data.tmp$ETref
        B = data.tmp$Data
        cos.sim.index <- sum(A*B)/sqrt(sum(A^2)*sum(B^2)) # cosine similarity function
        F.mat.list[[j]][i,15] <- cos.sim.index
        print("Feature #15: over")
        
        
        ## Get the average of weather parameters
        y1.w.loc <-  which(hms(data.tmp$time) == hms("10:00:00"))
        y2.w.loc <-  which(hms(data.tmp$time) == hms("15:00:00"))
        weather.DAT[i,1] <- p.avg_RH <- mean(data.tmp$RH[y1.w.loc : y2.w.loc], na.rm = TRUE)
        weather.DAT[i,2] <- p.avg_Temp <- mean(data.tmp$Temp[y1.w.loc : y2.w.loc], na.rm = TRUE)
        weather.DAT[i,3] <- p.avg_VPD <- mean(data.tmp$VPD[y1.w.loc : y2.w.loc], na.rm = TRUE)
        weather.DAT[i,4] <- p.avg_SR <- mean(data.tmp$SR[y1.w.loc : y2.w.loc], na.rm = TRUE)
        weather.DAT[i,5] <- p.avg_WS <- mean(data.tmp$WS[y1.w.loc : y2.w.loc], na.rm = TRUE)
        
        y3.w.loc <-  which(hms(data.tmp$time) == hms("00:00:00"))
        y4.w.loc <-  which(hms(data.tmp$time) == hms("06:00:00"))
        weather.DAT[i,6] <- mdn.avg_RH <- mean(data.tmp$RH[y3.w.loc : y4.w.loc], na.rm = TRUE)
        weather.DAT[i,7] <- mdn.avg_Temp <- mean(data.tmp$Temp[y3.w.loc : y4.w.loc], na.rm = TRUE)
        weather.DAT[i,8] <- mdn.avg_VPD <- mean(data.tmp$VPD[y3.w.loc : y4.w.loc], na.rm = TRUE)
        weather.DAT[i,9] <- mdn.avg_SR <- mean(data.tmp$SR[y3.w.loc : y4.w.loc], na.rm = TRUE)
        weather.DAT[i,10] <- mdn.avg_WS <- mean(data.tmp$WS[y3.w.loc : y4.w.loc], na.rm = TRUE)
        
        y5.w.loc <-  which(hms(data.tmp$time) == hms("19:45:00"))
        y6.w.loc <-  which(hms(data.tmp$time) == hms("23:45:00"))
        weather.DAT[i,11] <- evn.avg_RH <- mean(data.tmp$RH[y5.w.loc : y6.w.loc], na.rm = TRUE)
        weather.DAT[i,12] <- evn.avg_Temp <- mean(data.tmp$Temp[y5.w.loc : y6.w.loc], na.rm = TRUE)
        weather.DAT[i,13] <- evn.avg_VPD <- mean(data.tmp$VPD[y5.w.loc : y6.w.loc], na.rm = TRUE)
        weather.DAT[i,14] <- evn.avg_SR <- mean(data.tmp$SR[y5.w.loc : y6.w.loc], na.rm = TRUE)
        weather.DAT[i,15] <- evn.avg_WS <- mean(data.tmp$WS[y5.w.loc : y6.w.loc], na.rm = TRUE)
        
        print(i)
      }}
    
    # allFeatures <- as.data.frame(F.mat.list[[j]])
    
    print(paste0("Finished all features of col_index=", j))
    Sys.sleep(0.2)
  }
  
  b <- Sys.time()
  
  print(paste0("Time taken for feature extraction: ", round((b-a), 2), " minutes" ))
  
  list(allFeatures = F.mat.list, weather.DAT = weather.DAT)
  
}
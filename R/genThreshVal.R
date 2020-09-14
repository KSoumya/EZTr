genThreshVal <- function(x, y) {
  unq.dts <- x; baseTW.df <- y
  threshModelMAT <- data.frame(matrix(nrow = length(unq.dts), ncol=8))
  colnames(threshModelMAT) <- c("Slope", "Intercept", 
                                "Q_0", "Q_25", "Q_50", "Q_75", "Q_100", "solarRAD_Mean")
  
  for(i in 1:length(unq.dts)){
    
    print(paste0("Enter date: ", i, " = ", unq.dts[i]))
    
    tmp <- baseTW.df[baseTW.df$date == unq.dts[i], ]
    
    mod.tmp <- tmp[ , -c(1:4)]
    
    mod.tmp[mod.tmp == 999] <- NA
    
    var.geno <- apply(mod.tmp[,-1], 2, function(x){var(x[!is.na(x)])}) 
    
    max.var <- max(var.geno); 
    range.var <- (median(var.geno)+3*sd(var.geno))
    ifelse(max.var > range.var, mod.tmp <- mod.tmp[ ,-(which.max(var.geno)+1)],
           mod.tmp <- mod.tmp)
    
    tmp.mod.coefMAT <- data.frame(matrix(nrow = (ncol(mod.tmp)-1), ncol = 2))
    colnames(tmp.mod.coefMAT) <- c("Intercept", "Slope")
    
    quantileMAT <- as.data.frame(matrix(nrow = (ncol(mod.tmp)-1) , ncol = 5))
    colnames(quantileMAT) <- c("0%", "25%", "50%", "75%", "100%")
    
    for(j in 2:ncol(mod.tmp)){
      
      print(paste0("Enter geno: ",j, " = ", colnames(mod.tmp)[j], ", when i = ", i))
      
      if(!sum(is.na(mod.tmp[ ,j])) == length(mod.tmp[ ,j]))
      {
        tmp.model <- lm(ETref ~ mod.tmp[ ,j], data = mod.tmp, na.action=na.exclude)       
        
        tmp.mod.coefMAT[j-1, ] <- tmp.model$coefficients
        
      } else {tmp.mod.coefMAT[j-1, ] <- rep(NA, ncol(tmp.mod.coefMAT))}
      
      quantileMAT [j-1, ] <- quantile(mod.tmp[ ,j-1], na.rm = TRUE)
      
      print(paste0("Exit geno: ",j, " = ", colnames(mod.tmp)[j], ", when i = ", i))
    }
    
    threshModelMAT[i, 1]<- median(tmp.mod.coefMAT$Slope, na.rm = TRUE)
    threshModelMAT[i, 2]<- median(tmp.mod.coefMAT$Intercept, na.rm = TRUE)
    threshModelMAT[i, 3]<- median(quantileMAT$`0%`, na.rm = TRUE)
    threshModelMAT[i, 4]<- median(quantileMAT$`25%`, na.rm = TRUE)
    threshModelMAT[i, 5]<- median(quantileMAT$`50%`, na.rm = TRUE)
    threshModelMAT[i, 6]<- median(quantileMAT$`75%`, na.rm = TRUE)
    threshModelMAT[i, 7]<- median(quantileMAT$`100%`, na.rm = TRUE)
    threshModelMAT[i, 8]<- median(tmp$solarRAD, na.rm = TRUE)
    
    print(paste0("Exit date - ", i, " = ", unq.dts[i]))
  }
  
  return(threshModelMAT)
}

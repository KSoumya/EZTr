dataPART <- function(x) {
  
  ETr_smth_IP <- x
  
  by_unit <- ETr_smth_IP[-c(1:8) ,c(1, 6:ncol(ETr_smth_IP))]
  
  geno.ETr <- as.data.frame(t(by_unit))
  
  geno.ETr <- geno.ETr[-1, ]
  
  ETR_baseDF <- as.data.frame(apply(geno.ETr, 2, function(x) {as.numeric(as.character(x))}))
  
  colnames(ETR_baseDF) <- by_unit$unit
  
  rownames(ETR_baseDF) <- colnames(ETr_smth_IP)[6:ncol(ETr_smth_IP)]
  
  # Check date format in rownames and run accordingly
  ETR_baseDF$TS <- ymd_hms(rownames(ETR_baseDF)); 
  ETR_baseDF$date <- date(ETR_baseDF$TS); 
  ETR_baseDF$time <- strftime(ETR_baseDF$TS, format="%H:%M:%S", tz="UTC")
  ETR_baseDF$solarRAD <- as.numeric(as.character(ETr_smth_IP[4, 
                                                             6:ncol(ETr_smth_IP)])) ## Need to convert row to vector using c()
  
  ETR_baseDF <- ETR_baseDF[ ,c((ncol(ETR_baseDF)-3), (ncol(ETR_baseDF)-2), 
                               (ncol(ETR_baseDF)-1), ncol(ETR_baseDF), 
                               1:(ncol(ETR_baseDF)-4))]
  
  
  # Step2: Partitioning ETr data into TWs
  ETr_base_TW1 <- ETR_baseDF[ETR_baseDF$time >= "06:30:00" & ETR_baseDF$time < "18:30:00", ]
  ETr_base_TW2 <- ETR_baseDF[!ETR_baseDF$time %in% ETr_base_TW1$time, ]
  
  list(P1 = ETr_base_TW1[,-c(1:4)], P2 = ETr_base_TW2[,-c(1:4)])
}
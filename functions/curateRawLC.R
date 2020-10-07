curateRawLC <- function(x, y) {
  
  base.d <- x[,-1]
  
  data.RmOut <- base.d 
  
  for ( i in 1:(ncol(base.d)) ) {
    data <- base.d[, i]
    data[data < 0] <- NA # if less than 0, it is wrong.
    data <- as.numeric(unlist(data))
    
    ol <- boxplot(data, plot = FALSE)$out
    tf <- which(data %in% ol)
    
    data.RmOut[tf, i] <- NA
    
    # jpeg(paste0(opPATH, "Plot_BOX.OL/Rplot", colnames((base.d)[i]), ".jpeg"))
    # plot(data, cex = 1, lwd = 0.8, ylim=c(60000,90000))
    # points(x = (1:length(data))[-tf], y = data[-tf], ylim=c(60000,90000),
    #        col = "red", pch = 20, cex = 0.8)
    # dev.off()
  }
  
  rmOut <- as.data.frame(t(data.RmOut))
  
  rmOut.DF <- as.data.frame(cbind(y, rmOut))
  colnames(rmOut.DF) <- colnames(LC.MAT.raw)
  
  # write.csv(rmOut.DF, paste0(opPATH, "LC.MAT.BOX.OL_s1.csv"))
  
  
  ####### Imputation #######
  data.Imp <- data.frame(apply(data.RmOut, 2, na.approx, na.rm = F))
  
  data.Imp.df <- as.data.frame(t(data.Imp))
  
  imputed.DF <- as.data.frame(cbind(y, data.Imp.df))
  colnames(imputed.DF) <- colnames(LC.MAT.raw)
  
  # write.csv(imputed.DF, paste0(opPATH, "LC.MAT.BOX.OL_s1_IMPUTED.csv"))
  
  
  # There still could remain columns with maximum missing, hence run the below script
  na.list <- as.numeric(apply(imputed.DF[,-c(1:5)], 1, FUN = function(x) {sum(is.na(x))}))
  
  # keep sectors which have more than 30% of values
  na.list.G.Locs <- which(na.list > ceiling(0.3*dim(imputed.DF[,-c(1:5)])[2])) 
  
  # replace 'na.list.G.Locs' only with 0 if NA so that ETr can be extracted
  if(length(na.list.G.Locs) > 0)
  {
    imputed.DF.tmp <- imputed.DF
    imputed.DF.tmp <- imputed.DF.tmp[-(na.list.G.Locs+5), ]
  } else {imputed.DF.tmp <- imputed.DF}
  
  imputed.DF.final <- imputed.DF.tmp
  
  interp.ip <- imputed.DF.final[ ,6:ncol(imputed.DF.final)]
  
  interp.df <- as.data.frame(t(interp.ip))
  
  interp.df.op <- as.data.frame(apply(interp.df, 1, na.aggregate.default))
  
  imputed.DF.final[ ,6:ncol(imputed.DF.final)] <- interp.df.op
  
  return(imputed.DF.final)
}
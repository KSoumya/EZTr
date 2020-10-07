generateETrRatio <- function(x) {
  
  wthr.ETref.ETobs.Ratio <- as.data.frame(t(x[8:nrow(x), 6:ncol(x)]))
  
  for(j in 2:ncol(wthr.ETref.ETobs.Ratio))
  {
    for(i in 1:nrow(wthr.ETref.ETobs.Ratio))
    {
      if (is.na(wthr.ETref.ETobs.Ratio$ETref[i]) || wthr.ETref.ETobs.Ratio$ETref[i] == 0.0)
        wthr.ETref.ETobs.Ratio[i, j] <- 999
      else
        wthr.ETref.ETobs.Ratio[i, j] <- wthr.ETref.ETobs.Ratio[i, j]/wthr.ETref.ETobs.Ratio$ETref[i]
    }
    
    wthr.ETref.ETobs.Ratio[, j] <- na.approx(wthr.ETref.ETobs.Ratio[, j])
  }
  
  
  ET_ratio_mat <- x
  
  ET_ratio_mat[9:nrow(ET_ratio_mat), 6:ncol(ET_ratio_mat) ] <- round(t(wthr.ETref.ETobs.Ratio[ ,-1]), 3)
  
  return(ET_ratio_mat)
  
}
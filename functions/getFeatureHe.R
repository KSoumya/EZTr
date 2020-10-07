getFeatureHe <- function(x, y, d, p) {
  
  allFeatures <- x; raw.trans <- y
  
  opPATH <- p; dates <- d
  
  for (i in 1:length(dates)){
    
    # create tmp DF to store feature values of all genos of 'each date'
    tmp <- as.data.frame(matrix(NA, nr = (nrow(raw.trans)-8), ncol = 15)) # 15 corresponds to the features
    colnames(tmp) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                       "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                       "auc.prop.07-19", "auc.night", "cos.sim.index")
    
    for(j in 1:length(allFeatures)){
      
      tmp[j, ] <- allFeatures[[j]][i,]}
    
    # tmp$g.ID <- colnames(spl.smth.df)
    
    print(i)
    
    features.DF <- as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], tmp))
    
    write.csv(features.DF, paste0(opPATH, "D.", dates[i],".csv"))
    
    
    ## Calculate heritability estimates of all the features for each date ##
    # Heritability analysis for each feature, for each day, as well as on average of feature values
    he.vec <- c()
    
    for(f in 1:(ncol(features.DF)-5))
      
    {
      he.tmp <- features.DF[ , c(3, (f+5))] # 3 is Genotype-column and f+5 is the feature column
      he.tmp$Genotype <- factor(he.tmp$Genotype)
      
      y <- he.tmp[,2]
      selector <- !is.na(y)
      y <- y[selector]
      line <- as.factor(as.character(he.tmp$Genotype[selector]))
      
      # model for He
      model <- tryCatch(gam(y ~ s(line, bs = "re"), method = "REML"),
                error = function(e) NULL)
      
      if(!is.null(model)){
        res <- gam.vcomp(model)
        vg <- res[[1]]^2
        ve <- res[[4]]^2
        he <- vg / (vg + ve)}
        
        # he.vec[f] <- round(he, 3)
        he.vec <- c(he.vec, round(he, 3))
      
      }
    
    
    F.He[i, ] <- he.vec
    
  }
  
  return(F.He)
}
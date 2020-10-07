
prepcsWthr <- function(x, y) {
  
  
  temp.mat <- matrix(x[ ,y], nrow = (24*60), 
                     ncol = nrow(x)/(24*60), byrow = FALSE)
  
  # remove outliers from each date
  for(i in 1:ncol(temp.mat))
  {
    ol <- boxplot(temp.mat[,i], plot = FALSE)$out
    
    temp.mat[temp.mat[,i] %in% ol, i] <- NA
  }
  
  OP.vec <- as.vector(as.matrix(temp.mat))
  
  OP.vec <- na.aggregate.default(OP.vec)
  
  return(OP.vec)
}
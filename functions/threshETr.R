threshETr <- function(x, y, z) {
  
  thresh.val <- x; ratio.df <- y; etr.df <- z
  
  t.index <- c()
  
  ETr.tw.tmp <- etr.df
  
  for(c in 1:ncol(etr.df)){
    
    t.index <- which(abs(ratio.df[ ,c]) >= thresh.val)
    ETr.tw.tmp[t.index, c] <- NA
    
  }
  
  return(ETr.tw.tmp)
}
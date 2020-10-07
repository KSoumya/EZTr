filterLCExtremeCols <- function(x, y) {
  
  des.stats.ip <- as.data.frame(t(x[,-c(1:5)]))
  colnames(des.stats.ip) <- x$unit
  
  des.stats <- sapply(des.stats.ip, quantile)
  
  q0 <- which(des.stats[1, ] %in% boxplot(des.stats[1, ], plot = FALSE)$out)
  q25 <- which(des.stats[2, ] %in% boxplot(des.stats[2, ], plot = FALSE)$out)
  q50 <- which(des.stats[3, ] %in% boxplot(des.stats[3, ], plot = FALSE)$out)
  q75 <- which(des.stats[4, ] %in% boxplot(des.stats[4, ], plot = FALSE)$out)
  q100 <- which(des.stats[5, ] %in% boxplot(des.stats[5, ], plot = FALSE)$out)
  
  err.cols <- intersect(intersect(intersect(intersect(q25, q50), q75), q100), q0)
  
  err.sec.nm <- colnames(des.stats)[err.cols]
  
  err.sec.meta <- y[y$unit %in% err.sec.nm, ]
  
  list(err.sec.META = err.sec.meta, err.sec.NM = err.sec.nm)
}
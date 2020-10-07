filterETrExtremeCols <- function(x, y) {
  des.stats2.ip <- as.data.frame(t(x[ ,6:ncol(x)]))
  colnames(des.stats2.ip) <- x$unit
  
  des.stats2 <- as.data.frame(sapply(des.stats2.ip, quantile, na.rm = TRUE))
  
  des.stats2 <- as.data.frame(apply(des.stats2, 1, function(x) {as.numeric(as.character(x))}))
  rownames(des.stats2) <- x$unit
  
  q0 <- which(des.stats2[, 1] %in% boxplot.stats(des.stats2[, 1])$out)
  q25 <- which(des.stats2[, 2] %in% boxplot.stats(des.stats2[, 2])$out)
  q50 <- which(des.stats2[, 3] %in% boxplot.stats(des.stats2[ ,3])$out)
  q75 <- which(des.stats2[, 4] %in% boxplot.stats(des.stats2[ ,4])$out)
  q100 <- which(des.stats2[, 5] %in% boxplot.stats(des.stats2[ ,5])$out)
  
  err.cols <- intersect(intersect(intersect(intersect(q25, q50), q75), q100), q0)
  
  err.sec.nm <- rownames(des.stats2)[err.cols]
  
  err.sec.meta <- y[y$unit %in% err.sec.nm, ]
  
  list(ETr_err.sec.META = err.sec.meta, ETr_err.sec.NM = err.sec.nm)
}
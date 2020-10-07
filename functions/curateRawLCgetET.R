
curateRawLCgetET <- function(x, y) {
  
  x = LC.MAT.f; y = meta.LCDF
  
  # ifelse(str_detect(y$old_unit, ":"), 
  #        y$old_unit <- sub(":", "_", y$old_unit), y$old_unit<-y$old_unit)
  
  base.d <- x[,-1]
  
  data.RmOut <- base.d 


# stat.mat <- as.data.frame(apply(base.d, 2, function(x) {quantile(x, na.rm = TRUE)}))
  allStats <- as.data.frame(matrix(NA, nrow = 11, ncol = ncol(base.d)))
  rownames(allStats) <- c("mean", "sd", "median", "trimmed", "mad", "min", 
                          "max", "range", "skew", "kurtosis", "se")
  colnames(allStats) <- colnames(base.d)
  
  for(i in 1:ncol(base.d))
  {
    tmp <- c(describe(base.d[ ,i])[-c(1:2)])
    allStats[,i] <- c(unlist(psych::describe(base.d[ ,i])[-c(1:2)]))
  }
  
  
  # Function to get lower and upper limits
  maxL <- function(x) {(mean(as.numeric(x))+1.5*sd(as.numeric(x)))}
  minL <- function(x) {(mean(as.numeric(x))-1.5*sd(as.numeric(x)))}
  
  # Sector-IDs with values beyond limits
  min.OL0<- which(allStats[6, ] < minL(allStats[6, ]))
  max.OL100 <- which(allStats[7, ] > maxL(allStats[7, ]))
  max.OLskw <- which(allStats[9, ] > maxL(allStats[9, ]))
  min.OLskw <- which(allStats[9, ] < minL(allStats[9, ]))
  max.OLkurt <- which(allStats[10, ] > maxL(allStats[10, ]))
  min.OLkurt <- which(allStats[10, ] < minL(allStats[10, ]))
  max.OLse <- which(allStats[11, ] > maxL(allStats[11, ]))
  min.OLse <- which(allStats[11, ] < minL(allStats[11, ]))
  
  olSECs <- unique(c(min.OL0, max.OL100, max.OLskw, min.OLskw, 
                     max.OLkurt, min.OLkurt, max.OLse, min.OLse))
  # olSECnms <- as.character(LC.MAT.raw$old_unit[olSECs])

# Remove outliers from data  
base.d <- base.d[ ,-olSECs]
data.RmOut <- data.RmOut[ ,-olSECs]
y.tmp <- y[y$unit %in% colnames(base.d), ]

for ( i in 1:ncol(base.d)) {
  data <- as.numeric(c(base.d[, i]))
  data[data < 0] <- NA # if less than 0, it is wrong.
  res.dwt <- modwt(data, filter = "haar", n.level = 3, boundary = "periodic", fast = FALSE)
  bp.s1 <- boxplot.stats(res.dwt@W$W1); tf1 <- res.dwt@W$W1 %in% bp.s1$out#; tf1 <- MyFun(tf1)
  bp.s2 <- boxplot.stats(res.dwt@W$W2); tf2 <- res.dwt@W$W2 %in% bp.s2$out#; tf2 <- MyFun(tf2)
  bp.s3 <- boxplot.stats(res.dwt@W$W3); tf3 <- res.dwt@W$W3 %in% bp.s3$out#; tf3 <- MyFun(tf3)
  tf <- tf1 | tf2 | tf3
  data.RmOut[tf, i] <- NA

  # jpeg(paste0(opPATH, "Plot_dwt.OL/", y.tmp$old_unit[i], ".jpeg"))
  # plot(data, cex = 1, lwd = 0.8)
  # points(x = (1:length(data))[!tf], y = data[!tf], col = "red", pch = 20, cex = 0.8)
  # dev.off()
}
# write.csv(data.RmOut, file = "./RE-ANALYSIS_29.06/S1_2 results/allSECs.dwt.OL_s1.csv")


# Remove outlier again
data.RmOut.2 <- as.data.frame(data.RmOut)
for (i in 1:ncol(base.d) ) {
  data <- data.RmOut[, i]
  # plot(data)
  # points(data.im, col = "red", pch = 20, cex = .8)
  data.im <- na.locf(data, na.rm = FALSE)
  DIF <- diff(data.im)
  Threshold <- quantile(DIF, prob = 1 - (4 / length(DIF)), na.rm = T)
  #sum(DIF > Threshold, na.rm = T)
  tf <- DIF > Threshold
  tf[is.na(tf)] <- FALSE # NA should be FALSE
  tf <- c(FALSE, tf) # NA should be FALSE
  num <- (1:length(tf))[tf]
  num.rm <- as.numeric(c(sapply(num, "+", 0:12)))
  data.RmOut.2[num.rm, i] <- NA

  # jpeg(paste0(opPATH, "Plot_EMPRCL.OL/", y.tmp$old_unit[i], ".jpeg"))
  # plot(data)
  # points(num.rm, data[num.rm], col = "red", pch = 20, cex = 0.9)
  # dev.off()

  Sys.sleep(0.3)
}

ext <- which(!rownames(data.RmOut.2) %in% rownames(data.RmOut))
if (length(ext)>0)
  {data.RmOut.2 <- data.RmOut.2[-c(ext),]}else
    {data.RmOut.2 <- data.RmOut.2}

# ##### Just for revision #####
# data.RmOut.2 <- data.RmOut
# ##########################

# Imputation
data.Imp <- data.frame(apply(data.RmOut.2, 2, na.approx, na.rm = F))
rownames(data.Imp) <- rownames(data.RmOut.2)
data.final <- as.data.frame(apply(data.Imp, 2, na.aggregate))
rownames(data.final) <- rownames(data.RmOut.2)

# ET values
data.et <- data.final
data.et[1, ]<-0

for(i in 1:ncol(data.final))
{
  data.et[-1,i] <- -diff(data.final[, i])
  
  # jpeg(paste0(opPATH, "Plot_ET.final/", y.tmp$old_unit[i], ".jpeg"))
  # plot(-diff(data.final[, i]), type = "p")
  # # points(-diff(data.final[, i]), type = "l")
  # dev.off() 
}

rownames(data.et) <- x$TS
etMeta <- as.data.frame(cbind(y.tmp, t(data.et)))

return(etMeta)

}

getETr <- function(x) {
  data.final.err.rmvd <- as.data.frame(t(x))
  
  # ET values
  data.et <- as.data.frame(apply(data.final.err.rmvd[-c(1:5), ], 2, 
                                 function(x) {as.numeric(as.character(x))}))
  
  data.et.tmp <- data.et
  
  data.et[1, ] <- rep(0.0, ncol(data.et))
  
  for(i in 1:ncol(data.et))
  {
    data.et[-1,i] <- -diff(data.et.tmp[, i])
  }
  
  # Remove outliers from data.et
  data.RmOut2 <- data.et
  
  et.obs <- as.data.frame(t(data.et))
  
  ETr_Meta <- impData.errSEC.rmvd
  
  ETr_Meta[ ,6:ncol(ETr_Meta)] <- et.obs
  
  list(obsETr_core = et.obs, obsETr_meta = ETr_Meta)
}
ordFiltETr <- function(x, y) {
  
  ETr_filt_TW1 <- x; ETr_smth_IP <- y
  
  ETr_filt_TW1_ord <- ETr_filt_TW1[,order(match(colnames(ETr_filt_TW1), 
                                                colnames(ETr_smth_IP)[6:ncol(ETr_smth_IP)]))]
  ETr_filt_TW1_ord$TS <- rownames(ETr_filt_TW1_ord)
  ETr_filt_TW1_ord$TS <- ymd_hms(ETr_filt_TW1_ord$TS)
  ETr_filt_TW1_ord <- ETr_filt_TW1_ord[ ,c(ncol(ETr_filt_TW1_ord), 1:(ncol(ETr_filt_TW1_ord)-1))]
  
  return(ETr_filt_TW1_ord)
  
}
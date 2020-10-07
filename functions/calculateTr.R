	calculateTr <- function(x, y, z, d) {
  
  ETr_filt_imputed_FILE <- x; pe.df.ETr <- y;
  
  LAI.mat <- z; unq.dts <- d
  
  ETr_smth.mat <- ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), 
                                        6:ncol(ETr_filt_imputed_FILE)]
  Trans.mat <- ETr_smth.mat #Empty Trans Matrix
  
  TR.mat <- Trans.mat #Empty Trans-Rate Matrix
  
  for(i in 1:nrow(LAI.mat)){
    
    la.df.tmp <- pe.df.ETr[pe.df.ETr$old_unit == rownames(LAI.mat)[i], ]
    
    la.df.med <- la.df.tmp %>% group_by(date) %>% 
      dplyr::summarise(Med_LeafArea3D = median(LeafArea3D))
    
    date.mat <- data.frame(date = unq.dts, val = rep(NA, length(unq.dts)))
    
    for(j in 1:nrow(date.mat))
    {
      ifelse(date.mat$date[j] %in% la.df.med$date, 
             date.mat$val[j] <- la.df.med$Med_LeafArea3D[la.df.med$date == date.mat$date[j]],
             date.mat$val[j] <- NA)  
    }
    
    date.mat$val <- na.aggregate.default(date.mat$val)
    # date.mat$val <- na.spline(date.mat$val)
    
    LAI.all.dates[i, ] <- date.mat$val
    
    sec.lai.tmp <- rep(date.mat$val, each = 96)
    
    # Calculate LAI #
    # LAI.mat[i, ] <- ((((sec.lai.tmp/100) - 29.9)/0.36)*(1/0.26)/10000)
    LAI.mat[i, ] <- ((((sec.lai.tmp/100))/0.36)*(1/0.26)/10000)
    # LAI.mat[i, ] <- (((sec.lai.tmp*1.0114) - 0.0842)/0.26)/10^6
    
    # Calculate Transpiration #
    Trans.mat[i, ] <- (1-(1-exp(-0.463*LAI.mat[i, ])))*ETr_smth.mat[i, ] 
    
    # plot(x=1:length(Trans.mat[i, ]), y=Trans.mat[i, ])
    # lines(x=1:length(Trans.mat[i, ]), y=ETr_smth.mat[i, ], 
    #      type = "l", lty = 1, col="red")
    
    # Calculate Transpiration Rate #
    # TR.mat[i, ] <- (Trans.mat[i, ]/ (LAI.mat[i, ]*0.26*10^4))
    # TR.mat[i, ] <- (Trans.mat[i, ]/ (sec.lai.tmp/100))
        
    # plot(x=1:length(TR.mat[i, ]), y=TR.mat[i, ])
    # lines(x=1:length(TR.mat[i, ]), y=ETr_smth.mat[i, ], 
    #       type = "l", lty = 1, col="red")
    
  }
  
  list(Trans.mat = Trans.mat, LA3D_TS = LAI.all.dates, LAI.mat = LAI.mat)
}
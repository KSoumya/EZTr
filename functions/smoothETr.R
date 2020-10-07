smoothETr <- function(x) {
  
  e21op<-x
  
  e21ip<-(e21op)*1000 ## was not working for values ranged between 0-1, hence multiplied with 10^3
  
  e21ip[e21ip < -100]<--1 ## replace extreme -ve values (irrigation drainage) with -1
  
  for(i in 1:nrow(e21ip))
  {
    nlTS<- as.data.frame(nonLinearNoiseReduction(e21ip[i, ], 3, 12))
    names(nlTS)<-"x"
    
    # Use O/P of 1 as I/P to Spline Smoothing
    
    spldf<-smooth.spline(nlTS)
    
    e21op[i, ]<-spldf$y
    
  }
  
  plot(x=1:ncol(e21ip), y=(e21ip[55, ])/1000, col="black", type="p", main = "An example plot")
  lines(x=1:ncol(e21ip), y=(e21op[55, ])/1000, col="red")
  
  return(e21op)
  
}
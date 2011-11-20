plot.blm <- function(x,...){

  X2 <- gof(x)

  if(length(X2$observed)<10){
    
    barplot(rbind(X2$expected,X2$observed),beside=TRUE,
            legend=c("Expected","Observed"),xlab="Covariate Profile",ylab="",...)
  }
  else{
    Y = model.frame(x@formula,x@data)[,1]
    plot.predictive(predict(x),Y,...)
  }

}

plot.lexpit<- function(x,...){

  X2 <- gof(x)

  if(length(X2$observed)<10){
    
    barplot(rbind(X2$expected,X2$observed),beside=TRUE,
            legend=c("Expected","Observed"),xlab="Covariate Profile",ylab="",...)
  }
  else{
    Y = model.frame(x@formula.linear,x@data)[,1]
    plot.predictive(predict(x),Y,...)
  }

}

plot.predictive <- function(predicted,y,...){

  p <- mean(y)
  pre.bar <- mean(predicted)

  plot(y=sort(predicted),x=(1:length(predicted))/length(predicted)*100,las=1,
       ylab="Predicted risk",
       xlab="Risk percentile",type="l",
       ...)

  abline(h=p,col="red")

  deciles <- cut(predicted, breaks = quantile(predicted, prob = seq(0, 
        1, 1/10)), include.low = T)
  
  p.bar <- tapply(predicted,deciles,mean)
  decile.midpoints <- seq(.05,.95,1/10) 

  points(y=p.bar,x=decile.midpoints*100)
  
}

setMethod("plot",signature(x="blm"),plot.blm)
setMethod("plot",signature(x="lexpit"),plot.lexpit)

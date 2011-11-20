hosmerlem = function(y,predicted,groups=10,plot=FALSE,risk.labels=TRUE,sig=3,...) {
  
  group <- cut(predicted,breaks=quantile(predicted,prob=seq(0,1,1/groups)),include.low=T)
  
  O = tapply(y,group,sum)
  E = tapply(predicted,group,sum)

  chisq = sum((O - E)^2/E)

  P = 1 - pchisq(chisq, groups - 2)

  if(plot){

    midpoints <- quantile(predicted,seq(0+(1/(2*groups)),1-(1/(2*groups)),length=groups))

    plot(y=O,x=E,ylab="Observed",xlab="Expected",pch=19)

    if(risk.labels){
      text(y=O,x=E,labels=round(midpoints,sig),...)
   }


    abline(a=0,b=1)
    }
    
  return(list(chisq=chisq,p.value=P,O=O,E=E,group=table(group)))
}

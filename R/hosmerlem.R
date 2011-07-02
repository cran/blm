hosmerlem = function(f,fit,data,groups=10) {

  X <- model.matrix(f,data)
  Y <- model.frame(f,data)[,1]

  predicted <- fitted(fit)
  o <- order(predicted)
  
  group <- sort(rep(1:groups,length=length(predicted)))
  
  Y <- Y[o]
  predicted <- predicted[o]

  obs = xtabs(cbind(1 - Y, Y) ~ group)

  expect = xtabs(cbind(1 - predicted, predicted) ~ group)

  chisq = sum((obs - expect)^2/expect)

  P = 1 - pchisq(chisq, groups - 2)

  return(list(chisq=chisq,p.value=P))
}

setMethod("coef","blm",function(object,...) {object@fit$par})
setMethod("coef","lexpit",function(object,...) {object@fit$par})
